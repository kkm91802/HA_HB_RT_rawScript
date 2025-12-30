#include "EnergyEnvelope.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static Detection energyDet[MAX_DETECTIONS];
static int energyCount = 0;

static void addEnergyDetection(const char *type, const char *timeStr, float amplitude) {
    if (energyCount >= MAX_DETECTIONS) return;
    Detection *d = &energyDet[energyCount++];
    strcpy(d->polarity, type);
    int h, m, s, ms;
    sscanf(timeStr, "%d:%d:%d:%d", &h, &m, &s, &ms);
    d->hour = h; d->min = m; d->sec = s; d->ms = ms;
    d->amplitude = amplitude;
}


// ─────────────────────────────────────────────
// Simple DFT (for clarity) — replace with FFTW3 for performance
// ─────────────────────────────────────────────
static void simpleDFT(const double *in, double *out_real, double *out_imag, int N)
{
    for (int k = 0; k < N; k++) {
        double sum_r = 0.0, sum_i = 0.0;
        for (int n = 0; n < N; n++) {
            double angle = -2.0 * M_PI * k * n / N;
            sum_r += in[n] * cos(angle);
            sum_i += in[n] * sin(angle);
        }
        out_real[k] = sum_r;
        out_imag[k] = sum_i;
    }
}

// ─────────────────────────────────────────────
// CSV Loader
// ─────────────────────────────────────────────
int loadEnvelopeCSV(const char *filename, EnvelopeSample *data)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("open"); return -1; }

    char line[256];
    int count = 0;
    fgets(line, sizeof(line), fp); // Skip header

    while (fgets(line, sizeof(line), fp) && count < MAX_SAMPLES_ENVELOPE) {
        if (sscanf(line, "%[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                   data[count].time,
                   &data[count].Ax, &data[count].Ay, &data[count].Az,
                   &data[count].Gx, &data[count].Gy, &data[count].Gz) == 7)
            count++;
    }

    fclose(fp);
    return count;
}

// ─────────────────────────────────────────────
// Teager Energy Operator
// ─────────────────────────────────────────────
void computeTeagerEnergy(const double *x, double *E, int n)
{
    for (int i = 1; i < n - 1; i++)
        E[i] = x[i] * x[i] - x[i - 1] * x[i + 1];
    E[0] = E[n - 1] = 0.0;
}

// ─────────────────────────────────────────────
// Adaptive FCFR computation
// ─────────────────────────────────────────────
static int computeFCFR(const double *E2, int n, double fs, double *max_freq, double *fcfr_val)
{
    double max_val = 0.0, mean = 0.0;
    int max_idx = 0;

    for (int i = 1; i < n / 2; i++) {
        double val = fabs(E2[i]);
        mean += val;
        if (val > max_val) {
            max_val = val;
            max_idx = i;
        }
    }

    mean /= (n / 2);
    if (mean < EPS) mean = EPS;

    *fcfr_val = max_val / mean;
    *max_freq = (double)max_idx * fs / n;

    // ─────────────────────────────────────────────
    // Adaptive FCFR Threshold (Warm-up then freeze)
    // ─────────────────────────────────────────────
    static double fcfr_mean = 0.0;
    static double fcfr_M2 = 0.0;
    static int fcfr_count = 0;
    static double adaptive_threshold = FCFR_THRESHOLD; // fallback 3.0

    fcfr_count++;
    if (fcfr_count <= 30) { // first 30 windows = baseline phase
        double delta = *fcfr_val - fcfr_mean;
        fcfr_mean += delta / fcfr_count;
        fcfr_M2 += delta * (*fcfr_val - fcfr_mean);
        double fcfr_std = (fcfr_count > 1) ? sqrt(fcfr_M2 / (fcfr_count - 1)) : 0.0;

        adaptive_threshold = fcfr_mean + FCFR_WEIGHT * fcfr_std;
    }

    return (*fcfr_val > adaptive_threshold);
}

// ─────────────────────────────────────────────
// Compute Improved Envelope for one window
// ─────────────────────────────────────────────
void computeImprovedEnvelopeWindow(const EnvelopeSample *data, int startIdx, int n, double fs)
{
    static int lastDetectedIdx = -100000;  // remember last detection
    const int minSeparation = n / 2;       // ignore detections closer than half a window

    double *x = calloc(n, sizeof(double));
    double *E1 = calloc(n, sizeof(double));
    if (!x || !E1) { fprintf(stderr, "Memory error\n"); return; }

    for (int i = 0; i < n; i++)
        x[i] = data[startIdx + i].Ax;

    computeTeagerEnergy(x, E1, n);

    // FFT of E1
    double *real = calloc(FFT_SIZE, sizeof(double));
    double *imag = calloc(FFT_SIZE, sizeof(double));
    double *mag = calloc(FFT_SIZE, sizeof(double));
    for (int i = 0; i < n && i < FFT_SIZE; i++)
        real[i] = E1[i];

    simpleDFT(real, real, imag, FFT_SIZE);
    for (int i = 0; i < FFT_SIZE; i++)
        mag[i] = sqrt(real[i]*real[i] + imag[i]*imag[i]);

    // Second TEO with smoothing
    double *E2 = calloc(FFT_SIZE, sizeof(double));
    for (int i = 1; i < FFT_SIZE - 1; i++) {
        double m_prev = (mag[i - 1] + mag[i]) / 2.0;
        double m_next = (mag[i + 1] + mag[i]) / 2.0;
        E2[i] = mag[i]*mag[i] - m_prev * m_next;
    }

    // Compute FCFR (adaptive)
    double fcfr_val, max_freq;
    int detected = computeFCFR(E2, FFT_SIZE, fs, &max_freq, &fcfr_val);

    if (detected) {
        // find index of maximum in E1
        double maxE = -1e9; int maxIdx = 0;
        for (int i = 0; i < n; i++) {
            if (E1[i] > maxE) { maxE = E1[i]; maxIdx = i; }
        }

        // refine localization by searching original Ax around maxIdx
        int searchRadius = n / 8;
        if (searchRadius < 5) searchRadius = 5;
        int lo = startIdx + maxIdx - searchRadius;
        int hi = startIdx + maxIdx + searchRadius;
        if (lo < 0) lo = 0;
        if (hi >= MAX_SAMPLES_ENVELOPE) hi = MAX_SAMPLES_ENVELOPE - 1;

        double bestVal = data[startIdx + maxIdx].Ax;
        int bestIdx = startIdx + maxIdx;
        for (int idx = lo; idx <= hi; idx++) {
            double val = fabs(data[idx].Ax);
            if (val > fabs(bestVal)) {
                bestVal = data[idx].Ax;
                bestIdx = idx;
            }
        }

        // ─── Debounce + amplitude filter ───
        if (fabs(bestVal) >= 0.5 && (bestIdx - lastDetectedIdx > minSeparation)) {
            printf("Positive spike at %s (%.3f m/s²)\n",
                   data[bestIdx].time, data[bestIdx].Ax);
                       addEnergyDetection("Positive", data[bestIdx].time, data[bestIdx].Ax);
            lastDetectedIdx = bestIdx;
        }
    }

    free(x); free(E1); free(real); free(imag); free(mag); free(E2);
}

// ─────────────────────────────────────────────
// Sliding-Window ITEO Processing (90% overlap)
// ─────────────────────────────────────────────
void runImprovedEnvelopeSliding(const EnvelopeSample *data, int totalSamples, double fs)
{
    printf("\n--- Improved Energy Envelope (Sliding Window Analysis) ---\n");

    int step = WINDOW_SAMPLES / 10; // 90% overlap
    if (step < 1) step = 1;
    int nWindows = 1;
    if (totalSamples > WINDOW_SAMPLES) 
        nWindows = (totalSamples - WINDOW_SAMPLES) / step + 1;
    else
        nWindows = 0;

    for (int w = 0; w < nWindows; w++) {
        int start = w * step;
        computeImprovedEnvelopeWindow(data, start, WINDOW_SAMPLES, fs);
    }

    printf("\nAdaptive FCFR thresholding applied.\n");
    printf("----------------------------------------------------------\n");
}

int getEnergyDetections(Detection **array) {
    *array = energyDet;
    return energyCount;
}