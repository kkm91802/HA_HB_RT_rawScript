#include "DerivativeDetect.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

// ─────────────────────────────────────────────
// CSV Reader
// ─────────────────────────────────────────────
int loadCSV(const char *filename, DerivSample *data)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("open"); return -1; }

    char line[256];
    int count = 0;

    // Skip header if present
    fgets(line, sizeof(line), fp);

    while (fgets(line, sizeof(line), fp) && count < MAX_SAMPLES_DERIV) {
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
// Compute derivative energy and max magnitude per window
// ─────────────────────────────────────────────
static void computeDerivativeEnergy(const DerivSample *data,
                                    double *d_total, double *max_ax, int n)
{
    double a_prev = data[0].Ax;

    for (int k = 0; k < n / DERIV_WINDOW_SIZE; k++) {
        double sum = 0.0;
        double local_max = fabs(data[k * DERIV_WINDOW_SIZE].Ax);

        for (int i = 0; i < DERIV_WINDOW_SIZE; i++) {
            int idx = k * DERIV_WINDOW_SIZE + i;
            if (idx == 0) continue;

            double a_now = data[idx].Ax;
            double deriv = (a_now - a_prev);
            sum += fabs(deriv);

            if (fabs(a_now) > local_max)
                local_max = fabs(a_now);

            a_prev = a_now;
        }
        d_total[k] = sum / DERIV_WINDOW_SIZE;
        max_ax[k] = local_max;
    }
}

// ─────────────────────────────────────────────
// Compute mean and standard deviation
// ─────────────────────────────────────────────
static void computeStats(const double *x, int n, double *mean, double *std)
{
    double mu = 0.0, sigma = 0.0;
    for (int i = 0; i < n; i++) mu += x[i];
    mu /= n;
    for (int i = 0; i < n; i++) sigma += (x[i] - mu) * (x[i] - mu);
    sigma = sqrt(sigma / n);
    *mean = mu; *std = sigma;
}

// ─────────────────────────────────────────────
// Main detection pipeline (hysteresis + ±2-window search)
// ─────────────────────────────────────────────
void detectSpikes(const DerivSample *data, int n)
{
    int nWindows = n / DERIV_WINDOW_SIZE;
    double *d_total = calloc(nWindows, sizeof(double));
    double *max_ax = calloc(nWindows, sizeof(double));
    if (!d_total || !max_ax) return;

    computeDerivativeEnergy(data, d_total, max_ax, n);

    // Adaptive thresholds and hysteresis
    double mu, sigma;
    computeStats(d_total, nWindows, &mu, &sigma);

    double pos_enter = mu + 2.5 * sigma;
    double pos_exit  = mu + 0.5 * sigma;
    double neg_enter = -(mu + 2.5 * sigma);
    double neg_exit  = -(mu + 0.5 * sigma);

    printf("Processing %d samples...\n", n);

    int pos_active = 0, neg_active = 0;
    int start_idx = 0;
    double peak_val = 0.0;
    int peak_idx = 0;

    for (int k = 0; k < nWindows; k++) {

        // ───── Positive region detection with hysteresis ─────
        if (d_total[k] > pos_enter) {
            if (!pos_active) {
                pos_active = 1;
                start_idx = k;
                peak_val = max_ax[k];
                peak_idx = k;
            } else if (max_ax[k] > peak_val) {
                peak_val = max_ax[k];
                peak_idx = k;
            }
        } 
        else if (pos_active && d_total[k] < pos_exit) {
            pos_active = 0;
            int duration = k - start_idx;
            if (duration >= N_PERSIST) {

                // ─── centered search ±2 windows around detected peak ───
                int base = (peak_idx - 2) * DERIV_WINDOW_SIZE;
                if (base < 0) base = 0;

                int end = (peak_idx + 2) * DERIV_WINDOW_SIZE;
                if (end > n) end = n;

                double local_max = -1e9;
                int best_i = base;
                for (int i = base; i < end; i++) {
                    if (fabs(data[i].Ax) > local_max) {
                        local_max = fabs(data[i].Ax);
                        best_i = i;
                    }
                }
                printf("Positive spike at %s (%.3f m/s²)\n",
                       data[best_i].time, local_max);
                       addDerivativeDetection("Positive", data[best_i].time , local_max );
            }
        }

        // ───── Negative region detection with hysteresis ─────
        if (d_total[k] < neg_enter) {
            if (!neg_active) {
                neg_active = 1;
                start_idx = k;
                peak_val = max_ax[k];
                peak_idx = k;
            } else if (max_ax[k] > peak_val) {
                peak_val = max_ax[k];
                peak_idx = k;
            }
        }
        else if (neg_active && d_total[k] > neg_exit) {
            neg_active = 0;
            int duration = k - start_idx;
            if (duration >= N_PERSIST) {

                // ─── centered search ±1 windows around detected peak ───
                int base = (peak_idx - 1) * DERIV_WINDOW_SIZE;
                if (base < 0) base = 0;

                int end = (peak_idx + 1) * DERIV_WINDOW_SIZE;
                if (end > n) end = n;

                double local_max = -1e9;
                int best_i = base;
                for (int i = base; i < end; i++) {
                    if (fabs(data[i].Ax) > local_max) {
                        local_max = fabs(data[i].Ax);
                        best_i = i;
                    }
                }
                printf("Negative spike at %s (-%.3f m/s²)\n",
                       data[best_i].time, local_max);
                        addDerivativeDetection("Negetive", data[best_i].time , local_max );
            }
        }
    }

    printf("\nFinal thresholds: pos=%.3f neg=%.3f\n", pos_enter, neg_enter);
    free(d_total);
    free(max_ax);
}

static Detection derivDet[MAX_DETECTIONS];
static int derivCount = 0;

static void addDerivativeDetection(const char *type, const char *timeStr, float amplitude) {
    if (derivCount >= MAX_DETECTIONS) return;
    Detection *d = &derivDet[derivCount++];
    strcpy(d->polarity, type);
    int h, m, s, ms;
    sscanf(timeStr, "%d:%d:%d:%d", &h, &m, &s, &ms);
    d->hour = h; d->min = m; d->sec = s; d->ms = ms;
    d->amplitude = amplitude;
}



int getDerivativeDetections(Detection **array) {
    *array = derivDet;
    return derivCount;
}