// Butterworth.c
// Implementation of FFTW3-based automatic Butterworth filtering
// Requires linking with: -lfftw3 -lm

#define _GNU_SOURCE
#include "Butterworth.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

// -------------------- Internal Constants --------------------
#define OUTLIER_FACTOR 5
#define MIN_FC 0.5
#define DEFAULT_M 100

// -------------------- Internal Helpers --------------------
static int next_pow2(int v) { int p = 1; while (p < v) p <<= 1; return p; }

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x < y) ? -1 : (x > y);
}

// Parse time "HH:MM:SS:ms" -> seconds
static int parse_time_hhmmssms(const char *s, double *out_seconds) {
    if (!s || !out_seconds) return 0;
    int hh = 0, mm = 0, ss = 0;
    char fracbuf[16] = "0";
    const char *p = s;
    int n = 0;

    if (sscanf(p, "%d:%n", &hh, &n) != 1) return 0;
    p += n;
    if (sscanf(p, "%d:%n", &mm, &n) != 1) return 0;
    p += n;
    if (sscanf(p, "%d:%n", &ss, &n) != 1) return 0;
    p += n;

    int i = 0;
    while (*p && *p != ',' && *p != '\n' && *p != '\r' && i < 15) {
        if (*p >= '0' && *p <= '9') fracbuf[i++] = *p;
        else break;
        p++;
    }
    fracbuf[i] = '\0';

    double ms = 0.0;
    int len = strlen(fracbuf);
    if (len <= 3) ms = atof(fracbuf);
    else ms = atof(fracbuf) / pow(10.0, len - 3);

    *out_seconds = hh * 3600.0 + mm * 60.0 + ss + ms / 1000.0;
    return 1;
}

// -------------------- Sampling Frequency --------------------
double estimate_fs_from_timestamps(double *t, int count) {
    if (!t || count < 2) return 0.0;
    double *dt = malloc(sizeof(double) * (count - 1));
    if (!dt) return 0.0;

    int m = 0;
    for (int i = 0; i + 1 < count; i++) {
        double d = t[i + 1] - t[i];
        if (d < 0) d += 24.0 * 3600.0;
        if (d > 0) dt[m++] = d;
    }

    if (m == 0) { free(dt); return 0.0; }

    qsort(dt, m, sizeof(double), cmp_double);
    double median_dt = (m % 2) ? dt[m / 2] : 0.5 * (dt[m / 2 - 1] + dt[m / 2]);
    double fs = (median_dt > 0.0) ? 1.0 / median_dt : 0.0;

    free(dt);
    return fs;
}

// -------------------- FFT-based cutoff estimation --------------------
double estimate_fc_from_samples_fftw(const double *sig, int M, double fs, double energy_thresh) {
    if (!sig || M <= 2 || fs <= 0.0) return 0.0;

    int NFFT = next_pow2(M);
    double *in = fftw_malloc(sizeof(double) * NFFT);
    fftw_complex *out = fftw_malloc(sizeof(fftw_complex) * (NFFT / 2 + 1));
    fftw_plan plan = fftw_plan_dft_r2c_1d(NFFT, in, out, FFTW_ESTIMATE);

    double mean = 0;
    for (int i = 0; i < M; i++) mean += sig[i];
    mean /= M;

    for (int i = 0; i < M; i++) {
        double w = 0.5 * (1 - cos(2 * M_PI * i / (M - 1)));
        in[i] = (sig[i] - mean) * w;
    }
    for (int i = M; i < NFFT; i++) in[i] = 0;

    fftw_execute(plan);

    int K = NFFT / 2;
    double *P = malloc(sizeof(double) * (K + 1));
    double total = 0.0;
    for (int k = 0; k <= K; k++) {
        double re = out[k][0], im = out[k][1];
        double mag2 = re * re + im * im;
        if (k == 0 || (k == K && NFFT % 2 == 0))
            total += mag2;
        else
            total += 2 * mag2;
        P[k] = mag2;
    }

    double cum = 0.0, fc = 0.0;
    for (int k = 0; k <= K; k++) {
        double add = (k == 0 || (k == K && NFFT % 2 == 0)) ? P[k] : 2 * P[k];
        cum += add;
        if (cum / total >= energy_thresh) {
            fc = k * (fs / NFFT);
            break;
        }
    }

    if (fc < MIN_FC) fc = MIN_FC;
    if (fc > 0.4 * fs) fc = 0.4 * fs;

    free(P);
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return fc/32;
}

// -------------------- Auto order selection --------------------
int auto_select_order(double fc, double fs) {
    if (fs <= 0.0 || fc <= 0.0) return 4;  // Fallback safe value

    // --- Assumptions ---
    double As = 20.0;          // desired stopband attenuation in dB
    double fs_dash = 2.0 * fc; // stopband frequency (usually 2×fc)

    // --- Protect against invalid ratios ---
    if (fs_dash <= fc) fs_dash = fc * 1.5;  // ensure fs' > fc

    // --- Compute filter order ---
    double numerator = log10(pow(10.0, As / 10.0) - 1.0);
    double denominator = 2.0 * log10(fs_dash / fc);

    int n = (int)ceil(numerator / denominator); // round up to integer order

    // --- Clamp to safe practical range (1–10) ---
    if (n < 1) n = 1;
    if (n > 10) n = 10;

    return n ;
}

// -------------------- Butterworth filter design --------------------
static void designButterworthLPF(ButterworthLPF *F, int order, double fs, double fc) {
    int nSections = order / 2;
    double warped = tan(M_PI * fc / fs);
    F->nSections = nSections;

    for (int i = 0; i < nSections; i++) {
        double theta = M_PI * (2.0 * i + 1.0 + order - 1.0) / (2.0 * order);
        double sin_t = sin(theta);
        double a = 1.0 + sqrt(2.0) * warped * sin_t + warped * warped;
        double b0 = warped * warped / a;
        double b1 = 2.0 * b0;
        double b2 = b0;
        double a1 = 2.0 * (warped * warped - 1.0) / a;
        double a2 = (1.0 - sqrt(2.0) * warped * sin_t + warped * warped) / a;

        F->sec[i].b0 = b0; F->sec[i].b1 = b1; F->sec[i].b2 = b2;
        F->sec[i].a0 = 1.0; F->sec[i].a1 = a1; F->sec[i].a2 = a2;
        F->sec[i].x1 = F->sec[i].x2 = F->sec[i].y1 = F->sec[i].y2 = 0.0;
    }
}

static inline double applyBiquad(Biquad *q, double x) {
    double y = q->b0 * x + q->b1 * q->x1 + q->b2 * q->x2
             - q->a1 * q->y1 - q->a2 * q->y2;
    q->x2 = q->x1; q->x1 = x;
    q->y2 = q->y1; q->y1 = y;
    return y;
}

static double applyButterworth(ButterworthLPF *F, double x) {
    double y = x;
    for (int i = 0; i < F->nSections; i++)
        y = applyBiquad(&F->sec[i], y);
    return y;
}

// -------------------- Full processing pipeline --------------------
int run_auto_butterworth(const char *inputPath,
                         const char *outputPath,
                         int M_first_window,
                         double energy_thresh)
{
    if (!inputPath || !outputPath) return 1;
    if (M_first_window < 8) M_first_window = DEFAULT_M;
    if (energy_thresh <= 0.0 || energy_thresh >= 1.0) energy_thresh = 0.95;

    FILE *fin = fopen(inputPath, "r");
    if (!fin) { perror("open input"); return 2; }
    FILE *fout = fopen(outputPath, "w");
    if (!fout) { perror("open output"); fclose(fin); return 3; }

    char *line = NULL; size_t len = 0;
    if (getline(&line, &len, fin) == -1) { fclose(fin); fclose(fout); return 4; }
    fprintf(fout, "%s", line); // header

    // Allocate buffers
    int M = M_first_window;
    double *Ax = calloc(M, sizeof(double));
    double *Ay = calloc(M, sizeof(double));
    double *Az = calloc(M, sizeof(double));
    double *Gx = calloc(M, sizeof(double));
    double *Gy = calloc(M, sizeof(double));
    double *Gz = calloc(M, sizeof(double));
    double *timestamps = calloc(M, sizeof(double));
    char **timestr = calloc(M, sizeof(char*));

    int count = 0;
    while (count < M && getline(&line, &len, fin) != -1) {
        char ts[64]; double a1,a2,a3,g1,g2,g3;
        if (sscanf(line, " %63[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                   ts, &a1,&a2,&a3,&g1,&g2,&g3) == 7) {
            timestr[count] = strdup(ts);
            Ax[count]=a1; Ay[count]=a2; Az[count]=a3;
            Gx[count]=g1; Gy[count]=g2; Gz[count]=g3;
            double tsec;
            parse_time_hhmmssms(ts, &tsec);
            timestamps[count] = tsec;
            count++;
        }
    }

    double fs = estimate_fs_from_timestamps(timestamps, count);
    if (!(fs > 0)) fs = 50.0;

    double fc[6];
    fc[0] = estimate_fc_from_samples_fftw(Ax, count, fs, energy_thresh);
    fc[1] = estimate_fc_from_samples_fftw(Ay, count, fs, energy_thresh);
    fc[2] = estimate_fc_from_samples_fftw(Az, count, fs, energy_thresh);
    fc[3] = estimate_fc_from_samples_fftw(Gx, count, fs, energy_thresh);
    fc[4] = estimate_fc_from_samples_fftw(Gy, count, fs, energy_thresh);
    fc[5] = estimate_fc_from_samples_fftw(Gz, count, fs, energy_thresh);

    double fc_final = fc[0];
    for (int i = 1; i < 6; i++) if (fc[i] > fc_final) fc_final = fc[i];
    int order = auto_select_order(fc_final, fs);

    ButterworthLPF filters[6];
    for (int i = 0; i < 6; i++)
        designButterworthLPF(&filters[i], order, fs, fc_final);

    // Process first M samples
    for (int i = 0; i < count; i++) {
        double fx = applyButterworth(&filters[0], Ax[i]);
        double fy = applyButterworth(&filters[1], Ay[i]);
        double fz = applyButterworth(&filters[2], Az[i]);
        double gx = applyButterworth(&filters[3], Gx[i]);
        double gy = applyButterworth(&filters[4], Gy[i]);
        double gz = applyButterworth(&filters[5], Gz[i]);
        fprintf(fout, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                timestr[i], fx, fy, fz, gx, gy, gz);
    }

    // Continue filtering the rest of the file
    while (getline(&line, &len, fin) != -1) {
        char ts[64]; double a1,a2,a3,g1,g2,g3;
        if (sscanf(line, " %63[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                   ts, &a1,&a2,&a3,&g1,&g2,&g3) == 7) {
            double fx = applyButterworth(&filters[0], a1);
            double fy = applyButterworth(&filters[1], a2);
            double fz = applyButterworth(&filters[2], a3);
            double gx = applyButterworth(&filters[3], g1);
            double gy = applyButterworth(&filters[4], g2);
            double gz = applyButterworth(&filters[5], g3);
            fprintf(fout, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    ts, fx, fy, fz, gx, gy, gz);
        }
    }

    printf("Estimated fs = %.3f Hz, fc_final = %.3f Hz, order = %d\n", fs, fc_final, order);

    // Cleanup
    for (int i = 0; i < count; i++) free(timestr[i]);
    free(timestr); free(Ax); free(Ay); free(Az); free(Gx); free(Gy); free(Gz); free(timestamps);
    free(line); fclose(fin); fclose(fout);
    return 0;
}
