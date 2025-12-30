#ifndef BUTTERWORTH_H
#define BUTTERWORTH_H

#include <stdio.h>

// -------------------- Struct Definitions --------------------
#define MAX_ORDER 16

typedef struct {
    double b0, b1, b2;
    double a0, a1, a2;
    double x1, x2;
    double y1, y2;
} Biquad;

typedef struct {
    int nSections;
    Biquad sec[MAX_ORDER / 2];
} ButterworthLPF;

// -------------------- Public API --------------------

// Main pipeline function
int run_auto_butterworth(const char *inputPath,
                         const char *outputPath,
                         int M_first_window,   // number of samples for FFT-based fc
                         double energy_thresh); // e.g. 0.95

// Helper (optional): estimate fs directly from timestamps in an array
double estimate_fs_from_timestamps(double *t, int count);

// Helper (optional): estimate fc for one axis using FFT
double estimate_fc_from_samples_fftw(const double *sig, int M, double fs, double energy_thresh);

// Helper (optional): auto-select filter order from ratio 2*fc/fs
int auto_select_order(double fc, double fs);

#endif
