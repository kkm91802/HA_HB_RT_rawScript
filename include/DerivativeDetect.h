#ifndef DERIVATIVE_DETECT_H
#define DERIVATIVE_DETECT_H

#include <stdio.h>
#include "Detections.h"

#define MAX_SAMPLES_DERIV 100000
#define FS 50.0               // Sampling frequency (Hz)
#define T_S (1.0 / FS)
#define T_MIN 1.05          // Minimum duration (seconds)
// #ifdef WINDOW_SIZE
// #undef WINDOW_SIZE
// #endif
#define DERIV_WINDOW_SIZE 10      // Samples per window
#define N_PERSIST ((int)(T_MIN / (DERIV_WINDOW_SIZE * T_S) + 0.5))  // ≈18 windows

typedef struct {
    char  time[32];
    double Ax, Ay, Az;
    double Gx, Gy, Gz;
} DerivSample;

int  loadCSV(const char *filename, DerivSample *data);
void detectSpikes(const DerivSample *data, int n);
int getDerivativeDetections(Detection **array);

#endif


