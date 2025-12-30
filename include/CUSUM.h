#ifndef CUSUM_H
#define CUSUM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Detections.h"

// ========================== PARAMETERS ==========================
#define LINE_LEN 256
#define REF_SIZE 200             // reference window size
#define EXCLUDE_RECENT 10        // number of recent samples excluded from baseline
#define CUSUM_EPS 1e-9

#define DELTA_MAG 0.15           // mean shift magnitude (tune to sensitivity)
#define EVENT_LEN 120             // expected duration of a spike (samples)
#define ALPHA_BASE 0.05
#define ALPHA_FACTOR 1.2
#define RESET_DELAY 400          // samples to ignore after a detection

#define SMOOTH_WIN_CENTERED 7    // odd, centered smoothing window (no delay)
#define PEAK_PRE 80              // search window before crossing
#define PEAK_POST 120             // search window after crossing

// ========================== STRUCTURES ==========================
typedef struct {
    double data[REF_SIZE];
    int idx;
    int full;
} Window;

// ========================== FUNCTION DECLARATIONS ==========================
void add_to_window(Window *w, double value);
double window_mean(Window *w);
double window_var(Window *w, double mu);
double window_mean_exclude_recent(Window *w, int exclude);
void run_cusum_on_csv(const char *filename);

int getCUSUMDetections(Detection **array);

#endif
