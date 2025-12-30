#ifndef CWT_H
#define CWT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Detections.h"

#define MAX_LINES 2000000
#define MAX_TIMESTR_LEN 32

// --- Adjustable parameters ---
#define MIN_DUR_SEC 0.005     // 5 ms
#define MAX_DUR_SEC 0.050     // 50 ms
#define N_SCALES 8
#define MAD_SCALE_FACTOR 0.6745
#define THRESH_K 3.0          // threshold multiplier
#define MERGE_WINDOW_MS 175.0  // merge events closer than this
#define MAX_WAVELET_LEN 801   // must be odd and large enough

// --- Structs ---
typedef struct {
    char time_str[MAX_TIMESTR_LEN];
    double t_seconds;
    double ax, ay, az;
    double amp; // sqrt(ax^2+ay^2+az^2)
} Sample;

// --- Function prototypes ---
int run_cwt_to_csv(const char *filename);
int getCWTDetections(Detection **array);
#endif
