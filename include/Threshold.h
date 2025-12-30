#ifndef THRESHOLD_H
#define THRESHOLD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "Detections.h"


#define TIME_THRESHOLD 175  //No. of entries 
#define INIT_THRESHOLD_POS 0.0f
#define INIT_THRESHOLD_NEG 0.0
#define WINDOW_SIZE 100
#define MAX_SAMPLES 10000

// -------- Data Structures --------
typedef struct {
    char timeStr[16];   // e.g., "12:01:20:681"
    double timeSec;     // time converted to seconds
    float ax, ay, az, gx, gy, gz;
} IMUSample;

typedef struct {
    float posThreshold;
    float negThreshold;
    int timeThreshold;
    int lastSpikeIndex;
} AdaptiveThreshold;

// -------- Function Declarations --------
void initThreshold(AdaptiveThreshold *th, float posInit, float negInit, int timeTh);
void updateAdaptiveThreshold(AdaptiveThreshold *th, const float *window, int windowSize);
bool detectPositiveSpike(const float *signal, int i, AdaptiveThreshold *th);
bool detectNegativeSpike(const float *signal, int i, AdaptiveThreshold *th);

int readIMUCSV(const char *filename, IMUSample *data, int maxRows);
void computeMagnitude(const IMUSample *data, float *accMag, int n);

int run_Threshold(const char *input_csv);
int getThresholdDetections(Detection **array);

#endif
