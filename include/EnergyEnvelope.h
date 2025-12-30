#ifndef ENERGY_ENVELOPE_H
#define ENERGY_ENVELOPE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Detections.h"


// ─────────────────────────────────────────────
// Configuration Parameters
// ─────────────────────────────────────────────
#define FCFR_WEIGHT 0.8
#define MAX_SAMPLES_ENVELOPE 100000
#define FS 50.0                  // Sampling frequency (Hz)
#define FFT_SIZE 256            // Power of 2 ≥ window size (175 < 512)
#define FCFR_THRESHOLD 3.0       // Detection threshold
#define WINDOW_SAMPLES 256    // 3.5 s window for 50 Hz sampling
#define EPS 1e-9                 // Small epsilon for numerical safety

// ─────────────────────────────────────────────
// Data Structure for IMU Sample
// ─────────────────────────────────────────────
typedef struct {
    char time[32];
    double Ax, Ay, Az, Gx, Gy, Gz;
} EnvelopeSample;

// ─────────────────────────────────────────────
// Function Prototypes
// ─────────────────────────────────────────────
int loadEnvelopeCSV(const char *filename, EnvelopeSample *data);
void computeTeagerEnergy(const double *x, double *E, int n);
void computeImprovedEnvelopeWindow(const EnvelopeSample *data, int startIdx, int n, double fs);
void runImprovedEnvelopeSliding(const EnvelopeSample *data, int totalSamples, double fs);
int getEnergyDetections(Detection **array);

#endif
