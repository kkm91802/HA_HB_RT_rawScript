#ifndef POLLING_H
#define POLLING_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define MAX_SAMPLES_IMU 10000 
#define MAX_EVENTS 200
#define TIME_TOLERANCE_MS 300   // ±300 ms allowed difference
#define MIN_AGREE 3             // at least 3 algorithms must agree

/* A single detected spike event (shared type for all detectors) */
typedef struct {
    char time_str[32];    /* original timestamp string (e.g., "12:01:29:379") */
    double time_ms;       /* time in milliseconds since midnight (or since start); use consistent meaning */
    double magnitude;     /* amplitude (m/s^2) */
    int type;             /* +1 = positive spike, -1 = negative spike */
} SpikeEvent;

int run_cusum_detection(const char *filename, SpikeEvent *events, size_t max_events);

/* (Repeat prototypes for other detectors as you adapt them) */
int run_cwt_detection(const char *filename, SpikeEvent *events, size_t max_events);
int run_derivative_detection(const char *filename, SpikeEvent *events, size_t max_events);
int run_energy_detection(const char *filename, SpikeEvent *events, size_t max_events);
int run_threshold_detection(const char *filename, SpikeEvent *events, size_t max_events);
static long time_to_ms(const char *time_str);
static int events_match(const SpikeEvent *a, const SpikeEvent *b) ;
int run_polling_detection(const char *filename, SpikeEvent *final_events, size_t max_events);

#endif /* POLLING_H */
