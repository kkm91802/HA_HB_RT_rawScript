#include "Threshold.h"
#include "Detections.h"
// #include "Polling.h" 

static Detection threshDet[MAX_DETECTIONS];
static int threshCount = 0;

static void addThresholdDetection(const char *type, const char *timeStr, float amplitude) {
    if (threshCount >= MAX_DETECTIONS) return;
    Detection *d = &threshDet[threshCount++];
    strcpy(d->polarity, type);
    int h, m, s, ms;
    sscanf(timeStr, "%d:%d:%d:%d", &h, &m, &s, &ms);
    d->hour = h; d->min = m; d->sec = s; d->ms = ms;
    d->amplitude = amplitude;
}

// ------------------- Helper: Convert HH:MM:SS:ms → seconds -------------------
static double parseTimeToSeconds(const char *timeStr) {
    int h, m, s, ms;
    if (sscanf(timeStr, "%d:%d:%d:%d", &h, &m, &s, &ms) == 4) {
        return h * 3600.0 + m * 60.0 + s + (ms / 1000.0);
    }
    return 0.0;
}

// ------------------- Adaptive Threshold Functions -------------------
void initThreshold(AdaptiveThreshold *th, float posInit, float negInit, int timeTh) {
    th->posThreshold = posInit;
    th->negThreshold = negInit;
    th->timeThreshold = timeTh;
    th->lastSpikeIndex = -timeTh;
}

void updateAdaptiveThreshold(AdaptiveThreshold *th, const float *window, int windowSize) {
    float minVal = window[0], maxVal = window[0];
    for (int i = 1; i < windowSize; ++i) {
        if (window[i] < minVal) minVal = window[i];
        if (window[i] > maxVal) maxVal = window[i];
    }
    // th->posThreshold = 0.8f * th->posThreshold + 0.2f * maxVal;
    // th->negThreshold = 0.8f * th->negThreshold + 0.2f * minVal;
    th->posThreshold = (0.95f * th->posThreshold + 0.05f * maxVal) ;
    th->negThreshold = (0.95f * th->negThreshold + 0.05f * minVal) - 0.5f;

}

bool detectPositiveSpike(const float *signal, int i, AdaptiveThreshold *th) {
    if (i < 1) return false;
    if (signal[i] > th->posThreshold &&
        signal[i - 1] < signal[i] &&
        signal[i + 1] < signal[i]) {
        if (i - th->lastSpikeIndex > th->timeThreshold) {
            th->lastSpikeIndex = i;
            return true;
        }
    }
    return false;
}

bool detectNegativeSpike(const float *signal, int i, AdaptiveThreshold *th) {
    if (i < 1) return false;
    if (signal[i] < th->negThreshold &&
        signal[i - 1] > signal[i] &&
        signal[i + 1] > signal[i]) {
        if (i - th->lastSpikeIndex > th->timeThreshold) {
            th->lastSpikeIndex = i;
            return true;
        }
    }
    return false;
}

// ------------------- File Parsing & Magnitude Computation -------------------
int readIMUCSV(const char *filename, IMUSample *data, int maxRows) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        return -1;
    }

    char line[256];
    int count = 0;
    fgets(line, sizeof(line), fp);  // skip header

    while (fgets(line, sizeof(line), fp) && count < maxRows) {
        char timeStr[16];
        sscanf(line, "%15[^,],%f,%f,%f,%f,%f,%f",
               timeStr,
               &data[count].ax, &data[count].ay, &data[count].az,
               &data[count].gx, &data[count].gy, &data[count].gz);

        strcpy(data[count].timeStr, timeStr);
        data[count].timeSec = parseTimeToSeconds(timeStr);
        count++;
    }

    fclose(fp);
    return count;
}

void computeMagnitude(const IMUSample *data, float *accMag, int n) {
    for (int i = 0; i < n; i++) {
        // accMag[i] = sqrtf(data[i].ax * data[i].ax +
        //                   data[i].ay * data[i].ay +
        //                   data[i].az * data[i].az);
                accMag[i] = data[i].ax ;
    }
}

int run_Threshold(const char *input_csv)
{
    IMUSample dataX[MAX_SAMPLES];
    float accMag[MAX_SAMPLES];
    AdaptiveThreshold th;

    int numSamples = readIMUCSV(input_csv, dataX, MAX_SAMPLES);
    if (numSamples <= 0) {
        printf("No data read or file error in Threshold().\n");
        return -1;
    }

    computeMagnitude(dataX, accMag, numSamples);
    initThreshold(&th, INIT_THRESHOLD_POS, INIT_THRESHOLD_NEG, TIME_THRESHOLD);

    printf("Processing %d samples...\n", numSamples);

    for (int i = 1; i < numSamples - 1; i++) {

        if (i > WINDOW_SIZE && i % WINDOW_SIZE == 0)
            updateAdaptiveThreshold(&th, &accMag[i - WINDOW_SIZE], WINDOW_SIZE);

        if (detectPositiveSpike(accMag, i, &th)) {
            printf("Positive spike at %s (%.3f m/s²)\n", dataX[i].timeStr, accMag[i]);
            addThresholdDetection("Positive", dataX[i].timeStr, accMag[i]);
        }

        if (detectNegativeSpike(accMag, i, &th)) {
            printf("Negative spike at %s (%.3f m/s²)\n", dataX[i].timeStr, accMag[i]);
            addThresholdDetection("Negative", dataX[i].timeStr, accMag[i]);
        }
    }

    printf("\nFinal thresholds: pos=%.3f neg=%.3f\n",
           th.posThreshold, th.negThreshold);

    return 0;
}
 // Include SpikeEvent struct definition

// int run_threshold_detection(const char *filename, SpikeEvent *events, size_t max_events) {
//     // Load the IMU data
//     IMUSample data[MAX_SAMPLES_IMU];
//     int totalSamples = readIMUCSV(filename, data, MAX_SAMPLES_IMU);
//     if (totalSamples <= 0) {
//         fprintf(stderr, "No valid data loaded.\n");
//         return 0;
//     }

//     // Initialize threshold parameters
//     AdaptiveThreshold th;
//     initThreshold(&th, 0.5f, -0.5f, 10);  // Example initialization

//     // Compute the magnitude of the accelerometer signal (you can modify it to use all 3 axes or a custom calculation)
//     float accMag[MAX_SAMPLES_IMU];
//     computeMagnitude(data, accMag, totalSamples);

//     // Event counter for detected spikes
//     int event_count = 0;

//     // Loop through the signal and apply the threshold detection
//     for (int i = 1; i < totalSamples - 1 && event_count < max_events; i++) {
//         // Update the adaptive threshold with the current window
//         float window[WINDOW_SIZE];
//         for (int j = 0; j < WINDOW_SIZE; j++) {
//             if (i - j >= 0) {
//                 window[j] = accMag[i - j];
//             }
//         }
//         updateAdaptiveThreshold(&th, window, WINDOW_SIZE);

//         // Detect positive and negative spikes
//         if (detectPositiveSpike(accMag, i, &th)) {
//             // Store detected positive spike in events[] array
//             strncpy(events[event_count].time_str, data[i].timeStr, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = data[i].timeSec * 1000.0;  // Convert to milliseconds
//             events[event_count].magnitude = accMag[i];
//             events[event_count].type = +1;  // Positive spike
//             event_count++;
//         }

//         if (detectNegativeSpike(accMag, i, &th)) {
//             // Store detected negative spike in events[] array
//             strncpy(events[event_count].time_str, data[i].timeStr, sizeof(events[event_count].time_str) - 1);
//             events[event_count].time_str[sizeof(events[event_count].time_str) - 1] = '\0';
//             events[event_count].time_ms = data[i].timeSec * 1000.0;
//             events[event_count].magnitude = accMag[i];
//             events[event_count].type = -1;  // Negative spike
//             event_count++;
//         }
//     }

//     return event_count;  // Return the number of detected events
// }




int getThresholdDetections(Detection **array) {
    *array = threshDet;
    return threshCount;
}