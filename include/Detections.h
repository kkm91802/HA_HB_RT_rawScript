#ifndef DETECTIONS_H
#define DETECTIONS_H

#define MAX_DETECTIONS 1000

typedef struct {
    char polarity[10];   // "Positive" or "Negative"
    int hour, min, sec, ms;
    float amplitude;     // Ax value (m/s²)
} Detection;
void addCUSUMDetection(const char *type, const char *timeStr, float amp);
static void addEnergyDetection(const char *type, const char *timeStr, float amplitude);
static void addCWTDetection(const char *type, const char *timeStr, float amplitude);
static void addDerivativeDetection(const char *type, const char *timeStr, float amplitude);

#endif
