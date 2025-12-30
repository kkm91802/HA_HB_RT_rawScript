#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#define MAX_LINE_LENGTH 1024
#define MAX_ROWS 10000

typedef struct {
    char time_str[20];  // Store original time string
    double time;        // Time in seconds for calculations
    double ax, ay, az;
    double gx, gy, gz;
} DataPoint;

typedef struct {
    char time_str[20];  // Keep original time format
    double dax, day, daz;
    double dgx, dgy, dgz;
} Derivative;

typedef struct {
    char time_str[20];  // Keep original time format
    double ddax, dday, ddaz;
    double ddgx, ddgy, ddgz;
} DoubleDerivative;

// Function declarations
double timeStringToSeconds(const char* timeStr);
int readCSV(const char* filename, DataPoint data[]);
void calculateDerivatives(DataPoint data[], int data_count, Derivative derivatives[]);
void calculateDoubleDerivatives(DataPoint data[], Derivative derivatives[], int data_count, DoubleDerivative double_derivatives[]);
void writeDerivativesCSV(const char* filename, Derivative derivatives[], int count);
void writeDoubleDerivativesCSV(const char* filename, DoubleDerivative double_derivatives[], int count);
void analyzeTimeData(DataPoint data[], int count);

#endif