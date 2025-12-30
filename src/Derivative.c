#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Derivative.h"

// Function to convert time string "HH:MM:SS:msmsms" to seconds
double timeStringToSeconds(const char* timeStr) {
    int hours, minutes, seconds, milliseconds;
    
    if (sscanf(timeStr, "%d:%d:%d:%d", &hours, &minutes, &seconds, &milliseconds) == 4) {
        return hours * 3600.0 + minutes * 60.0 + seconds + milliseconds / 1000.0;
    }
    
    return 0.0;
}

// Function to read CSV file with time string parsing
int readCSV(const char* filename, DataPoint data[]) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error: Could not open file %s\n", filename);
        return -1;
    }

    char line[MAX_LINE_LENGTH];
    int row_count = 0;
    int is_first_line = 1;

    while (fgets(line, sizeof(line), file) && row_count < MAX_ROWS) {
        // Skip header line
        if (is_first_line) {
            is_first_line = 0;
            continue;
        }

        // Remove newline character
        line[strcspn(line, "\n")] = 0;

        // Parse CSV line
        char* token = strtok(line, ",");
        if (token == NULL) continue;

        // Store original time string
        strncpy(data[row_count].time_str, token, sizeof(data[row_count].time_str) - 1);
        data[row_count].time_str[sizeof(data[row_count].time_str) - 1] = '\0';

        // Convert time string to seconds for calculations
        data[row_count].time = timeStringToSeconds(token);
        
        token = strtok(NULL, ",");
        data[row_count].ax = atof(token);
        
        token = strtok(NULL, ",");
        data[row_count].ay = atof(token);
        
        token = strtok(NULL, ",");
        data[row_count].az = atof(token);
        
        token = strtok(NULL, ",");
        data[row_count].gx = atof(token);
        
        token = strtok(NULL, ",");
        data[row_count].gy = atof(token);
        
        token = strtok(NULL, ",");
        data[row_count].gz = atof(token);

        row_count++;
    }

    fclose(file);
    return row_count;
}

// Function to calculate derivatives using proper time differences
void calculateDerivatives(DataPoint data[], int data_count, Derivative derivatives[]) {
    if (data_count < 2) {
        printf("Error: Need at least 2 data points to calculate derivatives\n");
        return;
    }

    // Calculate derivatives for each point
    for (int i = 0; i < data_count; i++) {
        // Copy original time string
        strncpy(derivatives[i].time_str, data[i].time_str, sizeof(derivatives[i].time_str) - 1);
        derivatives[i].time_str[sizeof(derivatives[i].time_str) - 1] = '\0';

        if (i == 0) {
            // For the first point, use forward difference
            double dt = data[i+1].time - data[i].time;
            if (dt > 0) {
                derivatives[i].dax = (data[i+1].ax - data[i].ax) / dt;
                derivatives[i].day = (data[i+1].ay - data[i].ay) / dt;
                derivatives[i].daz = (data[i+1].az - data[i].az) / dt;
                derivatives[i].dgx = (data[i+1].gx - data[i].gx) / dt;
                derivatives[i].dgy = (data[i+1].gy - data[i].gy) / dt;
                derivatives[i].dgz = (data[i+1].gz - data[i].gz) / dt;
            } else {
                derivatives[i].dax = derivatives[i].day = derivatives[i].daz = 0;
                derivatives[i].dgx = derivatives[i].dgy = derivatives[i].dgz = 0;
            }
        } else {
            // For all other points: (current - previous) / (current_time - previous_time)
            double dt = data[i].time - data[i-1].time;
            if (dt > 0) {
                derivatives[i].dax = (data[i].ax - data[i-1].ax) / dt;
                derivatives[i].day = (data[i].ay - data[i-1].ay) / dt;
                derivatives[i].daz = (data[i].az - data[i-1].az) / dt;
                derivatives[i].dgx = (data[i].gx - data[i-1].gx) / dt;
                derivatives[i].dgy = (data[i].gy - data[i-1].gy) / dt;
                derivatives[i].dgz = (data[i].gz - data[i-1].gz) / dt;
            } else {
                derivatives[i].dax = derivatives[i].day = derivatives[i].daz = 0;
                derivatives[i].dgx = derivatives[i].dgy = derivatives[i].dgz = 0;
            }
        }
    }
}

// Function to calculate double derivatives using proper time differences
void calculateDoubleDerivatives(DataPoint data[], Derivative derivatives[], int data_count, DoubleDerivative double_derivatives[]) {
    if (data_count < 2) {
        printf("Error: Need at least 2 data points to calculate double derivatives\n");
        return;
    }

    // Calculate double derivatives for each point
    for (int i = 0; i < data_count; i++) {
        // Copy original time string
        strncpy(double_derivatives[i].time_str, data[i].time_str, 
                sizeof(double_derivatives[i].time_str) - 1);
        double_derivatives[i].time_str[sizeof(double_derivatives[i].time_str) - 1] = '\0';

        if (i == 0) {
            // For the first point, use forward difference
            double dt = data[i+1].time - data[i].time;
            if (dt > 0) {
                double_derivatives[i].ddax = (derivatives[i+1].dax - derivatives[i].dax) / dt;
                double_derivatives[i].dday = (derivatives[i+1].day - derivatives[i].day) / dt;
                double_derivatives[i].ddaz = (derivatives[i+1].daz - derivatives[i].daz) / dt;
                double_derivatives[i].ddgx = (derivatives[i+1].dgx - derivatives[i].dgx) / dt;
                double_derivatives[i].ddgy = (derivatives[i+1].dgy - derivatives[i].dgy) / dt;
                double_derivatives[i].ddgz = (derivatives[i+1].dgz - derivatives[i].dgz) / dt;
            } else {
                double_derivatives[i].ddax = double_derivatives[i].dday = double_derivatives[i].ddaz = 0;
                double_derivatives[i].ddgx = double_derivatives[i].ddgy = double_derivatives[i].ddgz = 0;
            }
        } else {
            // For all other points: (current_derivative - previous_derivative) / (current_time - previous_time)
            double dt = data[i].time - data[i-1].time;
            if (dt > 0) {
                double_derivatives[i].ddax = (derivatives[i].dax - derivatives[i-1].dax) / dt;
                double_derivatives[i].dday = (derivatives[i].day - derivatives[i-1].day) / dt;
                double_derivatives[i].ddaz = (derivatives[i].daz - derivatives[i-1].daz) / dt;
                double_derivatives[i].ddgx = (derivatives[i].dgx - derivatives[i-1].dgx) / dt;
                double_derivatives[i].ddgy = (derivatives[i].dgy - derivatives[i-1].dgy) / dt;
                double_derivatives[i].ddgz = (derivatives[i].dgz - derivatives[i-1].dgz) / dt;
            } else {
                double_derivatives[i].ddax = double_derivatives[i].dday = double_derivatives[i].ddaz = 0;
                double_derivatives[i].ddgx = double_derivatives[i].ddgy = double_derivatives[i].ddgz = 0;
            }
        }
    }
}

// Function to write first derivatives to CSV file
void writeDerivativesCSV(const char* filename, Derivative derivatives[], int count) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not create output file %s\n", filename);
        return;
    }

    // Write header
    fprintf(file, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");

    // Write data with original time format
    for (int i = 0; i < count; i++) {
        fprintf(file, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                derivatives[i].time_str,
                derivatives[i].dax,
                derivatives[i].day,
                derivatives[i].daz,
                derivatives[i].dgx,
                derivatives[i].dgy,
                derivatives[i].dgz);
    }

    fclose(file);
    printf("First derivatives written to %s\n", filename);
}

// Function to write double derivatives to CSV file
void writeDoubleDerivativesCSV(const char* filename, DoubleDerivative double_derivatives[], int count) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error: Could not create output file %s\n", filename);
        return;
    }

    // Write header
    fprintf(file, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");

    // Write data with original time format
    for (int i = 0; i < count; i++) {
        fprintf(file, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                double_derivatives[i].time_str,
                double_derivatives[i].ddax,
                double_derivatives[i].dday,
                double_derivatives[i].ddaz,
                double_derivatives[i].ddgx,
                double_derivatives[i].ddgy,
                double_derivatives[i].ddgz);
    }

    fclose(file);
    printf("Double derivatives written to %s\n", filename);
}

// Function to print time analysis
void analyzeTimeData(DataPoint data[], int count) {
    printf("\nTime Analysis:\n");
    printf("First time point: %s (%.6f seconds)\n", data[0].time_str, data[0].time);
    printf("Last time point: %s (%.6f seconds)\n", data[count-1].time_str, data[count-1].time);
    printf("Total duration: %.6f seconds\n", data[count-1].time - data[0].time);
    
    if (count > 1) {
        printf("\nTime differences between consecutive points:\n");
        for (int i = 1; i < count && i < 10; i++) {  // Show first 10 differences
            double dt = data[i].time - data[i-1].time;
            printf("Point %d to %d: %s -> %s (%.6f seconds, %.3f ms)\n", 
                   i-1, i, data[i-1].time_str, data[i].time_str, dt, dt * 1000);
        }
        if (count > 10) {
            printf("... (showing first 10 differences)\n");
        }
    }
}