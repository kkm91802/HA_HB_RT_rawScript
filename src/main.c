#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CSVConverter.h"
#include "Orientation.h"
#include "CalibrationMetrices.h"
#include "Calibration.h"
#include "CalibrationInitial.h"
#include "MassConverter.h"
#include "AutoCalibration.h"
#include "ButterworthLPF.h"
#include "Butterworth.h"
#include "InEKF.h"
#include "Derivative.h"
#include "Threshold.h"
#include "DerivativeDetect.h"
#include "EnergyEnvelope.h"
#include "CUSUM.h"
#include "CWT.h"
// #include "Polling.h"
#include "Detections.h"

#define NUM_MODULES 5

typedef struct {
    Detection *list;
    int count;
} DetectorOutput;


#define TIME_TOLERANCE_MS 300  // 0.3 second window

long detectionToMs(const Detection *d) {
    return ((d->hour * 3600L + d->min * 60L + d->sec) * 1000L) + d->ms;
}


static int isSameDetection(const Detection *a, const Detection *b) {
    if (strcmp(a->polarity, b->polarity) != 0)
        return 0;
    long diff = labs(detectionToMs(a) - detectionToMs(b));
    return diff <= TIME_TOLERANCE_MS;
}

int sameRoundedSecond(const Detection *a, const Detection *b) {
    if (strcmp(a->polarity, b->polarity) != 0)
        return 0;

    long ta = detectionToMs(a);
    long tb = detectionToMs(b);

    return labs(ta - tb) <= TIME_TOLERANCE_MS;  // ±500 ms window
}

void fuseDetections(void) {
    Detection *detLists[NUM_MODULES];
    int detCounts[NUM_MODULES];

    detCounts[0] = getCUSUMDetections(&detLists[0]);
    detCounts[1] = getEnergyDetections(&detLists[1]);
    detCounts[2] = getCWTDetections(&detLists[2]);
    detCounts[3] = getDerivativeDetections(&detLists[3]);
    detCounts[4] = getThresholdDetections(&detLists[4]);

    int fusedFlags[NUM_MODULES][1000] = {0}; // adjust 1000 to max per detector

    printf("\n--- Fusion Results ---\n");

    for (int i = 0; i < NUM_MODULES; i++) {
        for (int j = 0; j < detCounts[i]; j++) {
            if (fusedFlags[i][j]) continue;  // skip already fused

            Detection *base = &detLists[i][j];
            int votes = 1;
            float ampSum = base->amplitude;

            fusedFlags[i][j] = 1; // mark base as fused

            // Compare with all other detectors
            for (int k = 0; k < NUM_MODULES; k++) {
                if (k == i) continue;
                for (int m = 0; m < detCounts[k]; m++) {
                    if (!fusedFlags[k][m] && isSameDetection(base, &detLists[k][m])) {
                        votes++;
                        ampSum += detLists[k][m].amplitude;
                        fusedFlags[k][m] = 1; // mark match as fused
                        break; // count only first match per detector
                    }
                }
            }

            if (votes >= 3) {
                int hh = base->hour;
                int mm = base->min;
                int ss = base->sec;
                if (base->ms >= 500) {
                    ss++;
                    if (ss >= 60) { ss = 0; mm++; }
                    if (mm >= 60) { mm = 0; hh = (hh + 1) % 24; }
                }

                printf("%s spike at %02d:%02d:%02d  (consensus %d/%d, avg %.3f m/s²)\n",
                       base->polarity, hh, mm, ss, votes, NUM_MODULES, ampSum / votes);
            }
        }
    }
}

// void fuseDetections(void) {
//     Detection *detLists[5];
//     int detCounts[5];
//     detCounts[0] = getCUSUMDetections(&detLists[0]);
//     detCounts[1] = getEnergyDetections(&detLists[1]);
//     detCounts[2] = getCWTDetections(&detLists[2]);
//     detCounts[3] = getDerivativeDetections(&detLists[3]);
//     detCounts[4] = getThresholdDetections(&detLists[4]);

//     printf("\n--- Fusion Results ---\n");

//     // Iterate through all detections from all detectors
//     for (int i = 0; i < 5; i++) {
//         for (int j = 0; j < detCounts[i]; j++) {
//             Detection *base = &detLists[i][j];
//             int votes = 1;
//             float ampSum = base->amplitude;

//             // Compare with all other detectors
//             for (int k = i + 1; k < 5; k++) {
//                 for (int m = 0; m < detCounts[k]; m++) {
//                     if (isSameDetection(base, &detLists[k][m])) {
//                         votes++;
//                         ampSum += detLists[k][m].amplitude;
//                         break;
//                     }
//                 }
//             }

//             // Require at least 3 of 5 to agree
//             if (votes >= 3) {
//                 // Round milliseconds to nearest second
//                 int hh = base->hour;
//                 int mm = base->min;
//                 int ss = base->sec;
//                 if (base->ms >= 500) {
//                     ss++;
//                     if (ss >= 60) { ss = 0; mm++; }
//                     if (mm >= 60) { mm = 0; hh = (hh + 1) % 24; }
//                 }

//                 printf("%s spike at %02d:%02d:%02d  (consensus %d/5, avg %.3f m/s²)\n",
//                        base->polarity, hh, mm, ss, votes, ampSum / votes);
//             }
//         }
//     }
// }

#include <stdio.h>
#include "Detections.h"

void printAllDetections(void) {
    Detection *list;
    int count;

    printf("\n--- DEBUG: CUSUM Detections ---\n");
    count = getCUSUMDetections(&list);
    for (int i = 0; i < count; i++)
        printf("[%d] %s spike at %02d:%02d:%02d:%03d  (%.3f m/s²)\n",
               i, list[i].polarity, list[i].hour, list[i].min,
               list[i].sec, list[i].ms, list[i].amplitude);

    printf("\n--- DEBUG: Energy Envelope Detections ---\n");
    count = getEnergyDetections(&list);
    for (int i = 0; i < count; i++)
        printf("[%d] %s spike at %02d:%02d:%02d:%03d  (%.3f m/s²)\n",
               i, list[i].polarity, list[i].hour, list[i].min,
               list[i].sec, list[i].ms, list[i].amplitude);

    printf("\n--- DEBUG: CWT Detections ---\n");
    count = getCWTDetections(&list);
    for (int i = 0; i < count; i++)
        printf("[%d] %s spike at %02d:%02d:%02d:%03d  (%.3f m/s²)\n",
               i, list[i].polarity, list[i].hour, list[i].min,
               list[i].sec, list[i].ms, list[i].amplitude);

    printf("\n--- DEBUG: Derivative Detections ---\n");
    count = getDerivativeDetections(&list);
    for (int i = 0; i < count; i++)
        printf("[%d] %s spike at %02d:%02d:%02d:%03d  (%.3f m/s²)\n",
               i, list[i].polarity, list[i].hour, list[i].min,
               list[i].sec, list[i].ms, list[i].amplitude);

    printf("\n--- DEBUG: Threshold Detections ---\n");
    count = getThresholdDetections(&list);
    for (int i = 0; i < count; i++)
        printf("[%d] %s spike at %02d:%02d:%02d:%03d  (%.3f m/s²)\n",
               i, list[i].polarity, list[i].hour, list[i].min,
               list[i].sec, list[i].ms, list[i].amplitude);
}



// // Convert "HH:MM:SS" or "HH:MM:SS:MS" to milliseconds
// static long time_to_ms(const char *time_str) {
//     int hh = 0, mm = 0, ss = 0, ms = 0;
//     sscanf(time_str, "%d:%d:%d:%d", &hh, &mm, &ss, &ms);
//     return ((hh * 3600 + mm * 60 + ss) * 1000L) + ms;
// }

// // Check if two events are close in time and have the same sign
// static int events_match(const SpikeEvent *a, const SpikeEvent *b) {
//     return (a->type == b->type) &&
//            (fabs(a->time_ms - b->time_ms) <= TIME_TOLERANCE_MS);
// }

// int run_polling_detection(const char *filename, SpikeEvent *final_events, size_t max_events) {
//     // Storage for each algorithm’s events
//     SpikeEvent cusum[MAX_EVENTS], cwt[MAX_EVENTS],
//                deriv[MAX_EVENTS], energy[MAX_EVENTS],
//                thresh[MAX_EVENTS];

//     int o_cusum  = run_cusum_detection(filename, cusum, MAX_EVENTS);
//     int o_cwt    = run_cwt_detection(filename, cwt, MAX_EVENTS);
//     int o_deriv  = run_derivative_detection(filename, deriv, MAX_EVENTS);
//     int o_energy = run_energy_detection(filename, energy, MAX_EVENTS);
//     int o_thresh = run_threshold_detection(filename, thresh, MAX_EVENTS);

//     printf("\n=== Polling-based Consensus Detection (Positive & Negative) ===\n");

//     SpikeEvent *lists[5] = { cusum, cwt, deriv, energy, thresh };
//     int counts[5] = { o_cusum, o_cwt, o_deriv, o_energy, o_thresh };

//     int total_agree = 0;

//     // Iterate through all events from the reference algorithm (CUSUM)
//     for (int i = 0; i < o_cusum && total_agree < max_events; i++) {
//         int agree = 1; // include itself

//         for (int m = 1; m < 5; m++) {
//             for (int j = 0; j < counts[m]; j++) {
//                 if (events_match(&cusum[i], &lists[m][j])) {
//                     agree++;
//                     break;
//                 }
//             }
//         }

//         if (agree >= MIN_AGREE) {
//             printf("%s spike at %s (%.3f m/s²) — agreed by %d detectors\n",
//                    (cusum[i].type > 0) ? "Positive" : "Negative",
//                    cusum[i].time_str,
//                    cusum[i].magnitude,
//                    agree);
//             final_events[total_agree++] = cusum[i];
//         }
//     }

//     if (total_agree == 0)
//         printf("No common spikes detected within ±%d ms window.\n", TIME_TOLERANCE_MS);
//     else
//         printf("\nTotal %d consensus spikes detected.\n", total_agree);

//     return total_agree;
// }




int main(void) {
    
    /*CSVConversion Script*/
    const char *input = "/home/kjw2kor/Backups/Raw_Script/log/HA1_M0.log";
    const char *output = "/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv";

    if (convert_log_to_csv(input, output) == 0) {
        printf("Conversion completed successfully.\n");
    } else {
        printf("Conversion failed.\n");
    }

    /*Orientation Script*/
    process_csv_file("/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv");

    /*NEEDS DEBUGGING*/

    // /*Normalization*/
    // const char *csv_path = "/home/kjw2kor/shared_folder/csv/HA1_M0.csv";
    // normalize_csv(csv_path);

    /*CalibrationMetrices-Ka,Ma,ba*/
    const char *basepath = "/home/kjw2kor/Backups/Raw_Script/csv/";
    run_calibration(basepath);
    
    /*Calibration*/
    double X = 1.00;
    double Y = 1.00;
    double Z = 1.00;

    double Ka[3][3] = {
        {X, 0, 0},
        {0, Y, 0},
        {0, 0, Z}
    };
    apply_calibration("/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv", Sa, Ma, ba);

    /*Calibration_Initial*/
    const char *input_file = "/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv";
    const char *output_file = "/home/kjw2kor/Backups/Raw_Script/csv/CalibrationInitial.csv";

    IMUData *data = NULL;
    int row_count = read_IMU_data(input_file, &data);
    if (row_count <= 0) {
        fprintf(stderr, "Error: Failed to read IMU data.\n");
        return 1;
    }

    printf("Read %d rows from %s\n", row_count, input_file);

    float acc_bias[3], gyro_bias[3];
    compute_bias(data, row_count, acc_bias, gyro_bias);
    apply_calibration2(data, row_count, acc_bias, gyro_bias);

    if (write_calibrated_data(output_file, data, row_count)) {
        printf("Calibration complete. Output written to %s\n", output_file);
    }

    free(data);

    /*Mass CSVConverter*/
    int N;
    printf("Enter number of log files to convert (N): ");
    if (scanf("%d", &N) != 1 || N <= 0) {
        printf("Invalid input.\n");
        return 1;
    }

    // 🔹 Define input and output directories separately
    const char *input_dir  = "/home/kjw2kor/Backups/Raw_Script/log/Calibration_n/";     // Path where .log files are stored
    const char *output_dir = "/home/kjw2kor/Backups/Raw_Script/csv/Calibration_n/"; // Path where .csv files will be saved

    // Create output directory if it doesn’t exist
    char mkdir_cmd[512];
    snprintf(mkdir_cmd, sizeof(mkdir_cmd), "mkdir -p %s", output_dir);
    system(mkdir_cmd);

    // 🔹 Loop through and process all N files
    for (int i = 1; i <= N; i++) {
        char input_path[256], output_path[256];
        snprintf(input_path, sizeof(input_path), "%s%d.log", input_dir, i);
        snprintf(output_path, sizeof(output_path), "%s%d.csv", output_dir, i);

        MassConverter(input_path, output_path);
    }

    printf("\nConversion complete for %d files.\n", N);
    printf("CSV files saved in: %s\n", output_dir);

    /*AutoCalibration*/ 
    if(run_autocalibration("/home/kjw2kor/Backups/Raw_Script/csv/Calibration_n/", 1, 96, 200, 0.02) != 0){
        printf("Calibration failed!\n");
        return 1;
    }
    // const char *input_csv="/home/kjw2kor/shared_folder/csv/HA1_M0.csv";
    // const char *output_csv="/home/kjw2kor/shared_folder/csv/Calibration2.csv";
    apply_calibration_to_csv("/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv","/home/kjw2kor/Backups/Raw_Script/csv/AutoCalibration.csv");

    /*BUTTERWORTH LPF*/
    //  // Parameters
    // const char *input_csv = "/home/kjw2kor/shared_folder/csv/CalibrationInitial.csv";
    // const char *output_csv = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // double fs = 100.0;   // Sampling frequency
    // double fc = 10.0;    // Cutoff frequency
    // int order = 4;       // Filter order

    // // Design filter
    // ButterworthLPF filter;
    // butterworth_coefficients(&filter, order, fs, fc);

    // printf("Filter coefficients:\n");
    // for (int i = 0; i <= filter.order; i++) {
    //     printf("b[%d]=%f  a[%d]=%f\n", i, filter.b[i], i, filter.a[i]);
    // }

    // // Filter CSV
    // filter_csv(input_csv, output_csv, &filter);

    // if (argc < 6) {
    //     printf("Usage: %s <input.csv> <output.csv> <order> <fs> <fc>\n", argv[0]);
    //     return 1;
    // }

    // const char *input_csv = "/home/kjw2kor/shared_folder/csv/CalibrationInitial.csv";
    // const char *output_csv = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // int order = 4;
    // double fs = 100;
    // double fc = 10;

    // if (order < 1 || order > MAX_ORDER) {
    //     fprintf(stderr, "Order must be between 1 and %d.\n", MAX_ORDER);
    //     return 1;
    // }

    // ButterworthLPF filter;
    // butterworth_coefficients(&filter, order, fs, fc);

    // filter_csv(input, output, &filter);
    // int order = 4;
    // double fc = 50.0;
    // double fs = 500.0;
    // double a[COEFF_COUNT], b[COEFF_COUNT];

    // butterworth_coefficent(order, fc, fs, a, b);
    // const char *x= "/home/kjw2kor/shared_folder/csv/CalibrationInitial.csv";
    // const char *y= "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // // int order = 4;

    // process_csv(x, y);


//  #define MAX_ROWS 10000

// const char *input_path  = "/home/kjw2kor/shared_folder/csv/CalibrationInitial.csv";
//     const char *output_path = "/home/kjw2kor/shared_folder/csv/Filtered.csv";

//     char time[MAX_ROWS][32];
//     double sensor_data[MAX_ROWS][MAX_COLS];
//     int num_cols = 0;

//     int rows = read_multicol_csv(input_path, time, sensor_data, MAX_ROWS, &num_cols);
//     if (rows <= 0) {
//         printf("Error: could not read input file.\n");
//         return 1;
//     }

//     printf("Read %d rows and %d columns\n", rows, num_cols);

//     // Example: Filter Ax column (index 0)
//     IIRFilter filt;
//     butterworth_init(&filt);
//     for (int i = 0; i < rows; i++)
//         sensor_data[i][0] = butterworth_process(&filt,sensor_data[i][0]);

//     if (write_multicol_csv(output_path, time, sensor_data, rows, num_cols) != 0) {
//         printf("Error writing output file.\n");
//         return 1;
//     }


// #define INPUT_PATH  "/home/kjw2kor/shared_folder/csv/CalibrationInitial.csv"
// #define OUTPUT_PATH "/home/kjw2kor/shared_folder/csv/Filtered.csv"
//     char time[MAX_ROWS][32];
//     double sensor_data[MAX_ROWS][MAX_COLS];
//     int num_cols = 0;

//     int num_rows = read_multicol_csv(INPUT_PATH, time, sensor_data, MAX_ROWS, &num_cols);
//     if (num_rows <= 0) {
//         printf("Error: Failed to read input CSV.\n");
//         return 1;
//     }

//     // Initialize filters for each numeric column
//     IIRFilter filters[num_cols];
//     for (int i = 0; i < num_cols; i++)
//         butterworth_init(&filters[i]);

//     // Apply filter to each numeric column
//     for (int j = 0; j < num_cols; j++) {
//         for (int i = 0; i < num_rows; i++) {
//             sensor_data[i][j] = butterworth_process(&filters[j], sensor_data[i][j]);
//         }
//     }

//     // Write filtered data
//     if (write_multicol_csv(OUTPUT_PATH, time, sensor_data, num_rows, num_cols) == 0)
//         printf("Filtered data written to: %s\n", OUTPUT_PATH);




// main.c




    // if (argc < 3) {
    //     printf("Usage: %s <input.csv> <output.csv> [window=100] [energy=0.95]\n", argv[0]);
    //     return 1;
    // }

    const char *input1 = "/home/kjw2kor/Backups/Raw_Script/csv/CalibrationInitial.csv";
    const char *output1 = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";
    int M = 100;
    double lpf_energy = 0.95;

    int result = run_auto_butterworth(input1, output1, M, lpf_energy);
    if (result == 0)
        printf("Filtered CSV saved to %s \n", output1);
    else
        printf("Error: processing failed (code %d)\n", result);

    // const char *inputfile = "/home/kjw2kor/shared_folder/csv/HA1_M0.csv";
    // const char *outputfile = "/home/kjw2kor/shared_folder/csv/Filtered2.csv";
    
    // // Simply call the main processing function
    // run_InEKF(inputfile, outputfile);

        printf("Starting InEKF debug...\n");
    
    // Test 1: Basic initialization
    printf("Test 1: Initializing filter...\n");
    InEKF filter;
    inekf_init_default(&filter);
    printf("Filter initialization successful\n");
    
    // Test 2: Single prediction step
    printf("Test 2: Testing prediction...\n");
    double gyro[3] = {0.1, 0.2, 0.3};
    double accel[3] = {0.0, 0.0, 9.8};
    double dt = 0.01;
    
    inekf_predict(&filter, gyro, accel, dt);
    printf("Prediction successful\n");
    
    // Test 3: ZUPT update
    printf("Test 3: Testing ZUPT update...\n");
    double zero_vel[3] = {0.0, 0.0, 0.0};
    inekf_zupt_update(&filter, zero_vel);
    printf("ZUPT update successful\n");
    
    printf("All tests passed!\n");
    const char *INPUT_CSV  = "/home/kjw2kor/Backups/Raw_Script/csv/HA1_M0.csv";
    const char *OUTPUT_CSV = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered2.csv";
    // const char *LOG_FILE   = "/home/kjw2kor/Raw_Script/log/Test.log";

    // int rc = run_InEKF(INPUT_CSV, OUTPUT_CSV, LOG_FILE);

    int rc = run_InEKF(INPUT_CSV, OUTPUT_CSV);
    if (rc != 0) {
        fprintf(stderr, "run_InEKF failed with code %d\n", rc);
        return rc;
    }
    // fprintf(stderr, "Filtering completed. Output: %s, Log: %s\n", OUTPUT_CSV, LOG_FILE);

    fprintf(stderr, "Filtering completed. Output: %s\n",OUTPUT_CSV);




    // const char* x = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // const char* y = "/home/kjw2kor/shared_folder/csv/Jerk.csv";

    // // Allocate memory for data
    // DataPoint* data1 = malloc(MAX_ROWS * sizeof(DataPoint));
    // Derivative* derivatives = malloc(MAX_ROWS * sizeof(Derivative));

    // if (data1 == NULL || derivatives == NULL) {
    //     printf("Error: Memory allocation failed\n");
    //     return 1;
    // }

    // // Read CSV file
    // int data_count = readCSV(x, data1);
    // if (data_count <= 0) {
    //     printf("Error: No data read from file\n");
    //     free(data1);
    //     free(derivatives);
    //     return 1;
    // }

    // printf("Read %d data points from %s\n", data_count, input_file);
    
    // // Analyze time data
    // analyzeTimeData(data1, data_count);

    // // Calculate derivatives
    // calculateDerivatives(data1, data_count, derivatives);

    // // Write results to CSV
    // writeDerivativesCSV(y, derivatives, data_count);

    // // Free memory
    // free(data1);
    // free(derivatives);

    const char* input_file2 = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";
    const char* output_derivatives_file = "/home/kjw2kor/Backups/Raw_Script/csv/Jerk.csv";
    const char* output_double_derivatives_file = "/home/kjw2kor/Backups/Raw_Script/csv/Snap.csv";

    // // Generate double derivatives filename
    // char output_double_derivatives_file[MAX_LINE_LENGTH];
    // snprintf(output_double_derivatives_file, sizeof(output_double_derivatives_file), 
    //          "double_%s", output_derivatives_file);

    // Allocate memory for data
    DataPoint* data1 = malloc(MAX_ROWS * sizeof(DataPoint));
    Derivative* derivatives = malloc(MAX_ROWS * sizeof(Derivative));
    DoubleDerivative* double_derivatives = malloc(MAX_ROWS * sizeof(DoubleDerivative));

    if (data1 == NULL || derivatives == NULL || double_derivatives == NULL) {
        printf("Error: Memory allocation failed\n");
        return 1;
    }

    // Read CSV file
    int data_count = readCSV(input_file2, data1);
    if (data_count <= 0) {
        printf("Error: No data read from file\n");
        free(data1);
        free(derivatives);
        free(double_derivatives);
        return 1;
    }

    printf("Read %d data points from %s\n", data_count, input_file2);
    
    // Analyze time data
    analyzeTimeData(data1, data_count);

    // Calculate first derivatives
    calculateDerivatives(data1, data_count, derivatives);

    // Calculate double derivatives
    calculateDoubleDerivatives(data1, derivatives, data_count, double_derivatives);

    // Write results to CSV files
    writeDerivativesCSV(output_derivatives_file, derivatives, data_count);
    writeDoubleDerivativesCSV(output_double_derivatives_file, double_derivatives, data_count);

    // Free memory
    free(data1);
    free(derivatives);
    free(double_derivatives);

    printf("Derivative calculations completed successfully.\n");
    printf("First derivatives saved to: %s\n", output_derivatives_file);
    printf("Double derivatives saved to: %s\n", output_double_derivatives_file);


    /*DYNAMIC THRESHOLD*/
    // IMUSample dataX[MAX_SAMPLES];
    // float accMag[MAX_SAMPLES];
    // AdaptiveThreshold th;

    // int numSamples = readIMUCSV("/home/kjw2kor/shared_folder/csv/Filtered.csv", dataX, MAX_SAMPLES);
    // if (numSamples <= 0) {
    //     printf("No data read or file error.\n");
    //     return 1;
    // }

    // computeMagnitude(dataX, accMag, numSamples);
    // initThreshold(&th, INIT_THRESHOLD_POS, INIT_THRESHOLD_NEG, TIME_THRESHOLD);

    // printf("Processing %d samples...\n", numSamples);

    // for (int i = 1; i < numSamples - 1; i++) {
    //     if (i > WINDOW_SIZE && i % WINDOW_SIZE == 0)
    //         updateAdaptiveThreshold(&th, &accMag[i - WINDOW_SIZE], WINDOW_SIZE);

    //     if (detectPositiveSpike(accMag, i, &th))
    //         printf("Positive spike at %s (%.3f m/s²)\n", dataX[i].timeStr, accMag[i]);
    //         addThresholdDetection("Positive", dataX[i].timeStr, accMag[i]);

    //     if (detectNegativeSpike(accMag, i, &th))
    //         addThresholdDetection("Negative", dataX[i].timeStr, accMag[i]);

    // }

    // printf("\nFinal thresholds: pos=%.3f neg=%.3f\n",
    //        th.posThreshold, th.negThreshold);
    
    const char *path3 = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";
    run_Threshold(path3);

           /*DERIVATIVE DETECTION*/
    //  IMUSampleD dataD[DERIVATIVE_MAX_SAMPLES];
    // float accMagD[DERIVATIVE_MAX_SAMPLES];
    // float derivD[DERIVATIVE_MAX_SAMPLES];

    // const char *filename = "/home/kjw2kor/shared_folder/csv/Filtered.csv";  // path to your CSV file
    // int numSamplesD = readIMUCSV_D(filename, dataD, DERIVATIVE_MAX_SAMPLES);

    // if (numSamplesD <= 0) {
    //     printf("Error: No data or failed to read %s\n", filename);
    //     return 1;
    // }

    // computeAccMagnitude_D(dataD, accMagD, numSamplesD);
    // computeDerivative_D(accMagD, derivD, numSamplesD);
    // detectDerivativeSpikes_D(dataD, accMagD, derivD, numSamplesD);
//      const char *inputFile = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
// DerivSample *derivData = malloc(sizeof(DerivSample) * MAX_SAMPLES_DERIV);
// if (!derivData) {
//     perror("Memory allocation failed");
//     return 1;
// }

// int n = loadCSV(inputFile, derivData);
// if (n <= 0) {
//     printf("No valid data found or failed to load file: %s\n", inputFile);
//     free(derivData);
//     return 1;
// }

// detectSpikes(derivData, n);

// free(derivData);

    const char *inputFile = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";

    DerivSample *samples = malloc(sizeof(DerivSample) * MAX_SAMPLES_DERIV);
    if (!samples) { perror("malloc"); return 1; }

    int n = loadCSV(inputFile, samples);
    if (n <= 0) {
        printf("Failed to read data from %s\n", inputFile);
        free(samples);
        return 1;
    }

    detectSpikes(samples, n);

    free(samples); 

    /*ENERGY ENVELOPE*/
    
     const char *env_input = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";

    EnvelopeSample *env_data = malloc(sizeof(EnvelopeSample) * MAX_SAMPLES_ENVELOPE);
    if (!env_data) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    int n_env = loadEnvelopeCSV(env_input, env_data);
    if (n_env <= 0) {
        fprintf(stderr, "Error loading CSV file\n");
        free(env_data);
        return 1;
    }

    printf("Read %d samples for Energy Envelope analysis\n", n_env);
    runImprovedEnvelopeSliding(env_data, n_env, FS);

    free(env_data);
    /* CUSUM-Poisson Distribution */
    // char **times = NULL;
    // double *ax_data = NULL;

    // // Load CSV data from the file
    // int num_samples = load_csv("/home/kjw2kor/shared_folder/csv/Filtered.csv", &times, &ax_data);
    // if (num_samples == 0) {
    //     // If no samples were loaded, print an error and return
    //     printf("Error loading data\n");
    //     return 1;
    // }

    // printf("Successfully loaded %d samples\n", num_samples);

    // // Call the function to detect spikes using CUSUM
    // detect_spikes(ax_data, times, num_samples);

    // // Free dynamically allocated memory
    // if (times != NULL) {
    //     for (int i = 0; i < num_samples; i++) {
    //         free(times[i]);  // Free each allocated time string
    //     }
    //     free(times);  // Free the array of time string pointers
    // }

    // if (ax_data != NULL) {
    //     free(ax_data);  // Free the array of accelerometer data
    // }


    const char *path = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";
    run_cusum_on_csv(path);

     
    /*CWT DETECTOR*/
  const char *path1 = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";  
run_cwt_to_csv(path1);


/*POLLING*/
    const char *csv = "/home/kjw2kor/Backups/Raw_Script/csv/Filtered.csv";

    // CUSUM(csv);
    // EnergyEnvelope(csv);
    // CWT(csv);
    // DerivativeDetect(csv);
    // Threshold(csv);

    printAllDetections(); 

    fuseDetections();
    // const char *filename = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // SpikeEvent results[10000];
    // int o = run_polling_detection(filename, results, 10000);

    //     printf("Starting main()...\n");  // <-- Debug check

    // const char *filename = "/home/kjw2kor/shared_folder/csv/Filtered.csv";
    // SpikeEvent results[1000];

    // int o = run_polling_detection(filename, results, 1000);

    // printf("\nPolling detection finished. Total %d consensus spikes.\n", o);
    // fflush(stdout);  // force print immediately

    /*plotter.py*/
      int ret = system("python3 ../scripts/plotter.py");  //Add "." in order to get into outer folder
    if (ret == -1) {
        perror("Failed to run plotter.py");
    }

    /*plotter2.py*/
      int ret2 = system("python3 ../scripts/plotter2.py");  //Add "." in order to get into outer folder
    if (ret2 == -1) {
        perror("Failed to run plotter2.py");
    }

    /*plotter3.py*/
      int ret3 = system("python3 ../scripts/plotter3.py");  //Add "." in order to get into outer folder
    if (ret3 == -1) {
        perror("Failed to run plotter3.py");
    }

    /*plotter4.py*/
      int ret4 = system("python3 ../scripts/plotter4.py");  //Add "." in order to get into outer folder
    if (ret4 == -1) {
        perror("Failed to run plotter4.py");
    }

    

     /*plotter5.py*/
      int ret5 = system("python3 ../scripts/plotter5.py");  //Add "." in order to get into outer folder
    if (ret5 == -1) {
        perror("Failed to run plotter5.py");
    }

       /*plotter6.py*/
      int ret6 = system("python3 ../scripts/plotter6.py");  //Add "." in order to get into outer folder
    if (ret6 == -1) {
        perror("Failed to run plotter6.py");
    }
    
           /*plotter7.py*/
      int ret7 = system("python3 ../scripts/plotter7.py");  //Add "." in order to get into outer folder
    if (ret7 == -1) {
        perror("Failed to run plotter7.py");
    }
    
           /*plotter8.py*/
      int ret8 = system("python3 ../scripts/plotter8.py");  //Add "." in order to get into outer folder
    if (ret8 == -1) {
        perror("Failed to run plotter8.py");
    }
    
    
    return 0;
}
