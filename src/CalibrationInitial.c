#include "CalibrationInitial.h"

// Compute magnitude of a 3D vector
float compute_magnitude(float x, float y, float z) {
    return sqrtf(x * x + y * y + z * z);
}

// Count total data rows (excluding header)
int count_rows(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file for counting");
        return -1;
    }

    int count = 0;
    char line[LINE_BUFFER_SIZE];

    fgets(line, LINE_BUFFER_SIZE, file); // skip header

    while (fgets(line, LINE_BUFFER_SIZE, file)) {
        if (strlen(line) > 5)
            count++;
    }

    fclose(file);
    return count;
}

// Read IMU data dynamically
int read_IMU_data(const char *filename, IMUData **data) {
    int total_rows = count_rows(filename);
    if (total_rows <= 0)
        return -1;

    *data = (IMUData *)malloc(total_rows * sizeof(IMUData));
    if (!*data) {
        perror("Memory allocation failed");
        return -1;
    }

    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        free(*data);
        return -1;
    }

    char line[LINE_BUFFER_SIZE];
    fgets(line, LINE_BUFFER_SIZE, file); // skip header

    int row_count = 0;
    while (fgets(line, LINE_BUFFER_SIZE, file) && row_count < total_rows) {
        char time_str[32];
        float Ax, Ay, Az, Gx, Gy, Gz;

        // Parse the line
        if (sscanf(line, "%31[^,],%f,%f,%f,%f,%f,%f",
                   time_str, &Ax, &Ay, &Az, &Gx, &Gy, &Gz) == 7) {

            IMUData temp;
            strcpy(temp.time, time_str);
            temp.Ax = Ax;
            temp.Ay = Ay;
            temp.Az = Az;
            temp.Gx = Gx;
            temp.Gy = Gy;
            temp.Gz = Gz;
            temp.is_stationary = (row_count < CALIBRATION_SAMPLES) ? 1 : 0; // ✅ First N samples stationary

            (*data)[row_count++] = temp;
        }
    }

    fclose(file);
    return row_count;
}

// Compute biases from first N stationary samples
void compute_bias(const IMUData *data, int row_count, float acc_bias[3], float gyro_bias[3]) {
    int stationary_count = (row_count < CALIBRATION_SAMPLES) ? row_count : CALIBRATION_SAMPLES;

    acc_bias[0] = acc_bias[1] = acc_bias[2] = 0;
    gyro_bias[0] = gyro_bias[1] = gyro_bias[2] = 0;

    for (int i = 0; i < stationary_count; i++) {
        acc_bias[0] += data[i].Ax;
        acc_bias[1] += data[i].Ay;
        acc_bias[2] += data[i].Az;

        gyro_bias[0] += data[i].Gx;
        gyro_bias[1] += data[i].Gy;
        gyro_bias[2] += data[i].Gz;
    }

    for (int i = 0; i < 3; i++) {
        acc_bias[i] /= stationary_count;
        gyro_bias[i] /= stationary_count;
    }

    // Adjust for gravity if in m/s² range
    if (fabsf(acc_bias[2]) > 5.0f)
        acc_bias[2] -= 9.81f;

    printf("Using first %d samples for calibration.\n", stationary_count);
    printf("ACC Bias: %.4f, %.4f, %.4f\n", acc_bias[0], acc_bias[1], acc_bias[2]);
    printf("GYRO Bias: %.4f, %.4f, %.4f\n", gyro_bias[0], gyro_bias[1], gyro_bias[2]);
}

// Apply calibration (bias removal)
void apply_calibration2(IMUData *data, int row_count, const float acc_bias[3], const float gyro_bias[3]) {
    for (int i = 0; i < row_count; i++) {
        data[i].Ax -= acc_bias[0];
        data[i].Ay -= acc_bias[1];
        data[i].Az -= acc_bias[2];

        data[i].Gx -= gyro_bias[0];
        data[i].Gy -= gyro_bias[1];
        data[i].Gz -= gyro_bias[2];
    }
}

// Write calibrated data to output CSV
int write_calibrated_data(const char *output_filename, const IMUData *data, int row_count) {
    FILE *out = fopen(output_filename, "w");
    if (!out) {
        perror("Failed to open output file");
        return 0;
    }

    fprintf(out, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");
    for (int i = 0; i < row_count; i++) {
        fprintf(out, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                data[i].time,
                data[i].Ax, data[i].Ay, data[i].Az,
                data[i].Gx, data[i].Gy, data[i].Gz );
        //data[i].is_stationary ? "Stationary" : "Moving"
    }

    fclose(out);
    return 1;
}
