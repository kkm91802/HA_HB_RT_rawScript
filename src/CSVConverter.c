#include "CSVConverter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Convert signed binary (up to 32 bits) from int value
static int signed_bin_to_dec(uint32_t value, int bit_len) {
    if (value & (1 << (bit_len - 1))) {
        return (int)value - (1 << bit_len);
    } else {
        return (int)value;
    }
}

// Convert hex string (2 chars) to integer
static int hex_to_int(const char *hex) {
    int val;
    sscanf(hex, "%x", &val);
    return val;
}

// Main conversion logic
int convert_log_to_csv(const char *input_path, const char *output_path) {
    FILE *infile = fopen(input_path, "r");
    FILE *outfile = fopen(output_path, "w");

    if (!infile || !outfile) {
        printf("Error opening file: %s or %s\n", input_path, output_path);
        return 1;
    }

    fprintf(outfile, "Time,Ax,Ay,Az,Gx,Gy,Gz\n");

    char line[LINE_SIZE];
    double latest_accel[3] = {0}, latest_gyro[3] = {0};
    int accel_ready = 0, gyro_ready = 0;
    char latest_time[64];

    while (fgets(line, sizeof(line), infile)) {
        char *parts[20];
        int count = 0;

        char *token = strtok(line, " \t\n");
        while (token && count < 20) {
            parts[count++] = token;
            token = strtok(NULL, " \t\n");
        }
        if (count < 14) continue;

        char *timestamp = parts[0];
        char *can_id = parts[3];

        char binary[65] = {0};
        for (int i = 6; i < 14; i++) {
            int val = hex_to_int(parts[i]);
            for (int b = 7; b >= 0; b--) {
                strcat(binary, (val & (1 << b)) ? "1" : "0");
            }
        }

        if (strcmp(can_id, "0x2DA") == 0) {
            int start = 0;
            for (int i = 0; i < 3; i++) {
                char buf[20];
                strncpy(buf, binary + start, 16);
                buf[16] = '\0';
                int val = (int)strtol(buf, NULL, 2);
                val = signed_bin_to_dec(val, 16);
                latest_accel[i] = val / 100.0;
                start += 16;
            }
            latest_accel[0] = -latest_accel[0];
            latest_accel[2] = -latest_accel[2];

            strncpy(latest_time, timestamp, strlen(timestamp) - 1);
            latest_time[strlen(timestamp) - 1] = '\0';
            accel_ready = 1;
        }
        else if (strcmp(can_id, "0x2DB") == 0) {
            int start = 0;
            for (int i = 0; i < 3; i++) {
                char buf[25];
                strncpy(buf, binary + start, 19);
                buf[19] = '\0';
                int val = (int)strtol(buf, NULL, 2);
                val = signed_bin_to_dec(val, 19);
                latest_gyro[i] = val / 100.0;
                start += 19;
            }
            gyro_ready = 1;
        }

        if (accel_ready && gyro_ready) {
            fprintf(outfile, "%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
                latest_time,
                latest_accel[0], latest_accel[1], latest_accel[2],
                latest_gyro[0], latest_gyro[1], latest_gyro[2]);
            accel_ready = gyro_ready = 0;
        }
    }

    fclose(infile);
    fclose(outfile);

    printf("Converted %s → %s successfully.\n", input_path, output_path);
    return 0;
}

