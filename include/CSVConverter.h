#ifndef CSV_CONVERTER_H
#define CSV_CONVERTER_H

#include <stdint.h>

// Buffer size for reading lines
#ifdef LINE_SIZE
#undef LINE_SIZE
#endif
#define LINE_SIZE 256

// Converts a signed binary (up to 32 bits) to a decimal integer
static int signed_bin_to_dec(uint32_t value, int bit_len);

// Converts a 2-character hex string to integer
static int hex_to_int(const char *hex);

// Converts a log file to a CSV file (main utility function)
int convert_log_to_csv(const char *input_path, const char *output_path);

#endif
