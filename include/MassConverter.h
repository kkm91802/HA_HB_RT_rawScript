#ifndef MASS_CONVERTER_H
#define MASS_CONVERTER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef LINE_SIZE
#undef LINE_SIZE
#endif
#define LINE_SIZE 256

// ---- Function Prototypes ----

// Converts a signed binary (up to 32 bits) to decimal
static int signed_bin_to_dec(uint32_t value, int bit_len);

// Converts a 2-character hex string to integer
static int hex_to_int(const char *hex);

// Core conversion function
int MassConverter(const char *input_path, const char *output_path);

#endif // MASS_CONVERTER_H
