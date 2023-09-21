#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "qdbmp.h"

#define EXT ".bmp"
#define RGB 3

/* Checks if a given string is an int/float as atoi is unreliable on its own */
int is_positive_number(char *str) {
    int decimal = 0;
    if (!isdigit(str[0]) && str[0] != '.') return -1;
    if (str[0] == '.') decimal++;

    for (int i = 1; str[i] != '\0'; i++) {
        if (str[i] == '.') decimal++;
        if (!isdigit(str[i]) || decimal > 1) {
            if (str[i] == '.' && decimal <= 1) {
                continue;
            } else {
                return -1;
            }
        }
    }

    return 0;
}

/* Parse args to check correct types and formats passed */
int parse_args(int argc, char *argv[], float *std_dev) {
    if (argc != 4) {
        fprintf(stderr,
                "Usage: %s <input_file> <output_file> <std_dev>",
                argv[0]);
        return -1;
    } else if (strcmp(strrchr(argv[1], '.'), EXT) != 0 ||
               strcmp(strrchr(argv[2], '.'), EXT) != 0) {
        fprintf(stderr,
            "<input_file> or <output_file> must be a string ending in '.bmp'");
        return -1;
    }else if (!is_positive_number(argv[3]) && (*std_dev = atoi(argv[3])) < 1) {
        fprintf(stderr,
            "<std_dev> must be an integer >= 1");
        return -1;
    }

    return 0;
}

/* frees the kernel items */
void free_kernels(float **kernel, float kernel_dim) {
    for (int i = 0; i < kernel_dim; i++) {
        free(kernel[i]);
    }
    free(kernel);
}

/* allocates rows and displacements, accounting for necessary overlaps */
void assign_rows(int nproc, int height, int width, int *counts, int *allc,
                 int *displs) {
    int rows_per_proc = height / nproc;
    int remainder = height % nproc;
    int rows, overlap, offset = 0;

    for (int i = 0; i < nproc; i++) {
        rows = rows_per_proc + (i < remainder ? 1 : 0);
        overlap = (nproc > 1 ? (i == 0 || i == nproc - 1 ? 1 : 2) : 0);
        counts[i] = (rows + overlap) * width * RGB;
        allc[i] = rows * width * RGB;
        displs[i] = offset;
        offset += counts[i];
    }
}

/* flattens a BMP image to a 1D vector to be read in increments of 3 (RGB) */
void flatten_image(BMP *bmp, UCHAR *flattened_img, int height, int width,
                   int top_offset, int bottom_offset) {
    UCHAR r, g, b;
    int accumulator = 0;

    for (int row = 0 + top_offset; row < height - bottom_offset; row++) {
        for (int col = 0; col < width; col++) {
            BMP_GetPixelRGB(bmp, col, row, &r, &g, &b);
            flattened_img[accumulator] = r;
            flattened_img[accumulator + 1] = g;
            flattened_img[accumulator + 2] = b;
            accumulator += 3;
        }
    }
}