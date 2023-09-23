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
int parse_args(int argc, char *argv[], double *std_dev) {
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
void free_kernels(double **kernel, double kernel_dim) {
    for (int i = 0; i < kernel_dim; i++) {
        free(kernel[i]);
    }
    free(kernel);
}

/* allocates rows and displacements, accounting for necessary overlaps */
void assign_rows(int nproc, int height, int width, double std_dev, int *counts,
                 int *allc, int *displs) {
    int rows_per_proc = height / nproc;
    int remainder = height % nproc;
    int accumulator = 0;

    for (int i = 0; i < nproc; i++) {
        // allocate base rows excluding overlaps
        allc[i] = rows_per_proc + (i < remainder ? 1 : 0);

        // allocate counts including overlaps
        counts[i] = (allc[i] + (nproc == 1 ? 0 : (i == 0 || i == nproc - 1 ?
            3 * std_dev: 6 * std_dev))) * width * RGB;

        // allocate displacements
        displs[i] = accumulator;

        accumulator += counts[i];
    }
}

/* flattens a BMP image to a 1D array, with three positions for every pixel */
/* as the RBG values */
void flatten_image(BMP *bmp, UCHAR *flat_img, int start, int stop,
                   int width, int accum_offset) {
    UCHAR r, g, b;
    int accumulator = accum_offset;

    for (int row = start; row < stop; row++) {
        for (int col = 0; col < width; col++) {
            BMP_GetPixelRGB(bmp, col, row, &r, &g, &b);
            flat_img[accumulator] = r;
            flat_img[accumulator + 1] = g;
            flat_img[accumulator + 2] = b;
            accumulator += 3;
        }
    }
}

/* reconstructs a BMP image from 1D array containing pixel data read in */
/* increments of 3 (RGB) */
void reconstruct_bmp(BMP *bmp, UCHAR *img, int height, int width) {
    UCHAR r, g, b;
    int accumulator = 0;

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            r = round(img[accumulator]);
            g = round(img[accumulator + 1]);
            b = round(img[accumulator + 2]);
            BMP_SetPixelRGB(bmp, col, row, r, g, b);
            accumulator += 3;
        }
    }
}