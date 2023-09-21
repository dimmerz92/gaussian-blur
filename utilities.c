#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "qdbmp.h"

#define EXT ".bmp"

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
