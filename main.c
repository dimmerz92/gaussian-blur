#include "headers.h"

int main(int argc, char **argv) {
    BMP *bmp, *local_bmp, *final_local_bmp, *final_bmp;
    clock_t start, finish;
    UCHAR *flat_img, *flat_frag_img, *final_flat_frag_img, *final_flat_img;
    int nproc, rank, height, width, size, items;
    int *counts, *allc, *displs;
    float std_dev, origin, kernel_dim, kernel_max, colour_max;
    float *kernel;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (parse_args(argc, argv, &std_dev) < 0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}