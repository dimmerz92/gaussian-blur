#include "headers.h"

int main(int argc, char **argv) {
    // declare all globals
    BMP *bmp, *local_bmp, *final_local_bmp, *final_bmp;
    clock_t start, finish;
    UCHAR *flat_img, *flat_frag_img, *final_flat_frag_img, *final_flat_img;
    int nproc, rank, height, width, local_height, size, items;
    int *counts, *allc, *displs;
    float std_dev, origin, kernel_dim, kernel_max, colour_max;
    float *kernel;

    // init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    // parse args
    if (parse_args(argc, argv, &std_dev) < 0) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // init kernal dim and origin
    kernel_dim = 2 * KERNEL_DIMENSION_SD * std_dev + 1;
    origin = KERNEL_DIMENSION_SD * std_dev;

    // allocate and init kernel
    if ((kernel = (float **)malloc(kernel_dim * sizeof(float *))) == NULL) {
        fprintf(stderr, "Malloc error: kernel\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < kernel_dim; i++) {
        if ((kernel[i] = (float *)malloc(kernel_dim * sizeof(float))) == NULL) {
            fprintf(stderr, "Malloc error: kernel[%d]\n", i);
            free(kernel);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    generateGaussianKernel(kernel, kernel_dim, std_dev, origin, &kernel_max,
                           &colour_max);

    // ROOT read and distribute image fragments
    if (rank == ROOT) {
        // start the clock
        start = clock();

        // read the input BMP and its metadata
        bmp = BMP_ReadFile(argv[1]);
        BMP_CHECK_ERROR(stderr, -1);
        width = BMP_GetWidth(bmp);
        height = BMP_GetHeight(bmp);
        
    }
}