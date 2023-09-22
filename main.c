#include "headers.h"
#define RGB 3

int main(int argc, char **argv) {
    // declare all globals
    BMP *bmp, *local_bmp, *final_local_bmp, *final_bmp;
    clock_t start, finish;
    UCHAR *flat_img, *flat_frag_img, *final_flat_frag_img, *final_flat_img;
    int nproc, rank, height, width, local_height, size, items, from, to;
    int *counts, *allc, *displs;
    double std_dev, origin, kernel_dim, kernel_max, colour_max;
    double **kernel;

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
    //origin = KERNEL_DIMENSION_SD * std_dev;
    origin = (kernel_dim - 1) / 2;

    // allocate and init kernel
    if ((kernel = (double **)malloc(kernel_dim * sizeof(double *))) == NULL) {
        fprintf(stderr, "Malloc error: kernel\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    for (int i = 0; i < kernel_dim; i++) {
        if ((kernel[i] = (double *)malloc(kernel_dim * sizeof(double))) == NULL) {
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

        // allocate vectors for counts, allocations, and displacements
        if ((counts = (int *)malloc(nproc * sizeof(int))) == NULL) {
            fprintf(stderr, "Malloc error: counts\n");
            free_kernels(kernel, kernel_dim);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if ((allc = (int *)malloc(nproc * sizeof(int))) == NULL) {
            fprintf(stderr, "Malloc error: allc\n");
            free_kernels(kernel, kernel_dim);
            free(counts);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if ((displs = (int *)malloc(nproc * sizeof(int))) == NULL) {
            fprintf(stderr, "Malloc error: displs\n");
            free_kernels(kernel, kernel_dim);
            free(counts);
            free(allc);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // allocate rows to processes
        assign_rows(nproc, height, width, std_dev, counts, allc, displs);

        // allocate vector to send image fragments
        size = counts[nproc - 1] + displs[nproc - 1];
        if ((flat_img = (UCHAR *)malloc(size * sizeof(UCHAR))) == NULL) {
            fprintf(stderr, "Malloc error: flat_img\n");
            free_kernels(kernel, kernel_dim);
            free(counts);
            free(allc);
            free(displs);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // flatten image to 1D vector for scattering
        int offset = 0;
        for (int i = 0; i < nproc; i++) {
            from = (i == 0 || offset - 3 * std_dev < 0 ? 0 : offset - 3 * std_dev);
            to = from + counts[i] / width / RGB;
            flatten_image(bmp, flat_img, from, to, width, displs[i]);
            offset += allc[i];
        }
        BMP_Free(bmp);
    }

    // broadcast main image height and width data
    MPI_Bcast(&height, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&width, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    // scatter process pixel counts
    MPI_Scatter(counts, 1, MPI_INT, &items, 1, MPI_INT, ROOT,
                MPI_COMM_WORLD);
    
    // allocate local image fragment vectors
    if ((flat_frag_img = (UCHAR *)malloc(items * sizeof(UCHAR))) == NULL) {
        fprintf(stderr, "Malloc error: flat_frag_img\n");
        free_kernels(kernel, kernel_dim);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // scatter the pixel data
    MPI_Scatterv(flat_img, counts, displs, MPI_UNSIGNED_CHAR,
                    flat_frag_img, items, MPI_UNSIGNED_CHAR, ROOT,
                    MPI_COMM_WORLD);
    if (rank == ROOT) free(flat_img);

    // convert image fragment to bmp for convolution
    local_height = items / width / RGB;
    local_bmp = BMP_Create(width, local_height, 24);
    reconstruct_bmp(local_bmp, flat_frag_img, local_height, width);
    free(flat_frag_img);

    // apply convolution to fragment
    final_local_bmp = BMP_Create(width, local_height, 24);
    applyConvolution(kernel, kernel_dim, origin, colour_max, local_bmp,
                        final_local_bmp);
    free_kernels(kernel, kernel_dim);
    BMP_Free(local_bmp);

    // determine offsets for overlaps
    int rows_per_proc = height / nproc;
    int remainder = height % nproc;
    if (nproc == 1) {
        from = 0;
        to = local_height;
    } else if (rank == ROOT) {
        from = 0;
        to = rows_per_proc + (rank < remainder ? 1 : 0);
    } else if (rank == nproc - 1) {
        from = 3 * std_dev;
        to = local_height;
    } else {
        from = 3 * std_dev;
        to = from + rows_per_proc + (rank < remainder ? 1 : 0);
    }

    // allocate return vector
    items = (to - from) * width * RGB;
    final_flat_frag_img = (UCHAR *)malloc(items * sizeof(UCHAR));
    if (final_flat_frag_img == NULL) {
        fprintf(stderr, "Malloc error: final_flat_frag_img\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // flatten the return image without overlaps
    flatten_image(final_local_bmp, final_flat_frag_img, from, to, width, 0);
    BMP_Free(final_local_bmp);

    // receive convoluted fragments and put the image back together
    if (rank == ROOT) {
        // reinitialise displs
        int offset = 0;
        for (int i = 0; i < nproc; i++) {
            counts[i] = allc[i] * width * RGB;
            displs[i] = offset;
            offset += counts[i];
        }

        size = height * width * RGB;
        if ((final_flat_img = (UCHAR *)malloc(size * sizeof(UCHAR))) == NULL) {
            fprintf(stderr, "Malloc error: final_flat_img\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // create final bmp image
        final_bmp = BMP_Create(width, height, 24);
    }

    // gather the image fragments
    if (rank < nproc) {
        MPI_Gatherv(final_flat_frag_img, items, MPI_UNSIGNED_CHAR,
                    final_flat_img, counts, displs, MPI_UNSIGNED_CHAR, ROOT,
                    MPI_COMM_WORLD);
        free(final_flat_frag_img);
        if (rank == ROOT) {
            free(counts);
            free(allc);
            free(displs);
        }
    }

    // reconstruct the final image and save it
    if (rank == ROOT) {
        reconstruct_bmp(final_bmp, final_flat_img, height, width);
        free(final_flat_img);
        BMP_WriteFile(final_bmp, argv[2]);
        BMP_Free(final_bmp);
        finish = clock();
        printf("Time taken: %f seconds\n",
            (float)(finish - start) / CLOCKS_PER_SEC);
    }

    MPI_Finalize();
}