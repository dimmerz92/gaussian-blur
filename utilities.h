int parse_args(int argc, char *argv[], float *std_dev);
void free_kernels(float **kernel, float kernel_dim);
void assign_rows(int nproc, int height, int width, int *counts, int *allc,
                 int *displs) ;