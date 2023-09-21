int parse_args(int argc, char *argv[], float *std_dev);
void free_kernels(float **kernel, float kernel_dim);
void assign_rows(int nproc, int height, int width, int *counts, int *allc,
                 int *displs);
void flatten_image(BMP *bmp, UCHAR *flattened_img, int height, int width,
                   int top_offset, int bottom_offset);