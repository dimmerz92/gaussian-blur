int parse_args(int argc, char *argv[], double *std_dev);
void free_kernels(double **kernel, double kernel_dim);
void assign_rows(int nproc, int height, int width, double std_dev, int *counts,
                 int *allc, int *displs);
void flatten_image(BMP *bmp, UCHAR *flat_img, int start, int stop,
                   int width, int accum_offset);
void reconstruct_bmp(BMP *bmp, UCHAR *img, int height, int width);