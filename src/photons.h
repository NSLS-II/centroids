#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>

extern int _debug_print_flag;
extern int _error_print_flag;

#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__);

#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__); 

template<typename DT> void swap(DT *a, DT *b);
void bubble_sort(double *vals, int *x, int *y, int n);

template<typename DataType> int process_bias(DataType *pixels, uint16_t *out, size_t X, size_t Y, 
        double *bias);

template<typename DataType> int process_photons(DataType *image, uint16_t *out, double *table, 
        size_t X, size_t Y, double *bias, int box_val);

template<typename DataType> int find_photons(DataType *image, uint16_t *out, 
        size_t X, size_t Y, DataType threshold, int box_val);

int process_pixel(double *pixels, int *x, int *y, int n, double *com_x, double *com_y, int *sum);

int find_photons_uint16(uint16_t *image, uint16_t *out, double *table, double *bias, 
        size_t X, size_t Y, uint16_t threshold, int box);

#endif
