#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>

extern int _debug_print_flag;
extern int _error_print_flag;

#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__);

#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__); 

template<typename DataType> int process_bias(DataType *pixels, DataType *out, size_t X, size_t Y, 
        double *bias);
template<typename DataType> int process_photons(DataType *image, DataType *out, double *table, 
        size_t X, size_t Y, double *bias);
template<typename DataType> int find_photons(DataType *image, DataType *out, 
        size_t X, size_t Y, DataType threshold);
int process_pixel(double *pixels, int box_val, double *com_x, double *com_y, int *sum);

int find_photons_uint16(uint16_t *image, uint16_t *out, double *table, double *bias, 
        size_t X, size_t Y, uint16_t threshold);
#endif
