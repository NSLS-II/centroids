#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>

#define CENT_PIXEL          0x8000

#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__);

#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__); 

/* ----------------------------------------------------------------------------*/
/**
 * \brief Parameters for controlling the photon counting algorythm
 *
 * \tparam DataType Data image datatype
 */
/* ----------------------------------------------------------------------------*/
template <typename DataType> struct centroid_params {
    int box;
    DataType threshold; 
};

template<typename DT> void swap(DT *a, DT *b);

void bubble_sort(double *vals, int *x, int *y, int n);

template<typename DataType> int process_image(DataType *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, centroid_params<DataType> params);

template<typename DataType> int process_bias(DataType *pixels, uint16_t *out, size_t X, size_t Y, 
        double *bias);

template<typename DataType> int process_photons(DataType *image, uint16_t *out, double *table, 
        size_t X, size_t Y, double *bias, centroid_params<DataType> params);

template<typename DataType> int find_photons(DataType *image, uint16_t *out, 
        size_t X, size_t Y, centroid_params<DataType> params);

int process_pixel(double *pixels, int *x, int *y, int n, double *com_x, double *com_y, int *sum);

#endif
