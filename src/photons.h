#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>

#define CENTROIDS_CENT_PIXEL                  0x8000
#define CENTROIDS_TABLE_COLS                  8

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
    int box_t;
    int box_n;
    int pixel_photon_num;
    int overlap_max;
    double sum_min;
    double sum_max;
    DataType threshold; 
};

/* ----------------------------------------------------------------------------*/
/**
 * \brief 
 *
 * @tparam DataType
 */
/* ----------------------------------------------------------------------------*/
template <typename DataType> struct photons {
    DataType *image;
    uint16_t *out;
    int x;
    int y;
};

enum {
    CENTROIDS_PARAMS_OK = 0,
    CENTROIDS_PARAMS_BAD = 1
};

template <typename DataType> int centroids_initialize_params(centroid_params<DataType> &params);

template<typename DataType> int centroids_process(DataType *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, size_t N, centroid_params<DataType> params);

template<typename DataType> int centroids_process_bias(DataType *pixels, uint16_t *out, double *bias, 
        size_t X, size_t Y);

template<typename DataType> int centroids_process_photons(photons<DataType> *photon_map,
        double *photon_table, int n_photons, centroid_params<DataType> params);

template<typename DataType> int centroids_find_photons(DataType *image, uint16_t *out, double *table,
        DataType **image_list, uint16_t **out_list,
        size_t X, size_t Y, centroid_params<DataType> params);

int _calculate_com(double *pixels, int *x, int *y, double *com_x, double *com_y, int n);
template<typename DT> void _swap(DT *a, DT *b);
void _bubble_sort(double *vals, int *x, int *y, int n);

#endif
