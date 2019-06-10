#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <memory>

#define CENTROIDS_CENT_PIXEL                  0x8000
#define CENTROIDS_TABLE_COLS                  8

#ifdef DEBUG_OUTPUT
#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__, __VA_ARGS__);
#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, __LINE__, __func__); 
#else
#define DEBUG_PRINT(fmt, ...) \
    do {} while (0)
#define DEBUG_COMMENT(fmt) \
    do {} while (0)
#endif




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
    size_t x;
    size_t y;
    size_t n;
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
    size_t x;
    size_t y;
};

template<typename T>
using PixelValues = std::vector<T>;
template<typename T>
using PixelValuesPtr = std::shared_ptr<PixelValues<T>>;

template <typename OutputType> struct photon_params {
    OutputType x;
    OutputType y;
    OutputType com_x;
    OutputType com_y;
    OutputType sum;
    OutputType bgnd;
    int box_sum;
    PixelValuesPtr<OutputType> values;
};

enum {
    CENTROIDS_PARAMS_OK = 0,
    CENTROIDS_PARAMS_BAD = 1
};

template <typename DataType> 
using PhotonTable = std::vector<photon_params<DataType>>;
template <typename DataType> 
using PhotonTablePtr = std::unique_ptr<PhotonTable<DataType>>;
template <typename DataType> 
using PhotonMap = std::vector<photons<DataType>>;
template <typename DataType> 
using PhotonMapPtr = std::unique_ptr<PhotonMap<DataType>>;

template <typename DataType> void centroids_initialize_params(centroid_params<DataType> &params);
template <typename DataType> int centroids_calculate_params(centroid_params<DataType> &params);

template<typename DataType, typename OutputType>
int centroids_process(DataType *image, uint16_t *out, 
        PhotonTablePtr<OutputType> &photon_table,
        size_t X, size_t Y, size_t N, centroid_params<DataType> params);

template<typename DataType> int centroids_process_bias(DataType *pixels, uint16_t *out, double *bias, 
        size_t X, size_t Y);

template<typename DataType, typename OutputType> 
size_t centroids_process_photons(
        PhotonMapPtr<DataType> &photon_map,
        PhotonTablePtr<OutputType> &photon_table,
        centroid_params<DataType> params);

template<typename DataType> 
size_t centroids_find_photons(DataType *image, uint16_t *out, 
        PhotonMapPtr<DataType> &photon_map,
        size_t X, size_t Y, centroid_params<DataType> params);

template <typename DT>
int _calculate_com(std::unique_ptr <DT[]> &pixels, 
        std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
        DT &com_x, DT &com_y, int n);

template<typename DT> 
void _swap(DT &a, DT &b);

template <typename DT>
void _bubble_sort(std::unique_ptr <DT[]> &vals, 
        std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
        int n);

#endif
