#ifndef _PHOTONS_H_
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <memory>

#define CENTROIDS_CENT_PIXEL                  0x8000
#define CENTROIDS_TABLE_COLS                  8

#ifdef DEBUG_OUTPUT
#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILENAME__, __LINE__, __func__, __VA_ARGS__);
#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, __FILENAME__, __LINE__, __func__); 
#else
#define DEBUG_PRINT(fmt, ...) \
    do {} while (0)
#define DEBUG_COMMENT(fmt) \
    do {} while (0)
#endif

/* ----------------------------------------------------------------------------*/
/**
 * \brief Data structre for pixel LUT 
 *
 * @tparam OT
 */
/* ----------------------------------------------------------------------------*/
template <typename OT>
struct centroids_pixel_lut
{
    std::unique_ptr<OT[]> data;
    OT start;
    OT step;
    size_t n_points; 
};

/* ----------------------------------------------------------------------------*/
/**
 * \brief Parameters for controlling the photon counting algorythm
 *
 * \tparam DT Data image datatype
 */
/* ----------------------------------------------------------------------------*/
template <typename DT, typename OT> 
struct centroid_params 
{
    int box;
    int box_t;
    int box_n;
    int pixel_photon_num;
    int overlap_max;
    double sum_min;
    double sum_max;
    DT threshold; 
    size_t x;
    size_t y;
    size_t n;
    int store_pixels;
    int fit_pixels;
};

/* ----------------------------------------------------------------------------*/
/**
 * \brief 
 *
 * @tparam DT
 */
/* ----------------------------------------------------------------------------*/
template <typename DT> struct photons {
    DT *image;
    uint16_t *out;
    size_t x;
    size_t y;
};

template<typename T>
using PixelValues = std::vector<T>;
template<typename T>
using PixelValuesPtr = std::shared_ptr<PixelValues<T>>;

//template <typename OT> struct photon_params {
//    OT x;
//    OT y;
//    OT com_x;
//    OT com_y;
//    OT sum;
//    OT bgnd;
//    int box_sum;
//    PixelValuesPtr<OT> values;
//};

enum {
    CENTROIDS_PARAMS_OK = 0,
    CENTROIDS_PARAMS_BAD = 1
};

enum {
    CENTROIDS_LUT_OK = 0,
    CENTROIDS_LUT_RANGE_LOW = 1,
    CENTROIDS_LUT_RANGE_HIGH = 2
};

template <typename DT> 
using PhotonTable = std::vector<DT>;
template <typename DT> 
using PhotonTablePtr = PhotonTable<DT>*;
template <typename DT> 
using PhotonMap = std::vector<photons<DT>>;
template <typename DT> 
using PhotonMapPtr = std::unique_ptr<PhotonMap<DT>>;

template <typename DT, typename OT> 
void centroids_initialize_params(centroid_params<DT, OT> *params);

template <typename DT, typename OT> 
int centroids_calculate_params(centroid_params<DT, OT> *params);

template <typename OT>
int centroids_init_pixel_lut(centroids_pixel_lut<OT> &lut,
                             OT start, OT stop, int points);

template <typename OT>
int centroids_calculate_pixel_lut(centroids_pixel_lut<OT> &lut,
                                  OT start, OT stop, int points);

template <typename OT>
int centroids_lookup_pixel_lut(centroids_pixel_lut<OT> &lut, 
                               OT ival, OT &oval);

template<typename DT, typename OT>
size_t centroids_process(DT *image, uint16_t *out, 
                         PhotonTablePtr<OT> &photon_table,
                         size_t X, size_t Y, size_t N, 
                         centroid_params<DT, OT> &params);

template<typename DT, typename OT>
size_t centroids_process_photons(PhotonMapPtr<DT> &photon_map,
                                 PhotonTablePtr<OT> &photon_table,
                                 centroids_pixel_lut<OT> &pixel_lut,
                                 centroid_params<DT, OT> &params);

template<typename DT, typename OT> 
size_t centroids_find_photons(DT *image, uint16_t *out, 
                              PhotonMapPtr<DT> &photon_map,
                              size_t X, size_t Y, 
                              centroid_params<DT, OT> &params);

template <typename DT>
int _calculate_com(std::unique_ptr <DT[]> &pixels, 
                   std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
                   DT &com_x, DT &com_y, int n);

template<typename DT> 
void centroids_swap(DT &a, DT &b);

template <typename DT>
void _bubble_sort(std::unique_ptr <DT[]> &vals, 
        std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
        int n);

#endif
