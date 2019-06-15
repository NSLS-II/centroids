//
// CENTROIDS : C++ implementation of single photon counting for CCDs
// Stuart B. Wilkins, Brookhaven National Laboratory
//
//
// BSD 3-Clause License
//
// Copyright (c) 2019, Brookhaven Science Associates
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
// THE POSSIBILITY OF SUCH DAMAGE.
//
#ifndef SRC_PHOTONS_H_
#define SRC_PHOTONS_H_

#include <stdio.h>
#include <stdint.h>
#include <lmmin.h>
#include <vector>
#include <memory>

#define CENTROIDS_CENT_PIXEL                  0x8000
#define CENTROIDS_TABLE_COLS                  8
#define CENTROIDS_FIT_PARAMS_N                5

#ifdef DEBUG_OUTPUT
#define DEBUG_PRINT(fmt, ...) \
  fprintf(stderr, "%s:%d:%s(): " fmt, \
          __FILENAME__, __LINE__, __func__, __VA_ARGS__);
#define DEBUG_COMMENT(fmt) \
  fprintf(stderr, "%s:%d:%s(): " fmt, \
          __FILENAME__, __LINE__, __func__);
#else
#define DEBUG_PRINT(fmt, ...) \
    do {} while (0)
#define DEBUG_COMMENT(fmt) \
    do {} while (0)
#endif

typedef struct {
    double *tx, *ty;
    double *y;
    double (*f)( double tx, double tz, const double *p );
} fit_data_struct;

/* -------------------------------------------------------------------------*/
/**
 * \brief Data structre for pixel LUT
 *
 * @tparam OT
 */
/* -------------------------------------------------------------------------*/
template <typename OT>
struct centroids_pixel_lut {
    std::unique_ptr<OT[]> data;
    OT start;
    OT step;
    size_t n_points;
};

/* -------------------------------------------------------------------------*/
/**
 * \brief Parameters for controlling the photon counting algorythm
 *
 * \tparam DT Data image datatype
 */
/* -------------------------------------------------------------------------*/
template <typename DT, typename OT>
struct centroid_params {
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
    lm_control_struct control;
};

/* -------------------------------------------------------------------------*/
/**
 * \brief
 *
 * @tparam DT
 */
/* -------------------------------------------------------------------------*/
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

enum {
    CENTROIDS_PARAMS_OK = 0,
    CENTROIDS_PARAMS_BAD = 1
};

enum {
    CENTROIDS_LUT_OK = 0,
    CENTROIDS_LUT_RANGE_LOW = 1,
    CENTROIDS_LUT_RANGE_HIGH = 2
};

enum {
    CENTROIDS_STORE_NONE = 0,
    CENTROIDS_STORE_SORTED = 1,
    CENTROIDS_STORE_UNSORTED = 2
};

enum {
    CENTROIDS_FIT_NONE = 0,
    CENTROIDS_FIT_LMMIN = 1
};

template <typename DT>
using PhotonTable = std::vector<DT>;
template <typename DT>
using PhotonMap = std::vector<photons<DT>>;

template <typename DT, typename OT>
void centroids_initialize_params(centroid_params<DT, OT> *params);

template <typename DT, typename OT>
int centroids_calculate_params(centroid_params<DT, OT> *params);

template <typename OT>
int centroids_init_pixel_lut(centroids_pixel_lut<OT> *lut,
                             OT start, OT stop, int points);

template <typename OT>
int centroids_calculate_pixel_lut(centroids_pixel_lut<OT> *lut,
                                  OT start, OT stop, int points);

template <typename OT>
int centroids_lookup_pixel_lut(const centroids_pixel_lut<OT> &lut,
                               const OT ival, OT *oval);

template<typename DT, typename OT>
size_t centroids_process(DT *image, uint16_t *out,
                         PhotonTable<OT> *photon_table,
                         const centroid_params<DT, OT> &params);

template<typename DT, typename OT>
size_t centroids_process_photons(PhotonMap<DT> *photon_map,
                                 PhotonTable<OT> *photon_table,
                                 const centroids_pixel_lut<OT> &pixel_lut,
                                 centroid_params<DT, OT> *params);

template<typename DT, typename OT>
size_t centroids_find_photons(DT *image, uint16_t *out,
                              PhotonMap<DT> *photon_map,
                              const centroid_params<DT, OT> &params);

template <typename DT>
int centroids_calculate_com(const std::unique_ptr <DT[]> &pixels,
                            const std::unique_ptr <DT[]> &x,
                            const std::unique_ptr <DT[]> &y,
                            DT *com_x, DT *com_y, const int n);

template<typename DT>
void centroids_swap(DT *a, DT *b);

template <typename DT>
void centroids_bubble_sort(DT *vals, DT *x, DT *y, const int n);

#endif  // SRC_PHOTONS_H_
