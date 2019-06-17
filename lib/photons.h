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
#ifndef LIB_PHOTONS_H_
#define LIB_PHOTONS_H_

#include <stdio.h>
#include <stdint.h>
#include <lmmin.h>
#include <vector>
#include <memory>

#include "centroids.h"

#ifndef __FILENAME__
#define __FILENAME__ __FILE__
#endif

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

#define UNUSED(expr) do { (void)(expr); } while (0)

const char *centroids_photon_table_names[] = {
    "Pixel X", "Pixel Y",
    "COM X", "COM Y", "COR COM X", "COR COM Y",
    "Int", "Bgnd", "Overlap",
    "Fit X", "Fit Y", "Fit Bgnd", "Fit Amp", "Fit Sigma",
    "Fit Err X", "Fit Err Y", "Fit Err Bgnd", "Fit Err Amp", "Fit Err Sigma",
    "Fit Fnorm", "Fit Outcome", "Fit StdErr"
};


template <typename OT>
struct fit_data_struct {
    OT *x;
    OT *y;
    OT *z;
    OT (*f)(OT x, OT y, const double *p );
};

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

template <typename DT>
using PhotonMap = std::vector<photons<DT>>;

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
size_t centroids_process_photons(PhotonMap<DT> *photon_map,
        PhotonTable<OT> *photon_table, const centroids_pixel_lut<OT> &pixel_lut,
        std::vector<DT> *photons, const centroid_params<DT, OT> &params);

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

template <typename OT>
OT centroids_std_error_estimate(OT *pixels, OT *xvals, OT *yvals,
        double *fit_params, const int N);

template <typename OT>
OT centroids_2dgauss_int(OT x, OT y, const double *p);

template <typename OT>
void centroids_evaluate_2dgauss(const double *par, int m_dat,
        const void *data, double *fvec, int *info);


#endif  // LIB_PHOTONS_H_
