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
#ifndef LIB_CENTROIDS_H_
#define LIB_CENTROIDS_H_

#include <lmmin.h>
#include <vector>

extern const char* CENTROIDS_GIT_REV;
extern const char* CENTROIDS_GIT_BRANCH;
extern const char* CENTROIDS_GIT_VERSION;

#define CENTROIDS_TABLE_COLS                  9
#define CENTROIDS_FIT_PARAMS_MAX              5
#define CENTROIDS_FIT_PARAMS_2D_N             5
#define CENTROIDS_FIT_PARAMS_1D_N             4
#define CENTROIDS_FIT_EXTRA_N                 3

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
    CENTROIDS_FIT_2D = 1,
    CENTROIDS_FIT_1D_X = 2,
    CENTROIDS_FIT_1D_Y = 4
};

extern const char *centroids_photon_table_names[];
extern const char *centroids_photon_table_names_fit2d[];
extern const char *centroids_photon_table_names_fit1dx[];
extern const char *centroids_photon_table_names_fit1dy[];

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
    int search_box;
    int search_box_t;
    int search_box_n;
    int pixel_photon_num;
    int pixel_bgnd_num;
    int com_photon_num;
    int overlap_max;
    double sum_min;
    double sum_max;
    DT threshold;
    size_t x;
    size_t y;
    size_t n;
    int return_pixels;
    bool return_map;
    int fit_pixels;
    lm_control_struct control;
};

template <typename DT>
using PhotonTable = std::vector<DT>;

template <typename DT, typename OT>
void centroids_initialize_params(centroid_params<DT, OT> *params);

template <typename DT, typename OT>
int centroids_calculate_params(centroid_params<DT, OT> *params);

template<typename DT, typename OT>
size_t centroids_process(DT *image, uint16_t *out,
                         PhotonTable<OT> *photon_table,
                         std::vector<DT> *photons,
                         const centroid_params<DT, OT> &params);
#endif  // LIB_CENTROIDS_H_
