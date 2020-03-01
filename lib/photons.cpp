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

#include <lmmin.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <memory>

#include "photons.h"


/* -------------------------------------------------------------------------*/
/**
 * \brief
 *
 * @tparam OT
 * param lut
 * param start
 * param stop
 * param points
 *
 * Returns
 */
/* -------------------------------------------------------------------------*/
template <typename OT>
int centroids_init_pixel_lut(centroids_pixel_lut<OT> *lut,
                             OT start, OT stop, size_t points) {
    lut->start = start;
    lut->n_points = points + 1;
    lut->step = (stop - start) / points;

    lut->data = std::unique_ptr<OT[]>(new OT[lut->n_points]);

    DEBUG_PRINT("Created lookup table start = %lf step = %lf steps = %lf\n",
            (double)lut->start, (double)lut->step, (double)lut->step);

    return CENTROIDS_LUT_OK;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief
 *
 * @tparam OT
 * param lut
 * param start
 * param stop
 * param points
 *
 * Returns
 */
/* -------------------------------------------------------------------------*/
template <typename OT>
int centroids_calculate_pixel_lut(centroids_pixel_lut<OT> *lut,
                                  OT start, OT stop, size_t points) {
    int rtn = centroids_init_pixel_lut(lut, start, stop, points);
    if (rtn) {
        return rtn;
    }

    for (size_t n = 0; n < lut->n_points; n++) {
        OT x = start + n * lut->step;
        lut->data[n] = x;
    }

    return CENTROIDS_LUT_OK;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief
 *
 * @tparam OT
 * param lut
 * param ival
 * param oval
 *
 * Returns
 */
/* -------------------------------------------------------------------------*/
template <typename OT>
int centroids_lookup_pixel_lut(const centroids_pixel_lut<OT> &lut,
                               const OT ival, OT *oval) {
    OT pos = (ival - lut.start) / lut.step;
    if (pos < 0) {
        DEBUG_PRINT("ival = %lf pos = %lf RANGE LOW\n",
                (double)ival, (double)pos);
        return CENTROIDS_LUT_RANGE_LOW;
    }

    size_t ipos = (size_t)pos;
    if (ipos > lut.n_points) {
        DEBUG_PRINT("ival = %lf ipos = %ld RANGE HIGH\n",
                (double)ival, ipos);
        return CENTROIDS_LUT_RANGE_HIGH;
    }

    *oval = lut.data[ipos];
    DEBUG_PRINT("ival = %lf ipos = %ld oval = %lf\n",
            (double)ival, ipos, (double)(*oval));

    return CENTROIDS_LUT_OK;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief Initialize the paramaters structure
 *
 * @tparam DT Type for input data to centroids
 * param params Parameter structre
 */
/* -------------------------------------------------------------------------*/
template <typename DT, typename OT>
void centroids_initialize_params(centroid_params<DT, OT> *params) {
    params->box = 2;
    params->box_n = 5;
    params->box_t = 25;

    params->search_box = 1;
    params->search_box_n = 3;
    params->search_box_t = 9;

    params->pixel_photon_num = 9;
    params->com_photon_num = 9;
    params->pixel_bgnd_num = 10;
    params->overlap_max = 0;
    params->sum_min = 0;
    params->sum_max = 10000;
    params->threshold = 100;
    params->return_pixels = CENTROIDS_STORE_NONE;
    params->return_map = false;
    params->fit_pixels = 0;

    params->control = lm_control_double;
    params->control.verbosity = 0;
    params->control.patience = 1000;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief Calculate values from initial input to parameters structure
 *
 * @tparam DT Type for input data to centroids
 * param params Parameters structure
 *
 * Returns 0 if compute was correct, non-zero if an error occured
 */
/* -------------------------------------------------------------------------*/
template <typename DT, typename OT>
int centroids_calculate_params(centroid_params<DT, OT> *params) {
    params->box_n = (params->box * 2) + 1;
    params->box_t = params->box_n * params->box_n;
    params->search_box_n = (params->search_box * 2) + 1;
    params->search_box_t = params->search_box_n * params->search_box_n;

    if ((params->pixel_photon_num <= 0)
        || (params->pixel_photon_num >= params->box_t)) {
        return CENTROIDS_PARAMS_BAD;
    }

    if ((params->pixel_bgnd_num <= 0)
        || (params->pixel_bgnd_num >= params->box_t)) {
        return CENTROIDS_PARAMS_BAD;
    }

    if ((params->com_photon_num <= 0)
        || (params->com_photon_num >= params->box_t)) {
        return CENTROIDS_PARAMS_BAD;
    }

    return CENTROIDS_PARAMS_OK;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief Swap the values of two variables
 *
 * @tparam DT
 * \param a swap variable
 * \param b swap variable
 */
/* -------------------------------------------------------------------------*/
template<typename DT>
void centroids_swap(DT *a, DT *b) {
    DT t = *a;
    *a = *b;
    *b = t;
}

/* -------------------------------------------------------------------------*/
/**
 * \brief Bubble sort values in decending order
 *
 * \param vals Pointer to array to sort
 * \param x Pointer to x values to return sorted by vals
 * \param y Pointer to y values to return sorted by vals
 * \param n Number of elements in array
 */
/* -------------------------------------------------------------------------*/
template <typename DT>
void centroids_bubble_sort(DT *vals, DT *x, DT *y, const int n) {
    for (int i = 0; i < (n-1); i++) {
        for (int j = 1; j < (n - i); j++) {
            if (vals[j-1] < vals[j]) {
                centroids_swap<DT>(&vals[j-1], &vals[j]);
                centroids_swap<DT>(&x[j-1], &x[j]);
                centroids_swap<DT>(&y[j-1], &y[j]);
            }
        }
    }
}

template <typename OT>
OT centroids_2dgauss_int(OT x, OT y, const double *p) {
    // p[0] = x0
    // p[1] = y0
    // p[2] = bgnd
    // p[3] = amplitude
    // p[4] = sigma
    double xmin = (x - 0.5 - p[0]) / (p[4] * sqrt(2));
    double xmax = (x + 0.5 - p[0]) / (p[4] * sqrt(2));
    double ymin = (y - 0.5 - p[1]) / (p[4] * sqrt(2));
    double ymax = (y + 0.5 - p[1]) / (p[4] * sqrt(2));
    double out = p[2] + (p[3]
            * (erf(xmin) - erf(xmax)) * (erf(ymin) - erf(ymax)));
    return out;
}

template <typename OT>
OT centroids_1dgauss_int(OT x, const double *p) {
    // p[0] = x0
    // p[1] = bgnd
    // p[2] = amplitude
    // p[3] = sigma

    double xmin = (x - 0.5 - p[0]) / (p[3] * sqrt(2));
    double xmax = (x + 0.5 - p[0]) / (p[3] * sqrt(2));
    double out = p[1] + (p[2]
            * (erf(xmax) - erf(xmin)));
    return out;
}

template <typename OT>
OT centroids_std_error_estimate_1d(OT *pixels,
        double *fit_params, const int N) {
    OT sum = 0;
    OT calc;

    int i = 0;
    for (int x = -N; x <= N; x++) {
        // Evaluate fit function at each point
        calc = centroids_1dgauss_int<OT>((OT)x, fit_params);

        DEBUG_PRINT("FIT actual = %lf fit = %lf residual = %lf\n",
                (double)pixels[i], (double)calc,
                (double)abs(pixels[i] - calc));

        sum += pow(pixels[i] - calc, 2);
        i++;
    }

    sum /= (OT)N;

    return pow(sum, 0.5);
}

template <typename OT>
OT centroids_std_error_estimate_2d(OT *pixels,
        double *fit_params, const int N) {
    OT sum = 0;
    OT calc;

    int i = 0;
    for (int y = -N; y <= N; y++) {
        for (int x = -N; x <= N; x++) {
            // Evaluate fit function at each point
            calc = centroids_2dgauss_int<OT>((OT)x, (OT)y, fit_params);

            DEBUG_PRINT("FIT actual = %lf fit = %lf residual = %lf\n",
                    (double)pixels[i], (double)calc,
                    (double)abs(pixels[i] - calc));

            sum += pow(pixels[i] - calc, 2);
            i++;
        }
    }

    sum /= (OT)N;

    return pow(sum, 0.5);
}

template <typename OT>
void centroids_evaluate_2dgauss(const double *p, int m_dat,
        const void *data, double *fvec, int *info ) {
    UNUSED(info);
    UNUSED(m_dat);
    fit_data_struct<OT> *D = (fit_data_struct<OT>*)data;

    int i = 0;
    for (int y = -D->box; y <= D->box; y++) {
        for (int x = -D->box; x <= D->box; x++) {
            double out = centroids_2dgauss_int((OT)x, (OT)y, p);
            fvec[i] = D->measured[i] - out;
            i++;
        }
    }
}

template <typename OT>
void centroids_evaluate_1dgauss(const double *p, int m_dat,
        const void *data, double *fvec, int *info ) {
    UNUSED(info);
    UNUSED(m_dat);
    fit_data_struct<OT> *D = (fit_data_struct<OT>*)data;

    int i = 0;
    for (int x = -D->box; x <= D->box; x++) {
        double out = centroids_1dgauss_int((OT)x, p);
        fvec[i] = D->measured[i] - out;
        i++;
    }
}

/* -------------------------------------------------------------------------*/
/**
 * \brief Process a single CCD image to find photons
 *
 * @tparam DT
 * \param image Pointer to image data
 * \param out Pointer to store image
 * \param table
 * \param X
 * \param Y
 * \param params
 *
 * \returns
 */
/* -------------------------------------------------------------------------*/
template<typename DT, typename OT>
size_t centroids_process(DT *image, uint16_t *out,
                         PhotonTable<OT> *photon_table,
                         std::vector<DT> *photons,
                         const centroid_params<DT, OT> &params) {
    size_t n_photons = 0;

    // Make the arrays for each image

    // Allocate a single image
    OT start = -1.1;
    OT stop = 1.1;
    size_t np = 1200;
    centroids_pixel_lut<OT> pixel_lut;
    if (centroids_calculate_pixel_lut<OT>(&pixel_lut, start, stop, np)
       != CENTROIDS_LUT_OK) {
        DEBUG_COMMENT("Failed to calculate LUT\n");
        return 0;
    }

#ifdef _OPENMP
    // Output openmp info
    DEBUG_PRINT("omp_get_num_procs() = %d\n", omp_get_num_procs());
    DEBUG_PRINT("omp_get_max_threads() = %d\n", omp_get_max_threads());
    DEBUG_PRINT("omp_get_thread_limit() = %d\n", omp_get_thread_limit());
#else
    DEBUG_COMMENT("No openmp support\n");
#endif

    #pragma omp parallel shared(n_photons)
    {
        PhotonTable<OT> local_photon_table;
        std::vector<DT> local_photons;
        size_t local_n_photons = 0;
        PhotonMap<DT> photon_map;
        uint16_t *out_p;

        if (!params.return_map) {
            out_p = new uint16_t[params.x * params.y];
        } else {
            out_p = out;
        }

        #pragma omp for
        for (size_t n = 0; n < params.n; n++) {
            // Get he pointers for the input and output images
            DT *image_p = image + (n * params.x * params.y);
            size_t fphotons, pphotons;

            if (params.return_map) {
               out_p = out + (n * params.x * params.y);
            }

            photon_map.clear();
            fphotons = centroids_find_photons<DT, OT>(
                    image_p, out_p, &photon_map, params);

            UNUSED(fphotons);

            DEBUG_PRINT("fphotons = %zu, photon_map.size() = %zu\n",
                    fphotons, photon_map.size());

            pphotons = centroids_process_photons<DT, OT>(
                    &photon_map, &local_photon_table, pixel_lut,
                    &local_photons, params);

            DEBUG_PRINT("Image %zu Found %zu photons, processed %zu photons\n",
                    n, fphotons, pphotons);

            local_n_photons += pphotons;
        }

        #pragma omp critical
        {
            photon_table->insert(photon_table->end(),
                    std::make_move_iterator(local_photon_table.begin()),
                    std::make_move_iterator(local_photon_table.end()));

            if (params.return_pixels != CENTROIDS_STORE_NONE) {
                photons->insert(photons->end(),
                        std::make_move_iterator(local_photons.begin()),
                        std::make_move_iterator(local_photons.end()));
            }

            DEBUG_PRINT("local_n_photons = %zu\n", local_n_photons);
            n_photons += local_n_photons;

            if (!params.return_map) {
                delete [] out_p;
            }
        }
    }

    DEBUG_PRINT("Found %zu photons\n", n_photons);
    DEBUG_PRINT("photon_map->size() = %zu\n", photon_table->size());

    return n_photons;
}

template <typename DT>
int centroids_calculate_com(const std::unique_ptr <DT[]> &pixels,
                            const std::unique_ptr <DT[]> &x,
                            const std::unique_ptr <DT[]> &y,
                            DT *com_x, DT *com_y, const int n) {
    DT sum = 0;

    *com_x = 0;
    *com_y = 0;

    for (int i = 0; i < n; i++) {
        sum += (DT)pixels[i];
        *com_x += (pixels[i] * (DT)x[i]);
        *com_y += (pixels[i] * (DT)y[i]);
    }

    *com_y /= sum;
    *com_x /= sum;

    return 0;
}

template<typename DT, typename OT>
void centroids_fit_photon(double *fit_params,
        fit_data_struct<OT> *fit_data, PhotonTable<OT> *photon_table,
        const bool fit_2d,
        const centroid_params<DT, OT> &params, int n_params) {
    lm_status_struct fit_status;
    double fit_params_err[CENTROIDS_FIT_PARAMS_MAX];

    if (fit_2d) {
        lmmin2(n_params, fit_params, fit_params_err,
                NULL, params.box_t, NULL, (const void*) fit_data,
                centroids_evaluate_2dgauss<OT>, &params.control,
                &fit_status);
    } else {
        lmmin2(n_params, fit_params, fit_params_err,
                NULL, params.box_t, NULL, (const void*) fit_data,
                centroids_evaluate_1dgauss<OT>, &params.control,
                &fit_status);
    }

#ifdef DEBUG_OUTPUT
    DEBUG_PRINT("FIT : status after %d evaluations:  %s\n",
            fit_status.nfev, lm_shortmsg[fit_status.outcome]);

    for (int i = 0; i < n_params; i++) {
        DEBUG_PRINT("FIT : par[%i] = %12g +- %12g\n",
                i, fit_params[i], fit_params_err[i]);
    }
#endif

    OT std_err;
    if (fit_2d) {
        std_err = centroids_std_error_estimate_2d<OT>(
             fit_data->measured, fit_params, params.box);
    } else {
        std_err = centroids_std_error_estimate_1d<OT>(
             fit_data->measured, fit_params, params.box);
    }

    DEBUG_PRINT("FIT std err = %lf\n", (double)std_err);

    // Return fit in absolute pixel coords
    fit_params[0] += fit_data->x;
    if (fit_2d) {
        fit_params[1] += fit_data->y;
    }

    photon_table->insert(photon_table->end(), &fit_params[0],
            &fit_params[n_params]);
    photon_table->insert(photon_table->end(), &fit_params_err[0],
            &fit_params_err[n_params]);

    photon_table->insert(photon_table->end(),
            {(OT)fit_status.fnorm, (OT)fit_status.outcome,
            std_err});
}

template<typename DT, typename OT>
size_t centroids_process_photons(PhotonMap<DT> *photon_map,
        PhotonTable<OT> *photon_table, const centroids_pixel_lut<OT> &pixel_lut,
        std::vector<DT> *photons, const centroid_params<DT, OT> &params) {
    size_t n_photons = 0;

    std::unique_ptr<OT[]> pixel_cluster(new OT[params.box_t]);
    std::unique_ptr<OT[]> pixel_cluster_fit(new OT[params.box_t]);
    std::unique_ptr<OT[]> pixel_cluster_fit_int(new OT[params.box_n]);
    std::unique_ptr<OT[]> xvals(new OT[params.box_t]);
    std::unique_ptr<OT[]> yvals(new OT[params.box_t]);

    double fit_params[CENTROIDS_FIT_PARAMS_MAX];

    // Setup and store structures for fitting of
    // pixel data using liblmfit

    fit_data_struct<OT> fit_data_2d =
        { pixel_cluster_fit.get(), 0, 0, params.box};

    fit_data_struct<OT> fit_data_1d =
        { pixel_cluster_fit_int.get(), 0, 0, params.box};

    // Loop over all pixel clusters
    for (auto photon = photon_map->begin();
         photon != photon_map->end(); photon+=params.box_t) {
        int box_sum = -params.box_t;
        OT comx  = 0;
        OT comy  = 0;
        OT ccomx = 0;
        OT ccomy = 0;
        OT sum   = 0;
        OT bgnd  = 0;

        // Store the pixel values in working structures
        for (int m = 0; m < params.box_t; m++) {
            box_sum += static_cast<int>(*(photon[m].out) & 0x7FFF);
            pixel_cluster[m] = static_cast<double>(*photon[m].image);
            pixel_cluster_fit[m] = pixel_cluster[m];
            xvals[m] = static_cast<int>(photon[m].x);
            yvals[m] = static_cast<int>(photon[m].y);
            DEBUG_PRINT("pixel(%d) = %lf x = %lf y = %lf\n", m,
                    (double)pixel_cluster[m],
                    (double)xvals[m], (double)yvals[m]);
        }

        // Check the box sum, if it is too great
        // then skip this pixel
        if (box_sum > params.overlap_max) {
            DEBUG_COMMENT("Skipping due to overlap.");
            continue;
        }

        // -----------------------------
        // Bubble sort pixel intensities
        // -----------------------------

        centroids_bubble_sort<OT>(pixel_cluster.get(),
                xvals.get(), yvals.get(),
                params.box_t);

        // ------------------
        // Process background
        // ------------------

        bgnd = 0;
        for (int n = params.pixel_bgnd_num; n < params.box_t; n++) {
            bgnd += pixel_cluster[n];
        }
        bgnd /= (params.box_t - params.pixel_bgnd_num);

        DEBUG_PRINT("bgnd = %lf (%d)\n",
                (double)bgnd, params.pixel_photon_num);

        for (int n = 0; n < params.box_t; n++) {
            pixel_cluster[n] -= bgnd;
        }

        // --------------------------------
        // Process integral (sum) of pixels
        // --------------------------------

        sum = 0;
        for (int n = 0; n < params.pixel_photon_num; n++) {
            sum += pixel_cluster[n];
        }
        DEBUG_PRINT("sum = %lf\n", (double)sum);

#ifdef DEBUG_OUTPUT
        for (int m = 0; m < params.box_t; m++) {
            DEBUG_PRINT("sorted pixel = %lf x = %lf y = %lf\n",
                    (double)pixel_cluster[m], (double)xvals[m],
                    (double)yvals[m]);
        }
#endif

        // ------------------------------------
        // Now check sum (total charge in spot)
        // ------------------------------------

        if ((sum < params.sum_min) || (sum > params.sum_max)) {
            DEBUG_COMMENT("Skipping due to total charge in box.");
            continue;
        }

        // -----------------------------------------------------------
        // Calculate the COM
        // TODO(stuwilkins) : We can place a limit here to improve COM
        // -----------------------------------------------------------

        centroids_calculate_com<OT>(pixel_cluster, xvals, yvals,
                &comx, &comy, params.com_photon_num);
        DEBUG_PRINT("COM = %lf, %lf (%lf, %lf)\n",
                (double)comx, (double)comy,
                (double)xvals[0], (double)yvals[0]);

        // Correct COM using LUT for pixels
        centroids_lookup_pixel_lut<OT>(pixel_lut, comx - xvals[0], &ccomx);
        centroids_lookup_pixel_lut<OT>(pixel_lut, comy - yvals[0], &ccomy);

        // Insert result into photon table
        photon_table->insert(photon_table->end(),
                {xvals[0], yvals[0],
                 comx, comy,
                 xvals[0] + ccomx, yvals[0] + ccomy,
                 sum, bgnd, (OT)box_sum});

        // --------------------------------
        // Now fit results if required (2D)
        // --------------------------------

        if (params.fit_pixels & CENTROIDS_FIT_2D) {
            // Set the initial guess values
            fit_params[0] = 0;
            fit_params[1] = 0;
            fit_params[2] = bgnd;
            fit_params[3] = pixel_cluster[0];
            fit_params[4] = 0.5;

            fit_data_2d.x = xvals[0];
            fit_data_2d.y = yvals[0];

            centroids_fit_photon<DT, OT>(fit_params, &fit_data_2d,
                    photon_table, true, params, CENTROIDS_FIT_PARAMS_2D_N);
        }

        // ---------------------------------
        // Now fit results if required (1DX)
        // ---------------------------------

        if (params.fit_pixels & CENTROIDS_FIT_1D_X) {
            // Integrate the values

            for (int j = 0; j < params.box_n; j++) {
                pixel_cluster_fit_int[j] = 0;
                for (int i = 0; i < params.box_n; i++) {
                    pixel_cluster_fit_int[j] +=
                        pixel_cluster_fit[i * params.box_n + j];
                }
            }

            // Set the initial guess values
            fit_params[0] = 0;
            fit_params[1] = bgnd * params.box_n;
            fit_params[2] = 500;
            fit_params[3] = 0.5;
            fit_data_1d.x = xvals[0];

            centroids_fit_photon<DT, OT>(fit_params, &fit_data_1d,
                    photon_table, false, params, CENTROIDS_FIT_PARAMS_1D_N);
        }

        // ---------------------------------
        // Now fit results if required (1DY)
        // ---------------------------------

        if (params.fit_pixels & CENTROIDS_FIT_1D_Y) {
            // Integrate the values

            for (int j = 0; j < params.box_n; j++) {
                pixel_cluster_fit_int[j] = 0;
                for (int i = 0; i < params.box_n; i++) {
                    pixel_cluster_fit_int[j] +=
                        pixel_cluster_fit[j * params.box_n + i];
                }
            }

            // Set the initial guess values
            fit_params[0] = 0;
            fit_params[1] = bgnd * params.box_n;
            fit_params[2] = 500;
            fit_params[3] = 0.5;
            fit_data_1d.x = yvals[0];

            centroids_fit_photon<DT, OT>(fit_params, &fit_data_1d,
                    photon_table, false, params, CENTROIDS_FIT_PARAMS_1D_N);
        }

        // --------------------------
        // Store pixel cluster values
        // --------------------------

        if (params.return_pixels == CENTROIDS_STORE_SORTED) {
            photons->insert(photons->end(),
                    &pixel_cluster[0], &pixel_cluster[params.box_t]);
        } else if (params.return_pixels == CENTROIDS_STORE_UNSORTED) {
            for (int m = 0; m < params.box_t; m++) {
                photons->push_back(*(photon[m].image));
            }
        }

        n_photons++;
    }

    DEBUG_PRINT("Vector Size = %zu\n", photon_table->size());
    DEBUG_PRINT("n_photons = %zu\n", n_photons);

    return n_photons;
}

template<typename DT, typename OT>
size_t centroids_find_photons(DT *image, uint16_t *out,
                              PhotonMap<DT> *photon_map,
                              const centroid_params<DT, OT> &params) {
    uint16_t *out_p = out;
    DT *in_p = image;
    size_t n_photons = 0;

    // Blank the output
    for (size_t i = 0; i < (params.x * params.y); i++) {
        out[i] = 0;
    }

    // Start loop through image, we need to skip by params.box
    for (size_t j = 0; j < (size_t)params.box; j++) {
        for (size_t i = 0; i < params.x; i++) {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for (size_t j = params.box; j < (params.y - (size_t)params.box); j++) {
        for (size_t i = 0; i < (size_t)params.box; i++) {
            *(out_p++) = 1;
            in_p++;
        }

        for (size_t i = params.box; i < (params.x - params.box); i++) {
            // This is the main routine.
            // Is this pixel above threshold
            if (*in_p >= params.threshold) {
                // Check if this is the highest pixel
                // Rewind by 1 params.box in X and 1 params.box in y

                DT *_in_p = in_p;
                _in_p -= (params.search_box + (params.x * params.search_box));

                int flag = 1;
                for (size_t l = 0; l < (size_t)params.search_box_n; l++) {
                    for (size_t k = 0; k < (size_t)params.search_box_n; k++) {
                        if (*_in_p > *in_p) {
                            // We are not the highest pixel, set the flag
                            flag = 0;
                            break;
                        }

                        _in_p++;
                    }

                    if (!flag) {
                        break;
                    }

                    _in_p += (params.x - params.search_box_n);
                }

                if (flag) {
                    // We are the highest pixel
                    // Mark the image for duplicates

                    // Mark the center pixel using the bitwise value
                    // Increment 1 for each surrounding area

                    DEBUG_PRINT("Pixel above threshold (%ld, %ld)\n", i, j);

                    // Now set the pointers for dealing with the box
                    uint16_t *_out_p = out_p;
                    _out_p = out_p - (params.box + (params.x * params.box));
                    _in_p = in_p - (params.box + (params.x * params.box));

                    // Loop over the box
                    for (int l = -params.box; l <= params.box; l++) {
                        for (int k = -params.box; k <= params.box; k++) {
                            // Store the photon map values
                            // {image, out, x, y}
                            DEBUG_PRINT("Storing pixel for (%d, %d) %p %p\n",
                                    k, l, (void*)_in_p, (void*)_out_p);
                            photon_map->push_back(
                                    {_in_p, _out_p, i + k, j + l});

                            // Increment the out value
                            (*_out_p)++;

                            // Increment pointers
                            _in_p++;
                            _out_p++;
                        }
                        _out_p += (params.x - params.box_n);
                        _in_p += (params.x - params.box_n);
                    }
                    n_photons++;
                }  // if(flag)
            }  // if(threshold)
            out_p++; in_p++;
        }

        for (size_t i = (params.x - params.box); i < params.x; i++) {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for (size_t j = (params.y - params.box); j < params.y ; j++) {
        for (size_t i = 0; i < params.x; i++) {
            *(out_p++) = 1;
            in_p++;
        }
    }
    return n_photons;
}

// Templates for common datatypes
template void centroids_initialize_params<uint16_t, double>(
        centroid_params<uint16_t, double> *params);
template int centroids_calculate_params<uint16_t, double>(
        centroid_params<uint16_t, double> *params);
template size_t centroids_process<uint16_t, double>(
        uint16_t *image, uint16_t *out, PhotonTable<double> *photon_table,
        std::vector<uint16_t> *photons,
        const centroid_params<uint16_t, double> &params);

template void centroids_initialize_params<uint16_t, float>(
        centroid_params<uint16_t, float> *params);
template int centroids_calculate_params<uint16_t, float>(
        centroid_params<uint16_t, float> *params);
template size_t centroids_process<uint16_t, float>(
        uint16_t *image, uint16_t *out, PhotonTable<float> *photon_table,
        std::vector<uint16_t> *photons,
        const centroid_params<uint16_t, float> &params);
