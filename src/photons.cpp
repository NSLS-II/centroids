#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "photons.h"

using namespace std;

template <typename DataType> void centroids_initialize_params(centroid_params<DataType> &params)
{
    params.box = 2;
    params.box_n = 5;
    params.box_t = 25;

    params.pixel_photon_num = 9;
    params.overlap_max = 0;
    params.sum_min = 0;
    params.sum_max = 10000;
    params.threshold = 100;
}

template <typename DataType> int centroids_calculate_params(centroid_params<DataType> &params)
{
    params.box_n = (params.box * 2) + 1;
    params.box_t = params.box_n * params.box_n;

    return CENTROIDS_PARAMS_OK;
}

/* ----------------------------------------------------------------------------*/
/**
 * \brief Swap the values of two variables
 *
 * @tparam DT
 * \param a swap variable
 * \param b swap variable
 */
/* ----------------------------------------------------------------------------*/
template<typename DT> 
void _swap(DT &a, DT &b)
{
    DT t = a;
    a = b;
    b = t;
}

/* ----------------------------------------------------------------------------*/
/**
 * \brief Bubble sort values in decending order
 *
 * \param vals Pointer to array to sort
 * \param x Pointer to x values to return sorted by vals
 * \param y Pointer to y values to return sorted by vals
 * \param n Number of elements in array
 */
/* ----------------------------------------------------------------------------*/
template <typename DT>
void _bubble_sort(std::unique_ptr <DT[]> &vals, 
        std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
        int n)
{
    for (int i=0;i<(n-1);i++) 
    {
        for (int j=1;j<(n - i);j++) 
        {
            if (vals[j-1] < vals[j]) 
            {
                _swap<DT>(vals[j-1], vals[j]);
                _swap<DT>(x[j-1], x[j]);
                _swap<DT>(y[j-1], y[j]);
            }
        }
    }
}

/* ----------------------------------------------------------------------------*/
/**
 * \brief Process a single CCD image to find photons
 *
 * @tparam DataType
 * \param image Pointer to image data
 * \param out Pointer to store image 
 * \param table 
 * \param bias
 * \param X
 * \param Y
 * \param params
 *
 * \returns   
 */
/* ----------------------------------------------------------------------------*/
template<typename DataType, typename OutputType>
int centroids_process(DataType *image, uint16_t *out, 
        PhotonTablePtr<OutputType> &photon_table,
        size_t X, size_t Y, size_t N, centroid_params<DataType> params)
{
    int nphotons = 0;
    DataType *image_p = image;
    uint16_t *out_p = out;

    // Make the arrays for each image
    PhotonMapPtr<DataType> photon_map(new PhotonMap<DataType>);

    // Allocate a single image

    DEBUG_COMMENT("Looping through images....\n");
    for(size_t n=0;n<N;n++)
    {
        size_t n_photons;

        n_photons = centroids_find_photons<DataType>(image_p, out_p, photon_map, X, Y, params); 
        
        DEBUG_PRINT("find_photons returns %ld (%ld)\n", (unsigned long)n_photons, (unsigned long)n);
        
        n_photons = centroids_process_photons<DataType, OutputType>(photon_map, photon_table, params);

        DEBUG_PRINT("process_photons returns %ld (%ld)\n", (unsigned long)n_photons, (unsigned long)n);

        image_p += (X * Y);
        out_p += (X * Y);
        nphotons += n_photons;
    }

    return nphotons;
}

template <typename DT>
int _calculate_com(std::unique_ptr <DT[]> &pixels, 
        std::unique_ptr <DT[]> &x, std::unique_ptr <DT[]> &y, 
        DT &com_x, DT &com_y, int n)
{
    DT sum = 0;

    com_x = 0;
    com_y = 0;

    for(int i=0;i<n;i++)
    {
        sum += (DT)pixels[i];
        com_x += (pixels[i] * (DT)x[i]);
        com_y += (pixels[i] * (DT)y[i]);
    }

    com_y /= sum;
    com_x /= sum;
    
    return 0;
}

template<typename DataType, typename OutputType>
size_t centroids_process_photons(
        PhotonMapPtr<DataType> &photon_map,
        PhotonTablePtr<OutputType> &photon_table,
        centroid_params<DataType> params)
{
    //int *xvals = new int [params.box_t];
    //int *yvals = new int [params.box_t];
    //OutputType *pixel_cluster = new double [params.box_t];
    std::unique_ptr<OutputType[]> pixel_cluster(new double[params.box_t]);
    std::unique_ptr<OutputType[]> xvals(new double[params.box_t]);
    std::unique_ptr<OutputType[]> yvals(new double[params.box_t]);

    for(auto photon = photon_map->begin(); 
        photon != photon_map->end(); 
        photon+=params.box_t)
    {
        int box_sum = -params.box_t;
        for(int m=0; m<params.box_t; m++)
        {
            box_sum += (int)(*(photon[m].out) & 0x7FFF);
            pixel_cluster[m] = (double)(*photon[m].image);
            xvals[m] = (int)photon[m].x; 
            yvals[m] = (int)photon[m].y; 
            DEBUG_PRINT("pixel(%d) = %lf x = %d y = %d\n", m,
                    (double)pixel_cluster[m], (int)xvals[m], (int)yvals[m]);
        }

        //size_t i = xvals[params.box_t / 2];
        //size_t j = yvals[params.box_t / 2];

        if(box_sum <= params.overlap_max)
        {
            double comx, comy;

            // Bubble sort values
            _bubble_sort<OutputType>(pixel_cluster, xvals, yvals, params.box_t);

            // Now we process background
            double bgnd = 0;
            for(int n=params.pixel_photon_num;n<params.box_t;n++)
            {
                bgnd += pixel_cluster[n];
            }
            bgnd /= (params.box_t - params.pixel_photon_num);
            DEBUG_PRINT("bgnd = %lf (%d)\n", (double)bgnd, params.pixel_photon_num);

            // Subtract Background
            for(int n=0;n<params.box_t;n++)
            {
                pixel_cluster[n] -= bgnd;
            }

            // Now process sum
            double sum = 0;
            for(int n=0;n<params.pixel_photon_num;n++)
            {
                sum += pixel_cluster[n];
            }
            DEBUG_PRINT("sum = %lf\n", (double)sum);

            for(int m=0; m<params.box_t; m++)
            {
                DEBUG_PRINT("sorted pixel = %lf x = %lf y = %lf\n", (double)pixel_cluster[m], 
                        (double)xvals[m], (double)yvals[m]);
            }

            // Now check sum
            if((sum >= params.sum_min)  &&
               (sum < params.sum_max))
            {

                // Calculate the COM
                _calculate_com<OutputType>(pixel_cluster, xvals, yvals,
                        comx, comy, params.box_t);
                DEBUG_PRINT("COM = %lf, %lf (%lf, %lf)\n", 
                        (double)comx, (double)comy, 
                        (double)xvals[0], (double)yvals[0]);

                photon_table->push_back( {xvals[0], yvals[0], 
                        comx, comy,
                        sum, bgnd, 
                        box_sum}
                        );

                //for(int n=0;n<params.box_t;n++)
                //{
                //    table_p[8 + n] = pixel_cluster[n];
                //}

                //table_p += (params.box_t + CENTROIDS_TABLE_COLS);

            } // sum is correct
        } // overlap is correct
    }

    DEBUG_PRINT("Vector Size = %ld\n", (unsigned long)photon_table->size());
    
    return photon_table->size();
}

template<typename DataType> 
size_t centroids_find_photons(DataType *image, uint16_t *out, 
        PhotonMapPtr<DataType> &photon_map,
        size_t X, size_t Y, centroid_params<DataType> params)
{
    uint16_t *out_p = out;
    DataType *in_p = image;

    // Blank the output
    for(size_t i=0;i<(X* Y);i++)
    {
        out[i] = 0;
    }

    // Start loop through image, we need to skip by params.box
    for(size_t j=0;j<(size_t)params.box;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for(size_t j=params.box;j<(Y-(size_t)params.box);j++){
        for(size_t i=0;i<(size_t)params.box;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }

        for(size_t i=params.box;i<(X-params.box);i++)
        {
            // This is the main routine. 
            // Is this pixel above threshold
            if(*in_p >= params.threshold)
            {
                // Check if this is the highest pixel
                // Rewind by 1 params.box in X and 1 params.box in Y
                DataType *_in_p = in_p;
                _in_p -= (params.box + (X * params.box));

                int flag = 1;
                for(size_t l=0;l<(size_t)params.box_n;l++)
                {
                    for(size_t k=0;k<(size_t)params.box_n;k++)
                    {
                        if(*_in_p > *in_p)
                        {
                            // We are not the highest pixel, set the flag
                            flag = 0;
                            break;
                        }

                        _in_p++;
                    }

                    if(!flag)
                    {
                        break;
                    }

                    _in_p += (X - params.box_n);
                }

                if(flag)
                {
                    // We are the highest pixel
                    // Mark the image for duplicates

                    // Mark the center pixel using the bitwise value CENTROIDS_CENT_PIXEL
                    // Increment 1 for each surrounding area
                    
                    DEBUG_PRINT("Pixel above threshold (%ld, %ld)\n", i, j);
                    *out_p |= CENTROIDS_CENT_PIXEL; 


                    // Now set the pointers for dealing with the box
                    uint16_t *_out_p = out_p;
                    _out_p = out_p - (params.box + (X * params.box));
                    _in_p = in_p - (params.box + (X * params.box));

                    // Loop over the box
                    for(int l=-params.box;l<=params.box;l++)
                    {
                        for(int k=-params.box;k<=params.box;k++)
                        {
                            // Store the photon map values
                            // {image, out, x, y}
                            DEBUG_PRINT("Storing pixel for (%d, %d) %p %p\n", 
                                    k, l, (void*)_in_p, (void*)_out_p);
                            photon_map->push_back({_in_p, _out_p, i + k, j + l});

                            // Increment the out value
                            (*_out_p)++;

                            // Increment pointers
                            _in_p++;
                            _out_p++;
                        }
                        _out_p += (X - params.box_n);
                        _in_p += (X - params.box_n);
                    }
                } // if(flag)
            } // if(threshold)
            out_p++; in_p++;
        }

        for(size_t i=(X-params.box);i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for(size_t j=Y-params.box;j<Y;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }
    return photon_map->size();
}

template<typename DataType> int centroids_process_bias(DataType *pixels, uint16_t *out, double *bias, 
        size_t X, size_t Y)
{
    DataType *pix_p = pixels;
    uint16_t *out_p = out;
    double *bias_p = bias;
    
    for(size_t j=0;j<Y;j++)
    {
        double _bias = 0;
        size_t n = 0;

        for(size_t i=0;i<X;i++)
        {
            if(*out_p == 0)
            {
                _bias += *pix_p;
                n++;
            }
            out_p++;
            pix_p++;
        }
        if(n)
        {
            _bias /= n;
        } else {
            _bias = 0;
        }

        for(size_t i=0;i<X;i++)
        {
            *(bias_p++) = _bias;
        }
    }

    return 0;
}


// Templates for common datatypes
template void centroids_initialize_params<uint16_t>(centroid_params<uint16_t> &params);
template int centroids_calculate_params<uint16_t>(centroid_params<uint16_t> &params);

template int centroids_process<uint16_t, double>(uint16_t *image, uint16_t *out, 
        PhotonTablePtr<double> &photon_table,
        size_t X, size_t Y, size_t N, centroid_params<uint16_t> params);
