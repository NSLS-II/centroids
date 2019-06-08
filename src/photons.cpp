#include <stdint.h>
#include <stdlib.h>
#include "photons.h"

using namespace std;

template <typename DataType> int centroids_initialize_params(centroid_params<DataType> &params)
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
template<typename DT> void _swap(DT *a, DT *b)
{
    DT t = *a;
    *a = *b;
    *b = t;
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
void _bubble_sort(double *vals, int *x, int *y, int n)
{
    for (int i=0;i<(n-1);i++) 
    {
        for (int j=1;j<(n - i);j++) 
        {
            if (vals[j-1] < vals[j]) 
            {
                _swap<double>(&vals[j-1], &vals[j]);
                _swap<int>(&x[j-1], &x[j]);
                _swap<int>(&y[j-1], &y[j]);
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
template<typename DataType> int centroids_process(DataType *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, size_t N, centroid_params<DataType> params)
{
    int nphotons = 0;
    DataType *image_p = image;
    uint16_t *out_p = out;
    double *table_p = table;
    double *bias_p = bias;

    // Make the arrays for each image
    photons<DataType> map[params.box_t * 10000];
    photons<DataType> *map_p = map;

    for(int n=0;n<N;n++)
    {
        int n_photons;

        n_photons = centroids_find_photons<DataType>(image_p, out_p, map_p, X, Y, params); 
        
        DEBUG_PRINT("find_photons returns %d (%d)\n", n_photons, n);
        
        n_photons = centroids_process_photons<DataType>(map_p, table_p, n_photons, params);

        image_p += (X * Y);
        out_p += (X * Y);
        bias_p += (X * Y);
        table_p += n_photons * (params.box_t + CENTROIDS_TABLE_COLS);
        nphotons += n_photons;
    }

    return nphotons;
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

int _calculate_com(double *pixels, int *x, int *y, double *com_x, double *com_y, int n)
{
    double *pixels_p = pixels;
    int *x_p = x;
    int *y_p = y;
    int sum = 0;

    *com_x = 0;
    *com_y = 0;

    for(int i=0;i<n;i++)
    {
        sum += (*pixels_p); 
        *com_x += (*pixels_p) * *(x_p++);
        *com_y += (*pixels_p) * *(y_p++);
        pixels_p++;
    }

    *com_y /= sum;
    *com_x /= sum;
    
    return 0;
}

template<typename DataType> int centroids_process_photons(photons<DataType> *photon_map,
        double *photon_table, int n_photons, centroid_params<DataType> params)
{
    int xvals[params.box_t];
    int yvals[params.box_t];
    double pixel_cluster[params.box_t];

    photons<DataType> *map_p = photon_map;
    double *table_p = photon_table;

    int table_n = 0;

    for(int n=0;n<n_photons;n++)
    {
        DEBUG_PRINT("n = %d\n", n);
        int box_sum = -params.box_t;
        for(int m=0; m<params.box_t; m++)
        {
            box_sum += (*(map_p[m].out) & 0x7FFF);
            pixel_cluster[m] = (double)*(map_p[m].image);
            DEBUG_PRINT("pixel = %d\n", (int)pixel_cluster[m]);
            xvals[m] = map_p[m].x; 
            yvals[m] = map_p[m].y; 
        }

        int i = xvals[params.box_t / 2];
        int j = yvals[params.box_t / 2];

        DEBUG_PRINT("i = %d, j = %d\n", i, j);
        DEBUG_PRINT("box_sum = %d\n", box_sum);

        if(box_sum <= params.overlap_max)
        {
            double comx, comy;

            _bubble_sort(pixel_cluster, xvals, yvals, params.box_t);

            // Now we process background
            double bgnd = 0;
            for(int n=params.pixel_photon_num;n<params.box_t;n++)
            {
                bgnd += pixel_cluster[n];
            }
            bgnd /= (params.box_t - params.pixel_photon_num);
            DEBUG_PRINT("Pixel bgnd = %lf\n", bgnd);

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

            DEBUG_PRINT("Pixel sum = %lf\n", sum);

            // Now check sum
            if((sum >= params.sum_min)  &&
               (sum < params.sum_max))
            {

                // Calculate the COM
                _calculate_com(pixel_cluster, xvals, yvals,
                        &comx, &comy, params.box_t);

                table_n++;
                table_p[0] = i;  
                table_p[1] = j;  
                table_p[2] = comx + i;
                table_p[3] = comy + j;
                table_p[4] = sum;
                table_p[5] = bgnd;
                table_p[6] = box_sum;
                table_p[7] = 0;

                for(int n=0;n<params.box_t;n++)
                {
                    table_p[8 + n] = pixel_cluster[n];
                }

                table_p += (params.box_t + CENTROIDS_TABLE_COLS);

            } // sum is correct
        } // overlap is correct
        map_p += params.box_t;
    }

    DEBUG_PRINT("table_n = %d\n", table_n);

    return table_n;
}

template<typename DataType> int centroids_find_photons(DataType *image, uint16_t *out, 
        photons<DataType> *photon_map, size_t X, size_t Y, centroid_params<DataType> params)
{
    uint16_t *out_p = out;
    DataType *in_p = image;
    photons<DataType> *map_p = photon_map;

    int n_photons = 0;

    // Blank the output
    for(size_t i=0;i<(X* Y);i++)
    {
        out[i] = 0;
    }

    // Start loop through image, we need to skip by params.box
    for(size_t j=0;j<params.box;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for(size_t j=params.box;j<(Y-params.box);j++){
        for(size_t i=0;i<params.box;i++)
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
                DEBUG_PRINT("Pixel above threshold (%d, %d)\n", (int)i, (int)j);
                // Check if this is the highest pixel
                // Rewind by 1 params.box in X and 1 params.box in Y
                DataType *_in_p = in_p;
                _in_p -= (params.box + (X * params.box));

                int flag = 1;
                for(size_t l=0;l<params.box_n;l++)
                {
                    for(size_t k=0;k<params.box_n;k++)
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
                    DEBUG_PRINT("Found pixel above thresh at (%ld, %ld)\n", i, j);

                    // We are the highest pixel
                    // Mark the image for duplicates

                    // Mark the center pixel using the bitwise value CENTROIDS_CENT_PIXEL
                    // Increment 1 for each surrounding area
                    *out_p |= CENTROIDS_CENT_PIXEL; 

                    // Increment photon counter
                    n_photons++;

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
                            map_p->x = i + k;
                            map_p->y = j + l;
                            map_p->image = _in_p;
                            map_p->out = _out_p;
                            DEBUG_PRINT("pixel = %d\n", (int)*_in_p);

                            // Increment the out value
                            (*_out_p)++;

                            // Increment pointers
                            _in_p++;
                            _out_p++;
                            map_p++;
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
    return n_photons;
}

// Templates for common datatypes
template int centroids_initialize_params<uint16_t>(centroid_params<uint16_t> &params);
template int centroids_process<uint16_t>(uint16_t *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, size_t N, centroid_params<uint16_t> params);
