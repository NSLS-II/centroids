#include <stdint.h>
#include <stdlib.h>
#include "photons.h"

using namespace std;

/* ----------------------------------------------------------------------------*/
/**
 * \brief Swap the values of two variables
 *
 * @tparam DT
 * \param a Swap variable
 * \param b Swap variable
 */
/* ----------------------------------------------------------------------------*/
template<typename DT> void swap(DT *a, DT *b)
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
void bubble_sort(double *vals, int *x, int *y, int n)
{
    for (int i=0;i<(n-1);i++) 
    {
        for (int j=1;j<(n - i);j++) 
        {
            if (vals[j-1] < vals[j]) 
            {
                swap<double>(&vals[j-1], &vals[j]);
                swap<int>(&x[j-1], &x[j]);
                swap<int>(&y[j-1], &y[j]);
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
template<typename DataType> int process_image(DataType *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, centroid_params<DataType> params)
{
	find_photons<DataType>(image, out, X, Y, params);
    process_bias<DataType>(image, out, X, Y, bias);
    return process_photons<DataType>(image, out, table, X, Y, bias, params);
}

template<typename DataType> int process_bias(DataType *pixels, uint16_t *out, size_t X, size_t Y, double *bias)
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

int process_pixel(double *pixels, int *x, int *y, int n, double *com_x, double *com_y, int *sum)
{
    double *pixels_p = pixels;
    int *x_p = x;
    int *y_p = y;

    *com_x = 0;
    *com_y = 0;
    *sum = 0;

    for(int i=0;i<n;i++)
    {
        *sum += (*pixels_p); 
        *com_x += (*pixels_p) * *(x_p++);
        *com_y += (*pixels_p) * *(y_p++);
        pixels_p++;
    }

    *com_y /= *sum;
    *com_x /= *sum;
    
    return 0;
}

template<typename DataType> int process_photons(DataType *image, uint16_t *out, double *table, 
        size_t X, size_t Y, double *bias, centroid_params<DataType> params)
{
    int box_val = params.box;
    int box_val_n = (box_val * 2) + 1;
    int box_val_t = box_val_n * box_val_n;

    int xvals[box_val_t];
    int yvals[box_val_t];
    double pixel_cluster[box_val_t];
    
    DataType *image_p = image;
    uint16_t *out_p = out;
    double *table_p = table;
    double *bias_p = bias;
    int table_n = 0;

    for(size_t j=0;j<Y;j++){
        for(size_t i=0;i<X;i++){
            if((*out_p) & CENT_PIXEL)
            {
                // This is the central hot pixel
                DataType *_ip = image_p;
                _ip -= box_val;
                _ip -= (X * box_val);

                uint16_t *_op = out_p;
                _op -= box_val;
                _op -= (X * box_val);

                double *_pix = pixel_cluster;
                int *_xvals = xvals;
                int *_yvals = yvals;

                bool flag = true;

                uint16_t box_sum = 0;
                for(int l=-box_val;l<=box_val;l++)
                {
                    for(int k=-box_val;k<=box_val;k++)
                    {
                        box_sum += (*(_op++) & 0x7FFF);
                        *(_pix++) = *(_ip++) - *bias_p;
                        *(_xvals++) = k;
                        *(_yvals++) = l;
                    }
                    _ip += (X - box_val_n);
                    _op += (X - box_val_n);
                }

                if(1) // setup for choosing if to process
                {
                    // We have a valid pixel.

                    double comx, comy;
                    int sum;

                    bubble_sort(pixel_cluster, xvals, yvals, box_val_t);
                    process_pixel(pixel_cluster, xvals, yvals, box_val_t, &comx, &comy, &sum);
    
                    table_n++;
                    table_p[0] = i;  
                    table_p[1] = j;  
                    table_p[2] = sum;
                    table_p[3] = comx + i;
                    table_p[4] = comy + j;
                    table_p[5] = box_sum - box_val_t;
                    for(int n=0;n<9;n++)
                    {
                        table_p[6] += pixel_cluster[n];
                    }

                    for(int n=0;n<box_val_t;n++)
                    {
                        table_p[7 + n] = pixel_cluster[n];
                    }

                    table_p += box_val_t + 7;
                }
            }

            image_p++;
            out_p++;
            bias_p++;
        } // fox X
    } // for Y

    return table_n;
}

template<typename DataType> int find_photons(DataType *image, uint16_t *out, size_t X, size_t Y, 
        centroid_params<DataType> params)
{
    DataType threshold = params.threshold;
    int box_val = params.box;
    int box_val_n = (box_val * 2) + 1;

    uint16_t *out_p = out;
    DataType *in_p = image;

    // Blank the output
    for(size_t i=0;i<(X* Y);i++)
    {
        out[i] = 0;
    }

    // Start loop through image, we need to skip by box_val
    for(size_t j=0;j<box_val;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for(size_t j=box_val;j<(Y-box_val);j++){
        for(size_t i=0;i<box_val;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
        for(size_t i=box_val;i<(X-box_val);i++)
        {
            // This is the main routine. 
            // Is this pixel above threshold
            if(*in_p >= threshold)
            {
                // Check if this is the highest pixel
                // Rewind by 1 box_val in X and 1 box_val in Y
                DataType *_p;
                
                _p = in_p;
                _p -= box_val;
                _p -= (X * box_val);

                int flag = 1;
                for(size_t l=0;l<box_val_n;l++)
                {
                    for(size_t k=0;k<box_val_n;k++)
                    {
                        if(*_p > *in_p)
                        {
                            // We are not the highest pixel, set the flag
                            flag = 0;
                        }

                        _p++;
                    }
                    _p += (X - box_val_n);
                }

                if(flag)
                {
                    // Mark the image for duplicates
                    
                    // Mark the center pixel using the bitwise value CENT_PIXEL
                    // Increment 1 for each surrounding area
                    if (!(*out_p & CENT_PIXEL)){ 
                        *out_p = CENT_PIXEL; 
                        _p = out_p;
                        _p -= box_val;
                        _p -= (X * box_val);

                        for(size_t l=0;l<box_val_n;l++)
                        {
                            for(size_t k=0;k<box_val_n;k++)
                            {
                                (*_p)++;
                                _p++;
                            }
                            _p += (X - box_val_n);
                        }
                    }
                }
            }

            out_p++; in_p++;

        }
        for(size_t i=(X-box_val);i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }

    for(size_t j=Y-box_val;j<Y;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = 1;
            in_p++;
        }
    }
	return 0;
}

// Templates for common datatypes
template int process_image<uint16_t>(uint16_t *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, centroid_params<uint16_t> params);
