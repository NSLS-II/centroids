#include <stdint.h>
#include <stdlib.h>
#include "photons.h"

using namespace std;

int find_photons_uint16(uint16_t *image, uint16_t *out, double *table, double *bias,
        size_t X, size_t Y, uint16_t threshold)
{
	find_photons<uint16_t>(image, out, X, Y, threshold);
    process_bias<uint16_t>(image, out, X, Y, bias);
    process_photons<uint16_t>(image, out, table, X, Y, bias);

    return 0;
}

template<typename DataType> int process_bias(DataType *pixels, DataType *out, size_t X, size_t Y, double *bias)
{
    DataType *pix_p = pixels;
    DataType *out_p = out;
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

int process_pixel(double *pixels, int box_val, double *com_x, double *com_y, int *sum)
{
    *com_x = 0;
    *com_y = 0;
    *sum = 0;

    double *pixels_p = pixels;

    for(int j=-box_val;j<=box_val;j++)
    {
        for(int i=-box_val;i<=box_val;i++)
        {
           *sum += (*pixels_p); 
           *com_x += (*pixels_p) * i;
           *com_y += (*pixels_p) * j;
           pixels_p++;
        }
    }

    *com_y /= *sum;
    *com_x /= *sum;
    
    return 0;
}

template<typename DataType> int process_photons(DataType *image, DataType *out, double *table, 
        size_t X, size_t Y, double *bias)
{
    // We habe found the highest pixel.
    // Now Process it
    
    unsigned int box_val = 1;
    unsigned int box_val_n = (box_val * 2) + 1;

    double *pixel_cluster = (double*)malloc(sizeof(double) * box_val_n * box_val_n);
    if(!pixel_cluster)
    {
        return 127;
    }
    
    // Loop over all photons

    DataType *image_p = image;
    DataType *out_p = out;
    double *table_p = table;
    double *bias_p = bias;

    for(size_t j=0;j<Y;j++){
        for(size_t i=0;i<X;i++){
            if((*out_p) & 0x100)
            {
                // This is the central hot pixel
                DataType *_ip = image_p;
                _ip -= box_val;
                _ip -= (X * box_val);

                DataType *_op = out_p;
                _op -= box_val;
                _op -= (X * box_val);

                double *_pix = pixel_cluster;

                bool flag = true;

                for(size_t l=0;l<box_val_n;l++)
                {
                    for(size_t k=0;k<box_val_n;k++)
                    {
                        if((*_op & 0xFF) != 1)
                        {
                            flag = false;
                        }
                        *_pix = *_ip - *bias_p;
                        _ip++;
                        _op++;
                        _pix++;
                    }
                    _ip += (X - box_val_n);
                    _op += (X - box_val_n);
                }

                if(flag)
                {
                    // We hve a valid pixel.

                    double comx, comy;
                    int sum;

                    process_pixel(pixel_cluster, box_val, &comx, &comy, &sum);

                    table_p[0] = i;  
                    table_p[1] = j;  
                    table_p[2] = sum;
                    table_p[3] = comx + i;
                    table_p[4] = comy + j;
                    table_p[5] = 0;
                    table_p += 6; 
                }
            }

            image_p++;
            out_p++;
            bias_p++;
        } // fox X
    }

    return 0;
}

template<typename DataType> int find_photons(DataType *image, DataType *out, size_t X, size_t Y, DataType threshold)
{
	DEBUG_COMMENT("Hello world...\n");

    unsigned int box_val = 1; // This is a 3 x 3 grid
    unsigned int box_val_n = (box_val * 2) + 1;

    DataType *out_p = out;
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
            *(out_p++) = (DataType)1;
            in_p++;
        }
    }

    for(size_t j=box_val;j<(Y-box_val);j++){
        for(size_t i=0;i<box_val;i++)
        {
            *(out_p++) = (DataType)1;
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
                    
                    // Mark the center pixel using the bitwise value 0x100
                    // Increment 1 for each surrounding area
                    if (!(*out_p && 0x100)){ 
                        *out_p = 0x100; 
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
            *(out_p++) = (DataType)1;
            in_p++;
        }
    }

    for(size_t j=Y-box_val;j<Y;j++)
    {
        for(size_t i=0;i<X;i++)
        {
            *(out_p++) = (DataType)1;
            in_p++;
        }
    }
	return 0;
}

