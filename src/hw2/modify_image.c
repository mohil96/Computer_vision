#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include <stdlib.h>
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/

    return get_pixel(im, round(x), round(y), c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    float a_w = (float)im.w / (float)w;
    float b_w = 0.5 * (a_w - 1);
    float a_h = (float)im.h / (float)h;
    float b_h = 0.5 * (a_h - 1);
    image resized_im = make_image(w, h, im.c);
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            float x_prime = a_w * x + b_w;
            float y_prime = a_h * y + b_h;
            for (int c = 0; c < resized_im.c; c++)
            {
                float v = nn_interpolate(im, x_prime, y_prime, c);
                set_pixel(resized_im, x, y, c, v);
            }
        }
    }

    return resized_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    int x_h = ceilf(x);
    int x_l = floorf(x);
    int y_h = ceilf(y);
    int y_l = floorf(y);
    float d1 = x - x_l;
    float d2 = x_h - x;
    float d3 = y - y_l;
    float d4 = y_h - y;
    float A1 = d2 * d4;
    float A2 = d1 * d4;
    float A3 = d2 * d3;
    float A4 = d1 * d3;
    float q = get_pixel(im, x_l, y_l, c) * A1 + get_pixel(im, x_h, y_l, c) * A2 + get_pixel(im, x_l, y_h, c) * A3 + get_pixel(im, x_h, y_h, c) * A4;
    return q;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    float a_w = (float)im.w / (float)w;
    float b_w = 0.5 * (a_w - 1);
    float a_h = (float)im.h / (float)h;
    float b_h = 0.5 * (a_h - 1);
    image resized_im = make_image(w, h, im.c);
    for (int x = 0; x < w; x++)
    {
        for (int y = 0; y < h; y++)
        {
            float x_prime = a_w * x + b_w;
            float y_prime = a_h * y + b_h;
            for (int c = 0; c < resized_im.c; c++)
            {
                float v = bilinear_interpolate(im, x_prime, y_prime, c);
                set_pixel(resized_im, x, y, c, v);
            }
        }
    }

    return resized_im;
}

/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    float sum = 0.0;
    for (int i = 0; i < im.w * im.h * im.c; i++)
    {
        sum = im.data[i] + sum;
    }
    for (int i = 0; i < im.w * im.h * im.c; i++)
    {
        if (sum > 0)
        {
            im.data[i] = im.data[i] / sum;
        }
        else
        {
            im.data[i] = 1.0 / im.w / im.h;
        }
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image box_filter = make_image(w, w, 1);
    l1_normalize(box_filter);
    return box_filter;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be
      preserved. Check the detailed algorithm given in the README.
    ************************************************************************/
    assert(im.c == filter.c || filter.c == 1);
    image convolved_image;
    if (preserve == 1)
    {
        convolved_image = make_image(im.w, im.h, im.c);
        for (int x = 0; x < im.w; x++)
        {
            for (int y = 0; y < im.h; y++)
            {
                for (int c = 0; c < im.c; c++)
                {
                    float pixel_value = 0.0;

                    for (int xf = 0; xf < filter.w; xf++)
                    {
                        for (int yf = 0; yf < filter.h; yf++)
                        {
                            if (im.c == filter.c)
                            {
                                pixel_value += get_pixel(filter, xf, yf, c) * get_pixel(im, x - filter.w / 2 + xf, y - filter.h / 2 + yf, c);
                            }
                            else
                            {
                                pixel_value += get_pixel(filter, xf, yf, 0) * get_pixel(im, x - filter.w / 2 + xf, y - filter.h / 2 + yf, c);
                            }
                        }
                    }
                    set_pixel(convolved_image, x, y, c, pixel_value);
                }
            }
        }
    }
    else
    {
        convolved_image = make_image(im.w, im.h, 1);
        for (int x = 0; x < im.w; x++)
        {
            for (int y = 0; y < im.h; y++)
            {
                float pixel_value = 0.0;
                for (int c = 0; c < im.c; c++)
                {
                    for (int xf = 0; xf < filter.w; xf++)
                    {
                        for (int yf = 0; yf < filter.h; yf++)
                        {
                            if (im.c == filter.c)
                            {
                                pixel_value += get_pixel(filter, xf, yf, c) * get_pixel(im, x - filter.w / 2 + xf, y - filter.h / 2 + yf, c);
                            }
                            else
                            {
                                pixel_value += get_pixel(filter, xf, yf, 0) * get_pixel(im, x - filter.w / 2 + xf, y - filter.h / 2 + yf, c);
                            }
                        }
                    }
                }
                set_pixel(convolved_image, x, y, 0, pixel_value);
            }
        }
    }
    return convolved_image;
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image highpass_filter = make_image(3, 3, 1);
    highpass_filter.data[0] = 0;
    highpass_filter.data[1] = -1;
    highpass_filter.data[2] = 0;
    highpass_filter.data[3] = -1;
    highpass_filter.data[4] = 4;
    highpass_filter.data[5] = -1;
    highpass_filter.data[6] = 0;
    highpass_filter.data[7] = -1;
    highpass_filter.data[8] = 0;

    return highpass_filter;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image sharpen_filter = make_image(3, 3, 1);
    sharpen_filter.data[0] = 0;
    sharpen_filter.data[1] = -1;
    sharpen_filter.data[2] = 0;
    sharpen_filter.data[3] = -1;
    sharpen_filter.data[4] = 5;
    sharpen_filter.data[5] = -1;
    sharpen_filter.data[6] = 0;
    sharpen_filter.data[7] = -1;
    sharpen_filter.data[8] = 0;

    return sharpen_filter;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image emboss_filter = make_image(3, 3, 1);
    emboss_filter.data[0] = -2;
    emboss_filter.data[1] = -1;
    emboss_filter.data[2] = 0;
    emboss_filter.data[3] = -1;
    emboss_filter.data[4] = 1;
    emboss_filter.data[5] = 1;
    emboss_filter.data[6] = 0;
    emboss_filter.data[7] = 1;
    emboss_filter.data[8] = 2;

    return emboss_filter;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should use preserve for sharpen and emboss filters as they affect the original colors to either enhance the image or style it.
// High pass filter is used for finding edges and hence does not require 3 color channels.

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: We need to perform post-processing (clamp_image()) on all the filters(high_pass,sharpen and emboss) so that all pixels are normalized between 0 and 1 after convolution.

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    int size = ceilf(sigma * 6);
    if (size % 2 == 0)
    {
        size += 1;
    }
    int offset = (int)(size / 2);
    image gaussian_filter = make_image(size, size, 1);
    for (int x = 0; x < gaussian_filter.w; x++)
    {
        for (int y = 0; y < gaussian_filter.h; y++)
        {

            float pixel_value = 1.0 / TWOPI / pow(sigma, 2) * exp(-(pow(x - offset, 2) + pow(y - offset, 2)) / (2 * pow(sigma, 2)));
            set_pixel(gaussian_filter, x, y, 0, pixel_value);
        }
    }
    l1_normalize(gaussian_filter);
    return gaussian_filter;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.c == b.c && a.w == b.w && a.h == b.h);
    image added_image = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w * a.h * a.c; i++)
    {
        added_image.data[i] = a.data[i] + b.data[i];
    }

    return added_image;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert(a.c == b.c && a.w == b.w && a.h == b.h);
    image subtracted_image = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w * a.h * a.c; i++)
    {
        subtracted_image.data[i] = a.data[i] - b.data[i];
    }

    return subtracted_image;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image gx_filter = make_image(3, 3, 1);
    gx_filter.data[0] = -1;
    gx_filter.data[1] = 0;
    gx_filter.data[2] = 1;
    gx_filter.data[3] = -2;
    gx_filter.data[4] = 0;
    gx_filter.data[5] = 2;
    gx_filter.data[6] = -1;
    gx_filter.data[7] = 0;
    gx_filter.data[8] = 1;

    return gx_filter;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image gy_filter = make_image(3, 3, 1);
    gy_filter.data[0] = -1;
    gy_filter.data[1] = -2;
    gy_filter.data[2] = -1;
    gy_filter.data[3] = 0;
    gy_filter.data[4] = 0;
    gy_filter.data[5] = 0;
    gy_filter.data[6] = 1;
    gy_filter.data[7] = 2;
    gy_filter.data[8] = 1;

    return gy_filter;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
    float min = im.data[0];
    float max = min;
    for (int i = 0; i < im.w * im.h * im.c; i++)
    {
        if (im.data[i] > max)
        {
            max = im.data[i];
        }
        if (im.data[i] < min)
        {
            min = im.data[i];
        }
    }
    if (max - min == 0)
    {
        for (int i = 0; i < im.w * im.h * im.c; i++)
        {
            im.data[i] = 0;
        }
    }
    else
    {
        for (int i = 0; i < im.w * im.h * im.c; i++)
        {
            im.data[i] = (im.data[i] - min) / (max - min);
        }
    }
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    image Gx = convolve_image(im, make_gx_filter(), 0);
    image Gy = convolve_image(im, make_gy_filter(), 0);
    image gradient_magnitude = make_image(im.w, im.h, 1);
    image gradient_direction = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w * im.h; i++)
    {
        gradient_magnitude.data[i] = sqrtf(pow(Gx.data[i], 2) + pow(Gy.data[i], 2));
        gradient_direction.data[i] = atan2f(Gy.data[i], Gx.data[i]);
    }
    sobelimg[0] = gradient_magnitude;
    sobelimg[1] = gradient_direction;
    return sobelimg;
}

image colorize_sobel(image im)
{
    // TODO
    /***********************************************************************
      Create a colorized version of the edges in image "im" using the
      algorithm described in the README.
    ************************************************************************/
    image gaussian = make_gaussian_filter(4);
    image sob_img = convolve_image(im, gaussian, 1);
    image *sobel_img = sobel_image(sob_img);
    image mag = sobel_img[0];
    image dir = sobel_img[1];
    image hsv_img = make_image(im.w, im.h, im.c);
    feature_normalize(mag);
    feature_normalize(dir);
    for (int i = 0; i < im.w * im.h; i++)
    {
        hsv_img.data[i] = dir.data[i];
        hsv_img.data[i + im.w * im.h] = mag.data[i];
        hsv_img.data[i + 2 * im.w * im.h] = mag.data[i];
    }
    hsv_to_rgb(hsv_img);

    return hsv_img;
}

// EXTRA CREDIT: Median filter

int compare(const void *a, const void *b)
{
    float x = *(const float *)a;
    float y = *(const float *)b;
    return (x > y) - (x < y);
}

image apply_median_filter(image im, int kernel_size)
{
    image flitered_image = make_image(im.w, im.h, im.c);
    image median_filter = make_box_filter(kernel_size);
    int edgex = floor(median_filter.w / 2);
    int edgey = floor(median_filter.h / 2);
    for (int c = 0; c < im.c; c++)
    {
        for (int x = edgex; x < im.w - edgex; x++)
        {
            for (int y = edgey; y < im.h - edgey; y++)
            {
                int i = 0;
                for (int xf = 0; xf < median_filter.w; xf++)
                {
                    for (int yf = 0; yf < median_filter.h; yf++)
                    {
                        median_filter.data[i] = get_pixel(im, x + xf - edgex, y + yf - edgey, c);
                        i++;
                    }
                }
                qsort(median_filter.data, median_filter.w * median_filter.h, sizeof(float), compare);
                set_pixel(flitered_image, x, y, c, median_filter.data[median_filter.w * median_filter.h / 2]);
            }
        }
    }
    return flitered_image;
}

// SUPER EXTRA CREDIT: Bilateral filter

image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
    image filtered_image = make_image(im.w, im.h, im.c);
    int kernel_size = 6 * sigma1;
    image gaussian_spatial = make_gaussian_filter(sigma1);
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            for (int c = 0; c < im.c; c++)
            {
                float normalized_sum = 0;
                float convoluted_value = 0;
                for (int i = 0; i < kernel_size; i++)
                {
                    for (int j = 0; j < kernel_size; j++)
                    {
                        int offset_x = x - ((int)(kernel_size / 2) - i);
                        int offset_y = y - ((int)(kernel_size / 2) - j);
                        float gaussian_color = 1.0 / TWOPI / pow(sigma2, 2) * exp(-(pow((im.data[offset_x + offset_y * im.w + c * im.w * im.h] - im.data[x + y * im.w + c * im.w * im.h]), 2)) / pow(sigma2, 2));
                        normalized_sum += gaussian_spatial.data[j + i * kernel_size] * gaussian_color;
                    }
                }
                for (int i = 0; i < kernel_size; i++)
                {
                    for (int j = 0; j < kernel_size; j++)
                    {
                        int offset_x = x - ((int)(kernel_size / 2) - i);
                        int offset_y = y - ((int)(kernel_size / 2) - j);
                        float gaussian_color = 1.0 / TWOPI / pow(sigma2, 2) * exp(-(pow((im.data[offset_x + offset_y * im.w + c * im.w * im.h] - im.data[x + y * im.w + c * im.w * im.h]), 2)) / pow(sigma2, 2));
                        convoluted_value += im.data[offset_x + offset_y * im.w + c * im.w * im.h] * ((gaussian_color * gaussian_spatial.data[j + i * kernel_size]) / normalized_sum);
                    }
                }
                filtered_image.data[x + y * im.w + c * im.w * im.h] = convoluted_value;
            }
        }
    }
    return filtered_image;
}
