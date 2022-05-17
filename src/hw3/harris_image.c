#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for (i = 0; i < n; ++i)
    {
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i % im.w;
    d.p.y = i / im.w;
    d.data = calloc(w * w * im.c, sizeof(float));
    d.n = w * w * im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for (c = 0; c < im.c; ++c)
    {
        float cval = im.data[c * im.w * im.h + i];
        for (dx = -w / 2; dx < (w + 1) / 2; ++dx)
        {
            for (dy = -w / 2; dy < (w + 1) / 2; ++dy)
            {
                float val = get_pixel(im, i % im.w + dx, i / im.w + dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for (i = -9; i < 10; ++i)
    {
        set_pixel(im, x + i, y, 0, 1);
        set_pixel(im, x, y + i, 0, 1);
        set_pixel(im, x + i, y, 1, 0);
        set_pixel(im, x, y + i, 1, 0);
        set_pixel(im, x + i, y, 2, 1);
        set_pixel(im, x, y + i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for (i = 0; i < n; ++i)
    {
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: make separable 1d Gaussian.

    int size = ceilf(sigma * 6);
    if (size % 2 == 0)
    {
        size += 1;
    }
    int offset = (int)(size / 2);
    int y = 0;
    image gaussian_filter = make_image(size, 1, 1);
    for (int x = 0; x < gaussian_filter.w; x++)
    {

        float pixel_value = 1.0 / TWOPI / pow(sigma, 2) * exp(-(pow(x - offset, 2) + pow(y - offset, 2)) / (2 * pow(sigma, 2)));
        set_pixel(gaussian_filter, x, y, 0, pixel_value);
    }
    l1_normalize(gaussian_filter);
    return gaussian_filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // TODO: use two convolutions with 1d gaussian filter.
    image gauss_c = make_1d_gaussian(sigma);
    image gauss_r = make_image(1, gauss_c.w, 1);
    for (int y = 0; y < gauss_r.h; y++)
    {

        float pixel_value = get_pixel(gauss_c, y, 0, 0);
        set_pixel(gauss_r, 0, y, 0, pixel_value);
    }
    image g1 = convolve_image(im, gauss_c, 1);

    return convolve_image(g1, gauss_r, 1);
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image Ix = convolve_image(im, make_gx_filter(), 0);
    image Iy = convolve_image(im, make_gy_filter(), 0);
    for (int i = 0; i < im.w * im.h; i++)
    {
        S.data[i] = Ix.data[i] * Ix.data[i];
        S.data[i + im.w * im.h] = Iy.data[i] * Iy.data[i];
        S.data[i + 2 * im.w * im.h] = Ix.data[i] * Iy.data[i];
    }
    S = smooth_image(S, sigma);
    free_image(Ix);
    free_image(Iy);

    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float alpha = 0.06;
    for (int x = 0; x < R.w; x++)
    {
        for (int y = 0; y < R.h; y++)
        {

            float det_value = get_pixel(S, x, y, 0) * get_pixel(S, x, y, 1) - get_pixel(S, x, y, 2) * get_pixel(S, x, y, 2);
            float trace_value = get_pixel(S, x, y, 0) + get_pixel(S, x, y, 1);
            float cornerness_value = det_value - alpha * trace_value * trace_value;
            set_pixel(R, x, y, 0, cornerness_value);
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    for (int x = 0; x < r.w; x++)
    {
        for (int y = 0; y < r.h; y++)
        {
            float center_pixel = get_pixel(r, x, y, 0);
            for (int x_n = -w; x_n <= w; x_n++)
            {
                for (int y_n = -w; y_n <= w; y_n++)
                {
                    float neighbour_pixel = get_pixel(im, x + x_n, y + y_n, 0);
                    if (neighbour_pixel > center_pixel)
                    {
                        set_pixel(r, x, y, 0, -999999);
                        goto end;
                    }
                }
            }
        end:;
        }
    }

    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    // TODO: count number of responses over threshold
    int count = 0; // change this
    for (int i = 0; i < Rnms.w * Rnms.h; i++)
    {
        if (Rnms.data[i] > thresh)
        {
            count++;
        }
    }

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    // TODO: fill in array *d with descriptors of corners, use describe_index.
    int j = 0;
    for (int i = 0; i < Rnms.w * Rnms.h; i++)
    {
        if (Rnms.data[i] > thresh)
        {
            d[j] = describe_index(im, i);
            j++;
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}