#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    assert(c >= 0 && c < im.c);
    if (x < 0)
    {
        x = 0;
    }
    if (x >= im.w)
    {
        x = im.w - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= im.h)
    {
        y = im.h - 1;
    }

    return im.data[im.h * im.w * c + y * im.w + x];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    if (x < 0 || y < 0 || c < 0 || x >= im.w || y >= im.h || c >= im.c)
    {
        return;
    }
    im.data[im.h * im.w * c + y * im.w + x] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float gray_pixel = 0.299 * get_pixel(im, x, y, 0) + 0.587 * get_pixel(im, x, y, 1) + 0.114 * get_pixel(im, x, y, 2);
            set_pixel(gray, x, y, 0, gray_pixel);
        }
    }

    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float shift_pixel = get_pixel(im, x, y, c);
            set_pixel(im, x, y, c, shift_pixel + v);
        }
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w * im.h * im.c; i++)
    {
        if (im.data[i] < 0.0)
        {
            im.data[i] = 0.0;
        }
        else if (im.data[i] > 1.0)
        {
            im.data[i] = 1.0;
        }
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float R = get_pixel(im, x, y, 0);
            float G = get_pixel(im, x, y, 1);
            float B = get_pixel(im, x, y, 2);
            float V = three_way_max(R, G, B);
            float m = three_way_min(R, G, B);
            float C = V - m;
            float S = (V == 0) ? 0 : (C / V);
            float H_dash, H;
            if (C != 0)
            {
                if (V == R)
                {
                    H_dash = (G - B) / C;
                }
                else if (V == G)
                {
                    H_dash = (B - R) / C + 2;
                }
                else if (V == B)
                {
                    H_dash = (R - G) / C + 4;
                }
                H = H_dash < 0 ? (H_dash / 6 + 1) : H_dash / 6;
            }
            else
            {
                H = 0;
            }

            set_pixel(im, x, y, 0, H);
            set_pixel(im, x, y, 1, S);
            set_pixel(im, x, y, 2, V);
        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float H = get_pixel(im, x, y, 0) * 6.0;
            float S = get_pixel(im, x, y, 1);
            float V = get_pixel(im, x, y, 2);
            int Hi = floor(H);
            float F = H - Hi;
            float P = V * (1 - S);
            float Q = V * (1 - F * S);
            float T = V * (1 - (1 - F) * S);

            if (Hi == 0)
            {
                set_pixel(im, x, y, 0, V);
                set_pixel(im, x, y, 1, T);
                set_pixel(im, x, y, 2, P);
            }
            else if (Hi == 1)
            {
                set_pixel(im, x, y, 0, Q);
                set_pixel(im, x, y, 1, V);
                set_pixel(im, x, y, 2, P);
            }
            else if (Hi == 2)
            {
                set_pixel(im, x, y, 0, P);
                set_pixel(im, x, y, 1, V);
                set_pixel(im, x, y, 2, T);
            }
            else if (Hi == 3)
            {
                set_pixel(im, x, y, 0, P);
                set_pixel(im, x, y, 1, Q);
                set_pixel(im, x, y, 2, V);
            }
            else if (Hi == 4)
            {
                set_pixel(im, x, y, 0, T);
                set_pixel(im, x, y, 1, P);
                set_pixel(im, x, y, 2, V);
            }
            else if (Hi == 5)
            {
                set_pixel(im, x, y, 0, V);
                set_pixel(im, x, y, 1, P);
                set_pixel(im, x, y, 2, Q);
            }
        }
    }
}

void scale_image(image im, int c, float v)
{
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float scale_pixel = get_pixel(im, x, y, c) * v;
            set_pixel(im, x, y, c, scale_pixel);
        }
    }
}
void rgb_to_hcl(image im)

{
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float R = get_pixel(im, x, y, 0) / 255;
            float G = get_pixel(im, x, y, 1) / 255;
            float B = get_pixel(im, x, y, 2) / 255;

            float R_prime = (R <= 0.04045) ? R / 12.92 : pow((R + 0.055) / 1.055, 2.4);
            float G_prime = (G <= 0.04045) ? G / 12.92 : pow((G + 0.055) / 1.055, 2.4);
            float B_prime = (B <= 0.04045) ? B / 12.92 : pow((B + 0.055) / 1.055, 2.4);

            float X = 0.4124564 * R_prime + 0.3575761 * G_prime + 0.1804375 * B_prime;
            float Y = 0.2126729 * R_prime + 0.7151522 * G_prime + 0.0721750 * B_prime;
            float Z = 0.0193339 * R_prime + 0.1191920 * G_prime + 0.9503041 * B_prime;

            float Yn = 1.0;

            float un = 0.2009;
            float vn = 0.4610;
            float u = 4 * X / (X + 15 * Y + 3 * Z);
            float v = 9 * Y / (X + 15 * Y + 3 * Z);

            float L = Y / Yn <= pow(6 / 29, 3) ? pow(29 / 3, 3) * Y / Yn : 116 * pow((Y / Yn), (1 / 3)) - 16;
            float U = 13 * L * (u - un);
            float V = 13 * L * (v - vn);

            float C = sqrt(pow(U, 2) + pow(V, 2));
            float H = atan2(V, U);

            set_pixel(im, x, y, 0, H);
            set_pixel(im, x, y, 1, C);
            set_pixel(im, x, y, 2, L);
        }
    }
}
void hcl_to_rgb(image im)
{
    for (int x = 0; x < im.w; x++)
    {
        for (int y = 0; y < im.h; y++)
        {
            float H = get_pixel(im, x, y, 0);
            float C = get_pixel(im, x, y, 1);
            float L = get_pixel(im, x, y, 2);

            float U = C * cos(H);
            float V = C * sin(H);
            float un = 0.2009;
            float vn = 0.4610;
            float Yn = 1.0;
            float u = U / 13 / L + un;
            float v = V / 13 / L + vn;

            float Y = L <= 8 ? Yn * L * pow(3 / 29, 3) : Yn * pow((L + 16) / 116, 3);
            float X = Y * 9 * U / 4 / V;
            float Z = Y * (12 - 3 * U - 20 * V) / 4 / V;

            float R_prime = 3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z;
            float G_prime = -0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z;
            float B_prime = 0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z;

            float R = R_prime <= 0.0031308 ? R_prime * 12.92 : 1.055 * pow(R_prime, 1 / 2.4) - 0.055;
            float G = G_prime <= 0.0031308 ? G_prime * 12.92 : 1.055 * pow(G_prime, 1 / 2.4) - 0.055;
            float B = B_prime <= 0.0031308 ? B_prime * 12.92 : 1.055 * pow(B_prime, 1 / 2.4) - 0.055;

            set_pixel(im, x, y, 0, R * 255);
            set_pixel(im, x, y, 1, G * 255);
            set_pixel(im, x, y, 2, B * 255);
        }
    }
}