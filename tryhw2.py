from uwimg import *

im = load_image("data/dogsmall.jpg")
a = nn_resize(im, im.w * 4, im.h * 4)
save_image(a, "dog4x-nn")

im = load_image("data/dogsmall.jpg")
a = bilinear_resize(im, im.w * 4, im.h * 4)
save_image(a, "dog4x-bl")

im = load_image("data/dog.jpg")
a = nn_resize(im, im.w // 7, im.h // 7)
save_image(a, "dog7th-bl")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w // 7, blur.h // 7)
save_image(thumb, "dogthumb")

im = load_image("data/dog.jpg")
f = make_highpass_filter()
blur = convolve_image(im, f, 0)
post = clamp_image(blur)
save_image(blur, "dog-highpass")

im = load_image("data/dog.jpg")
f = make_sharpen_filter()
blur = convolve_image(im, f, 1)
post = clamp_image(blur)
save_image(blur, "dog-sharpen")

im = load_image("data/dog.jpg")
f = make_emboss_filter()
blur = convolve_image(im, f, 1)
post = clamp_image(blur)
save_image(blur, "dog-emboss")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")

im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
save_image(mag, "magnitude")

im = load_image("data/dog.jpg")
color = colorize_sobel(im)
save_image(color, "dog_color")

im = load_image("figs/salt_petter_building.jpg")
f = apply_median_filter(im, 3)
save_image(f, "median2")

im = load_image("data/landscape.jpg")
f = apply_bilateral_filter(im, 2.5, 0.32)
save_image(f, "bilateral")
