---
layout: post
title: Solving the TV-l1 Denoising Model using Chambolle-Pock
short_description: "We will solve the ROF denoising model using Chambolle-Pock method"
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/dual_rof_denoise.png"
tags: chambolle pock tv l1 fenchel dual optimization matlab bilevel toolbox
active: blog
---
Image denoising models typically make assumptions on the recovered image, this is due to the ill-posedness of the problem at hand. It can be seen as an inverse problem, meaning, that there may be different images that can be recovered from a noise contaminated one.
All previous post involving [ROF denoising](https://david.villacis.net/blog/2018/12/14/primal-dual-methods-for-ROF-image-denoising/) made one of such strong assumtions. It assumed the noise in the image to have a *Gaussian* distribution along the image, moreover with the same mean and standar deviation. This assumption is rather strong since most of the noise in the real images come from errors in AD converters and rounding errors from microprocessors. Moreover, we can have an image contaminated with a different type of noise, e.g., *Poisson Noise* and the results will be less than satisfactory.

![poisson](https://im.snibgo.com/ns_ns_toes_nse_Poisson.png)

Indeed, the gaussian noise assumtions are violated and the ROF no longer provides a good solution. In order to solve this problem we have to make a different assumption, for a *Poisson Noise* assumption, thanks to the work by Nikolova et. al. that the best suited model it the TV-$$l_1$$ Image Denoising model.

$$
  J(x) = \|x-f\|_1 + \alpha\|\mathbb{K}x\|_{2,1}
$$
