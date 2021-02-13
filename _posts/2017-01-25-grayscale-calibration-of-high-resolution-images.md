---
layout: post
title: Grayscale Calibration of High Resolution Images
short_description: This is a description of the process of intensity calibration for the playing cards open dataset
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/gray_calibration.png"
tags: high resolution grayscale intensity calibration
active: blog
---

## Introduction
In this post I will describe in detail, the steps taken for building the Playing Cards Photographic Dataset that can be downloaded in this [link](http://www.fips.fi/photographic_dataset.php).

This dataset was build to provide researchers in the Image Processing and Inverse Problems communities a set of high resolution images corrupted with real noise. Traditionally denoising models have been tested against synthetic noise created with a known distribution which in my opinion can bias the performance of certain models. We chose lpaying cards given the fact that the strong discontinuities in the colors can be ideal for testing Total Variation related models, which will be described in later posts.

This dataset was built at the Industrial Mathematics Laboratory of the Department of Mathematics and Statistics of the University of Helsinki in colaboration with the Center of Mathematical Modelling (MODEMAT) at Escuela Politécnica Nacional in Quito, Ecuador.

### Camera Equipment
We used a PhaseOne XF medium-format camera equipped with an achromatic IQ260 digital block. The lens is Phase One Digital AF 120mm F4. The pixel size in the resulting 16bit TIFF image is 8964 x 6716.

### Lighting
The targets were lit with ve Olight X6 Marauder LED flashlights with luminous flux of 5000 lm. The lights were positioned at roughly equiangular arrangement. The distance of each light from the target was roughly 850 mm.
The lights were heating up quickly as they were used at maximum power. Cooling was enhanced with three regular household fans. A diuser was placed between the lights and the target to make the lighting more uniform and to reduce sharp shadows.

### Details of target
We placed playing cards on a horizontal surface. The camera was aimed directly down, so the optical axis was roughly vertical. In every picture there is a five-step grayscale target for calibration.

### Varying the noise level in the data
For each arrangement of cards, four different images were taken:


1. Best-quality photo (im clean.tif). Minimum ISO setting 200 and histogram approximately spanning the full dynamic range. The noise level is very small.
2. Medium-quality photo, rst type (im noise1.tif). Minimum ISO setting 200 and histogram approximately spanning a quarter the full dynamic range. This will introduce round-o errors from the AD converter.
3. Medium-quality photo, second type (im noise2.tif). Maximum ISO setting 3200 and histogram approximately spanning the full dynamic range. This will introduce photon-counting noise as well as some electronic noise.
4. Low-quality photo (im noise3.tif). Maximum ISO setting 3200 and histogram approximately spanning a quarter of the full dynamic range. The noise will be increased by added round-off errors.

<div id="portfolio" class="portfolio grid-container clearfix">

						<article class="portfolio-item pf-media pf-icons">
							<div class="portfolio-image">
								<a href="portfolio-single.html">
									<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/im_clean.png" alt="Clean Image">
								</a>
							</div>
							<div class="portfolio-desc">
								<h3><a>Clean Image</a></h3>
							</div>
						</article>

						<article class="portfolio-item pf-media pf-icons">
							<div class="portfolio-image">
								<a>
									<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/im_noise1.png" alt="Noise Type 1">
								</a>
							</div>
							<div class="portfolio-desc">
								<h3><a>Noise Type 1</a></h3>
							</div>
						</article>

						<article class="portfolio-item pf-media pf-icons">
							<div class="portfolio-image">
								<a>
									<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/im_noise2.png" alt="Noise Type 2">
								</a>
							</div>
							<div class="portfolio-desc">
								<h3><a>Noise Type 2</a></h3>
							</div>
						</article>

						<article class="portfolio-item pf-media pf-icons">
							<div class="portfolio-image">
								<a>
									<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/im_noise3.png" alt="Noise Type 3">
								</a>
							</div>
							<div class="portfolio-desc">
								<h3><a>Noise Type 3</a></h3>
							</div>
						</article>
</div>

## Intensity Calibration Process
As it can be seen in the previous Figure, we had 4 different images with different intensity values, the goal of the calibration process is to generate a new set of images based that span the complete dynamic range and preserves all the non-linear properties presented during the aquisition process. The entire process can be summarized as follows:


1. *Linear Intensity Adjustment*: We perform a linear adjustment so that all the images span all the dynamic range.
2. *Grayscale Non Linear Calibration*: This process uses the 5-level grayscale phantoms included in the image to include the nonlinearity generated during the aquisition process.

### Linear Intensity Adjustment
This process takes all of the 4 images, drops 1% of the highest and lowest values presented in the histogram, and linearly expands the histogram so it covers the full dynamic range $$[0-65535]$$.

<pre class="matlab-code">
%% Linear adjustment to make the image span all the dynamic range
function [ adjusted_image ] = image_adjust_linear( img, lowIn, highIn, lowOut, highOut )
  adjusted_image = (img < lowIn) .* lowOut;
  adjusted_image = adjusted_image + (img >= lowIn & img < highIn) .* ...
      (lowOut + (highOut - lowOut) .* ((img-lowIn)/(highIn - lowIn)));
  adjusted_image = adjusted_image + (img >= highIn).* highOut;
</pre>

### Grayscale Non-linear Calibration
In this step we take the information extracted from the 5-level grayscale phantom provided within the image to set a vector containing the corresponding intensity for each level.

<pre class="matlab-code">
%% Averaging the grayleves in each patch
grayLevels = [round(mean(patch1(:))),round(mean(patch2(:))),...
    round(mean(patch3(:))),round(mean(patch4(:))),round(mean(patch5(:)))];
</pre>

Using these levels from the corrupted image and from the clean image, we perform a spline interpolation in such a way that the non-linearities presenten in the clean image are preserved during the calibration process.

<pre class="matlab-code">
%% Spline correlation between the source and target grayscale
function [ adjusted_image ] = image_adjust_spline( img, grayLevels_src, grayLevels_ref )
    xx = 0:1:2^16;
    yy = spline(grayLevels_src,grayLevels_ref,xx);
    adjusted_image = round(yy(round(img+1)));
end
</pre>

This way we end up with the following calibrated images:

<div id="portfolio" class="portfolio grid-container clearfix">
	<article class="portfolio-item pf-media pf-icons">
		<div class="portfolio-image">
			<a href="portfolio-single.html">
				<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/corrected_im_clean.png" alt="Clean Image">
			</a>
		</div>
		<div class="portfolio-desc">
			<h3><a>Clean Image</a></h3>
		</div>
	</article>

	<article class="portfolio-item pf-media pf-icons">
		<div class="portfolio-image">
			<a>
				<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/corrected_im_noise1.png" alt="Noise Type 1">
			</a>
		</div>
		<div class="portfolio-desc">
			<h3><a>Noise Type 1</a></h3>
		</div>
	</article>

	<article class="portfolio-item pf-media pf-icons">
		<div class="portfolio-image">
			<a>
				<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/corrected_im_noise2.png" alt="Noise Type 2">
			</a>
		</div>
		<div class="portfolio-desc">
			<h3><a>Noise Type 2</a></h3>
		</div>
	</article>

	<article class="portfolio-item pf-media pf-icons">
		<div class="portfolio-image">
			<a>
				<img src="https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/grayscale_calibration/corrected_im_noise3.png" alt="Noise Type 3">
			</a>
		</div>
		<div class="portfolio-desc">
			<h3><a>Noise Type 3</a></h3>
		</div>
	</article>
</div>
