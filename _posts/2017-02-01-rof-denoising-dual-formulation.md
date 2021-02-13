---
layout: post
title: Dual ROF Image Denoising Model Formulation
short_description: We will propose a dual model formulation for the Rudin Osher Fatemi denoising model
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/dual_rof_denoise.png"
tags: Rudin Osher Fatemi dual denoising variational models fista forward backward
active: blog
---
In this post we will review a technique for solving the image denoising problem using the method proposed by Rudin, Osher and Fatemi also known as the ROF denoising model, more details on the model can be found [here](http://www.math-info.univ-paris5.fr/~lomn/Cours/ECE/PhysicaRudinOsher.pdf). In this case we asume that the image is corrupted by gaussian distributed noise, and if we assume the original image is a matrix of size $$n_1 \times n_2$$ in what follows we will treat this as a vector $$u \in \mathbb{R}^{n=n_1\times n_2}$$. Therefore the ROF primal formulation is as follows:

$$\min_{u\in \mathbb{R}^n} \frac{1}{2}\| u-f \|_2^2 + \lambda \| \mathbb{K}u \|_{2,1}$$

where $$f \in \mathbb{R}^m$$ is the corrupted image, $$\mathbb{K}$$ is a linear operator (in this case the gradient), $$\lambda$$ is a Tickonov Regularization term and  $$\| . \|_{2,1}$$ is the *Isotropic Total Variation* seminorm defined as:

$$
\| u \|_{2,1} = \sum_{j=1}^n \sqrt{u_j^2 + u_{n+j}^2}
$$

for a vector $$u\in\mathbb{R}^{2n}$$.

Several techniques have used to perform this optimization problem, the main challenge in facing this kind of optimization problems is the treatment of the non-differentiable Total Variation (TV) seminorm. Most of the available literature uses a smoothed version of this norm (Huber, Berkovier-Engelman regularizations) in order to obtain the optimality conditions and to characterize a gradient.

The approach that we will be taking is to formulate the Fechel ROF Dual problem, with the goal of having a formulation that doesn't require any smoothing.

## The Fenchel Dual Formulation

### Convex Conjugate
Let's recall the definition of a convex conjugate: Given a function $$f:\mathbb{R} \to \mathbb{R}$$ a general possibly non-convex function, we define the *convex conjugate* as:


$$
f^* (p) = sup_{u\in\mathbb{R}^n}(\langle u,p \rangle - f(u))
$$

### ROF Fechel Dual Formulation
If we separate the ROF functional in two terms: $$f(\mathbb{K}u) = \| u \|_{2,1}$$ and $$g(u) = \frac{1}{2}\| u-f \|_2^2$$, we can find the dual formulation of a problem of the type: $$\min_u f(\mathbb{K}u) + g(u)$$

$$
 \begin{align}
 \min_{u\in \mathbb{R}^n} f(\mathbb{K}x) + g(x) &= \min_{u\in \mathbb{R}^n} \sup_{p \in \mathbb{R}^m} \langle p, \mathbb{K}u \rangle - f^* (p)+ g(u)\\
 &= \max_{p \in \mathbb{R}^m} \inf_u  \langle p, \mathbb{K}u \rangle - f^* (p) + g(u)\\
 &= \max_{p \in \mathbb{R}^m} -f^* (p) - g^* (\mathbb{K}^* p)\\
 &= \min_{p \in \mathbb{R}^m} f^* (p) + g^* (-\mathbb{K}^* p)
 \end{align}
$$

#### Dual of the function $$g(u) = \frac{1}{2}\|u-f\|_2^2 $$

For this function we will apply the definition of the fenchel dual:


$$

  \begin{align}
  g(u) &= \frac{1}{2}\|u-f\|_2^2\\
  g^*(p) &= \sup_u \left\{\langle p,u \rangle - \frac{1}{2}\|u-f\|_2^2 \right\}\\
  \end{align}
$$

In order to find the supremum of such argument we find its critical point:


$$

  \begin{align}
  \partial_u \left( \langle p,u \rangle - \frac{1}{2} \|u-f\|_2^2 \right) &= 0\\
  p -(u-f) &= 0\\
  \end{align}
$$

From this optimality condition we can define the optimal point characterized as:


$$

  \begin{align}
  p &= u-f\\
  u &= p + f
  \end{align}
$$

This point can be used as the point where the function reaches a maximum, then we can describe the dual funtion $$g^*$$ as:


$$

  \begin{align}
  g^* (p) &= \langle p, p+f \rangle - \frac{1}{2} \|p+f-f\|_2^2,\\
  g^* (p)&= \langle p,p \rangle + \langle p,f \rangle - \frac{1}{2}\|p\|_2^2\\
  g^* (p)&= \frac{1}{2}\|p\|_2^2 + \langle p,f \rangle
  \end{align}
$$

#### Dual of the function $$f(u) = \lambda \| u \|_{2,1} $$

Applying the convex conjugate definition, we can see that $$f^* (p)$$ is the indicator function of the pointwise two-dimensional unit ball, given that $$p\in\mathbb{R}^{2m}$$:


$$
f^* (p) =
\begin{cases}
0,\; \max_{j=1,\dots,m} \sqrt{p_j^2 + p_{m+j}^2} \le \lambda,\\
\infty, \; \text{otherwise}
\end{cases}
$$

#### Placing it all together

Now we will use the dual functions calculated previously to get a representation of the *Dual ROF Model*, from the Fechel-Rockafellar duality we know that:


$$

  \begin{align}
  \min_{p \in \mathbb{R}^m}& f^* (p) + g^* (-\mathbb{K}^* p)\\
  \min_{p \in \mathbb{R}^m}& f^* (p) + \frac{1}{2} \|\mathbb{K}^* p\|_2^2 - \langle \mathbb{K}^* p, f \rangle\\
  \min& \left\{ \frac{1}{2} \|\mathbb{K}^* p\|_2^2 - \langle \mathbb{K}^* p, f \rangle : \sqrt{p_j^2 + p_{m+j}^2} \le \lambda \; \forall j = 1,\dots,m\right\}
  \end{align}
$$
