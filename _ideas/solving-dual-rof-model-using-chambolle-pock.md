---
layout: post
title: Solving the Dual ROF Denoising Model using the Chambolle Pock Algorithm
short_description: "We will solve the dual ROF model using the chambolle pock method"
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/dual_rof_denoise.png"
tags: chambolle pock backward rudin osher fatemi rof fenchel saddle point optimization python numpy scipy bilevel imaging toolbox
active: blog
---
# Solving the Dual ROF Denoising Model using the Chambolle Pock Algorithm
The Rudin-Osher-Fatemi (ROF) Model presented in ... can be solved using a saddle point formulation which can be written as:

$$ \min_x \max_y g(x) + \langle \mathbb{K} x,y\rangle - f^*(y) $$

$$ \min_x \max_y \frac{1}{2}\|x-f\|_2^2 + \langle \mathbb{K} x,y\rangle - \delta_{B_\lambda}(y) $$

This models presents a convenient separation between a differentiable term and a non-differentiable term $$\delta_{B_\lambda}(y)$$. This separation allows us to use a Forward-Barckward numerical approach to solve this problem.

## The Proximal Operator
The proximal operator is a key piece when dealing with non-differentiable functions. It consists on modifying the original function with a Moreau-Yosida regularization:

$$ f_{MY}(x) = f(x)+\frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

This regularization makes the original function to be strictly convex, which implies that it has a unique minimizer (this is not the case for convex functions). Indeed, this minimizer is known as the *proximal* of the function $$f$$:

$$ prox_{\tau\partial f}(x) = arg\min_{\bar{x}} f(x) + \frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

## Chambolle-Pock Method
This is a first order method, it means that it only uses information for the first derivative of the objective function. This method exploits the differentiability of the regular part of the function and the proximal operator of the non-differentiable part.

This method defines the following iterative scheme:

$$ x_{k+1} = prox_{\tau \partial g}(x_k-\tau\mathbb{K}^\top y_k)$$

$$ \bar{x}_{k+1} = 2x_{k+1} - x_k $$

$$ y_{k+1} = prox_{\sigma \partial f^*}(y_k+\sigma \mathbb{K}\bar{x}_{k+1}) $$

Therefore, for our specific case we need to calculate both the proximal of $$g$$ and $$f^*$$:

#### Proximal $$g(x) = \frac{1}{2}\|x-f\|_2^2$$
In this calculation we will use as described in the definition, a Moreau-Yosida regularization of $$g$$:

$$
prox_{\tau\partial g}(x) = arg\min_{\bar x} \frac{1}{2}\|x-f\|_2^2 + \frac{1}{2\tau}\|\bar{x}-x\|_2^2
$$

The necessary and sufficient condition for the optimal value $$\bar{x}$$ is:

$$
(\bar{x}-f)+\frac{1}{\tau}(\bar{x}-x) = 0
$$

$$
x = \tau(\bar{x}-f) + \bar{x}
$$

$$
x = (\tau+1)\bar{x}-\tau f
$$

$$
\bar{x} = \frac{x+\tau f}{\tau+1}
$$

Which yields:

$$
prox_{\tau\partial g}(x) = \frac{x_k + \tau(f-\mathbb{K}^\top y_k)}{\tau +1}
$$

#### Proximal $$f^*(y) = \delta_{B_\lambda}(y)$$

$$
prox_{\sigma\partial f^*}(x) = arg\min_{\bar{x}} \delta_{B_\lambda}(y) + \frac{1}{2\sigma}\|\bar{y}-y\|_2^2
$$

As can be seen this leads to a new optimization problem, this problem can be solved by analyzing the first order optimality condition of this problem:

$$
0 \in \partial \delta_{B_\lambda}(y) + (\bar{y}-y)
$$

Where $$ \partial \delta_{B_\lambda}(y) $$ is the convex subdifferendial of the indicator function of the $$\lambda$$-ball.
