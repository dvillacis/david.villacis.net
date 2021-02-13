---
layout: post
title: Solving the Dual ROF Denoising Model using a Forward-Backward Algorithm
short_description: "We will solve the dual ROF model using the forward backward method"
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/dual_rof_denoise.png"
tags: forward backward rudin osher fatemi rof fenchel dual optimization python numpy scipy coitoolbox
active: blog
---
# Solving the Dual ROF Denoising Model using the Forward-Backward Algorithm
The Rudin-Osher-Fatemi (ROF) Model presented in ... can be solved using a Fenchel-Rockafellar version which can be written as:

$$ \min_{y\in\mathbb{R}^m} \|f-\mathbb{K}^\top y\|_2^2 + \delta_{B_\lambda}(y) $$

This models presents a convenient separation between a differentiable term and a non-differentiable term $$\delta_{B_\lambda}(y)$$. This separation allows us to use a Forward-Barckward numerical approach to solve this problem.

## The Proximal Operator
The proximal operator is a key piece when dealing with non-differentiable functions. It consists on modifying the original function with a Moreau-Yosida regularization:

$$ f_{MY}(x) = f(x)+\frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

This regularization makes the original function to be strictly convex, which implies that it has a unique minimizer (this is not the case for convex functions). Indeed, this minimizer is known as the *proximal* of the function $$f$$:

$$ prox_{\tau\partial f}(x) = arg\min_{\bar{x}} f(x) + \frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

## Forward-Backward Method
This is a first order method, it means that it only uses information for the first derivative of the objective function. This method exploits the differentiability of the regular part of the function and the proximal operator of the non-differentiable part.

This method defines the following iterative scheme:

$$
y_{k+1} = prox_{\tau \partial f}(y_k - \tau \nabla g(y_k)).
$$

Therefore, for our specific case we need to calculate the gradient of the differentiable part which can be obtained using stadard vectorial calculus:

$$
 \nabla \|f-\mathbb{K}^\top y\|_2^2 = \mathbb{K}(f-\mathbb{K}^\top y)
$$

Now, it is required to calculate the proximal operator of the non-differentiable part:

$$
prox_{\tau\partial f}(x) = arg\min_{\bar{x}} \delta_{B_\lambda}(y) + \frac{1}{2\tau}\|\bar{y}-y\|_2^2
$$

As can be seen this leads to a new optimization problem, this problem can be solved by analyzing the first order optimality condition of this problem:

$$
0 \in \partial \delta_{B_\lambda}(y) + (\bar{y}-y)
$$

Where $$ \partial \delta_{B_\lambda}(y) $$ is the convex subdifferendial of the indicator function of the $$\lambda$$-ball.
