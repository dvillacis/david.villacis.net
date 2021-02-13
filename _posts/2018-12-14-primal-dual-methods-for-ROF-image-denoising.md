---
layout: post
title: Solving the ROF Denoising Model using Chambolle-Pock
short_description: "We will solve a saddle point formulation of the ROF model using chambolle-pock method"
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/chambolle_pock_denoising/cp_playing_cards.png"
tags: chambolle pock backward rudin osher fatemi rof fenchel saddle point optimization matlab
active: blog
---
Let us recall the variational model used for image denoising proposed by Rudin, Osher and Fatemi (ROF)

$$ \max_x \frac{1}{2}\|x-f\|^2+\alpha\|\mathbb{K}x\|_{2,1}, $$

where $$f\in\mathbb{R}^n$$ is the noise contaminated image, $$\mathbb{K}\in\mathbb{R}^{n\times m}$$ the discrete gradient operator and $$\|\cdot\|_{2,1}$$ is the anisotropic Total Variation (TV) seminorm.

This model as presented is very challenging to solve, mainly due to the nonsmoothnes of the total variation seminorm. In [this previous blog post](https://david.villacis.net/blog/2017/02/01/rof-denoising-dual-formulation/) we used a fenchel dual formulation to solve it by using a projected gradient descent algorithm as described [here](https://david.villacis.net/blog/2017/04/06/solving-dual-rof-using-projected-gradient-descent/). Now, our goal is to improve this results by relying in a *Saddle Point Formulation*. Indeed, such formulation involves primal and dual variables, rougly speaking it removes the complexity nonsmooth term (TV) by relying on its dual formulation while leaving the smooth part untouched.

Particularly, the Rudin-Osher-Fatemi (ROF) Model can be solved using a saddle point formulation which can be written as:

$$ \min_x \max_y g(x) + \langle \mathbb{K} x,y\rangle - f^*(y) $$

$$ \min_x \max_y \frac{1}{2}\|x-f\|_2^2 + \langle \mathbb{K} x,y\rangle - \delta_{B_\alpha}(y) $$

This separation allows us to use a primal-dual numerical approach to solve this problem. On order to do so, this method relies on the calculation of the proximal operators of the previously presented functions.

## The Proximal Operator
The proximal operator is a key piece when dealing with non-differentiable functions. It consists on modifying the original function with a Moreau-Yosida regularization:

$$ \hat{f}(x) = f(x)+\frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

This regularization makes the original function to be strictly convex, which implies that it has a unique minimizer (this is not the case for convex functions). Indeed, this minimizer is known as the *proximal* of the function $$f$$:

$$ prox_{\tau\partial f}(x) = arg\min_{\bar{x}} f(x) + \frac{1}{2\tau}\|\bar{x}-x\|_2^2 $$

Eventhough this equation shows a formula for calculating this proximal operator, we have to solve an optimization problem to find it. Therefore, this formulation is useful when the solution of this problem has a closed form that is easy to implement and can be done offline.


## Chambolle-Pock Method
This is a first order method, it means that it only uses information for the first derivative of the objective function. Both proximal operators for the smooth and nonsmooth parts are easy to calculate and implement.

This method defines the following iterative scheme:

$$ x_{k+1} = prox_{\tau \partial g}(x_k-\tau\mathbb{K}^\top y_k)$$

$$ \bar{x}_{k+1} = 2x_{k+1} - x_k $$

$$ y_{k+1} = prox_{\sigma \partial f^*}(y_k+\sigma \mathbb{K}\bar{x}_{k+1}) $$

Therefore, for our specific case we need to calculate both the proximal of $$g$$ and $$f^*$$:

#### Proximal $$g(x) = \frac{1}{2}\|x-f\|_2^2$$
In this calculation we will use as described in the definition, a Moreau-Yosida regularization of $$g$$:

$$
prox_{\tau\partial g}(x) = arg\min_{\bar x} \frac{1}{2}\|\bar{x}-f\|_2^2 + \frac{1}{2\tau}\|\bar{x}-x\|_2^2
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

#### Proximal $$f^*(y) = \delta_{B_\alpha}(y)$$

$$
prox_{\sigma\partial f^*}(y) = arg\min_{\bar{y}} \delta_{B_\alpha}(\bar{y}) + \frac{1}{2\sigma}\|\bar{y}-y\|_2^2
$$

As can be seen this leads to a new optimization problem, this problem can be solved by analyzing the first order optimality condition of this problem:

$$
0 \in \partial \delta_{B_\alpha}(\bar{y}) + (\bar{y}-y)
$$

Where $$ \partial \delta_{B_\alpha}(y) $$ is the convex subdifferendial of the indicator function of the $$\lambda$$-ball. Which can be characterized by the piwelwise orthogonal projection onto an $$l_2$$ ball of radius $$\alpha$$, and can be easily computed using this formula

$$
\tilde{y}_j = \frac{y_j}{\max(1,\alpha^{-1}\|y\|)},\;\forall j=1,\dots,n
$$

### Stopping Criteria
As usual with these primal-dual methos, it is customary to use a measure for convergence. In this case it corresponds to the primal-dual gap

$$
PD(x,y) = f(\mathbb{K}x)+g(x)+f^*(y)+g^*(-\mathbb{K}^*y)
$$

where $$g^*$$ is the classical fenchel conjugate of $$g$$ which we previously reviewed.

## MATLAB Code
A working implementation can be found in our newly created repository [Bilevel Toolbox](https://github.com/dvillacis/bilevel_toolbox). The main part of the method has the following form

```matlab
function [sol,gap] = solve_rof_cp_single_gaussian(f,param)

  % Start the counter
  t1 = tic;

  [M, N] = size(f);
  f = f(:);

  nabla = gradient_matrix(M,N);

  p = zeros(M*N*2,1);
  sol = f;
  sol_ = sol;
  L = sqrt(8);
  tau = 0.01;
  sigma = 1/tau/L^2;
  gap = [];

  % Auxiliary terms
  a = 1/(1+tau);
  b = 1/param.alpha;

  for k = 1:param.maxiter

    % Dual step calculation
    p = p + sigma*nabla*sol_;
    p = reshape(p,M*N,2);
    p = reshape(bsxfun(@rdivide,p,max(1, b*rssq(p,2))), M*N*2,1);

    %Primal Step calculation
    sol_ = sol;
    sol = sol - tau*nabla'*p;
    sol = a*(sol+tau*f);

    % Overrelaxation step
    sol_ = 2*sol -sol_;

    ga = compute_rof_pd_gap(nabla, sol, p, f, param.alpha, 0, M, N);

    gap = [gap, ga];

    if mod(k, param.check) == 0 && param.verbose > 1
      fprintf('rof_cp: iter = %4d, gap = %f\n', k, ga);
    end

    if ga < param.tol
      break;
    end

  end
  sol = reshape(sol,M,N);

  % Print summary
  if param.verbose>0
    fprintf(['\n ','ROF_CHAMBOLLE_POCK',':\n']);
    fprintf(' %i iterations\n', k);
    fprintf(' Primal-Dual Gap: %f \n', gap(end));
    fprintf(' Execution Time: %f \n\n', toc(t1));
  end

end
```
As for the primal dual gap it can be implemented by calculating the primal and dual values as follows

```matlab
function [gap, primal, dual] = compute_rof_pd_gap(nabla, u, p, f, lambda, M, N)

  div_p = nabla'*p;
  nabla_u = nabla*u;

  nu = sqrt(sum(reshape(nabla_u, M*N, 2).^2,2));

  primal = lambda * sum(nu) + 1/2*norm(u-f)^2;
  dual = -norm(div_p)^2/2 + f'*div_p;

  gap = primal-dual;
end
```

## Numerical Experiment
We can see that this method converges much faster than the dual version presented in a previous post. Furthermore, as we can calculate the projection operator pixelwise, it is possible to take advantage of parallel computing, in particular GPU computing to make this operation very efficient. A recent effort in implementing this can be read in the following [paper](https://lajc.epn.edu.ec/index.php/LAJC/article/view/133).

We can see the algorithm indeed denoises the image as expected

![cp_playing_cards](https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/chambolle_pock_denoising/cp_playing_cards.png)

And we can see a decrease in the primal dual gap, confirming it is indeed converging to a solution.

![cp_gap](https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/chambolle_pock_denoising/cp_gap.png)
