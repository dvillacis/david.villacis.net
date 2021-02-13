---
layout: post
title: Solving the Dual ROF Denoising Model using Projected Gradient Descent
short_description: "We will solve the dual ROF model using a projected gradient descent algorithm using python, numpy and scipy"
img_url: "https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/gradient_descent_dual_rof.png"
tags: projected gradient rudin osher fatemi rof fenchel dual optimization python numpy scipy coitoolbox
active: blog
---
In this blog post we will continue talking about the Image Denoise problem described in my [previous post](https://david.villacis.net/blog/2017/02/01/rof-denoising-dual-formulation/). In this post we described a dual formulation for the ROF Image Denoise model.

We can see that in this particular case the problem becomes easier to tackle, given that there are no non-differentiable terms. Now we have to find ways to solve this dual problem numerically, with this goal in mind we will propose a first approach using a gradient descent algorithm. If you are unfamiliar with this optimization techniques you can refer to this [book](https://link.springer.com/book/10.1007%2Fb106451), wich provides a comprenhensive foundation.

The main idea behind gradient based methods is to define an iterative scheme for solving an optimization problem in which an objective function decreses its value when its variable changes in the direction of the gradient of such function. In the case of the ROF model we can estimate its value at each iteration

$$ROF(u) = \frac{1}{2}\| u-f \|_2^2 + \lambda \| \mathbb{K}u \|_{2,1}$$

This value can be used to give a bound on the $$l_2$$ error $$ \frac{1}{2} \|u - u^* \| $$. We will be implementing this gap using a the following python function:

```python
def ROF_value(f,x,y,clambda):
    r""" Compute the ROF cost functional

    Parameters
    ----------
    f : numpy array
        Noisy input image
    x : numpy array
        Primal variable value
    y : numpy array
        Dual variable value
    clambda : float
        Tickonov regularization parameter
    """
    a = np.linalg.norm((f-x).flatten())**2/2
    b = np.sum(np.sqrt(np.sum(y**2,axis=2)).flatten())
    return a+clambda*b
```

Recalling a version of the dual formulation presented in the previous post, we know that the dual form of the ROF model is:

$$
\min \left\{ \frac{1}{2} \|\mathbb{K}^* p\|_2^2 - \langle \mathbb{K}^* p, f \rangle : \sqrt{p_j^2 + p_{m+j}^2} \le \lambda \; \forall j = 1,\dots,m\right\}
$$

If we look closely we can define a simpler form of this optimization problem by implementing a projection operator into the unit $$\mathbb{R}^n$$ ball:

```python
def prox_project(clambda, z):
    r""" Projection to the clambda-ball

    Parameters
    ----------
    clambda : float
        Radius of the ball
    z : numpy array
        data to be projected

    """
    nrm = np.sqrt(z[:,:,0]**2 + z[:,:,1]**2)
    fact = np.minimum(clambda, nrm)
    fact = np.divide(fact,nrm, out=np.zeros_like(fact), where=nrm!=0)

    y = np.zeros(z.shape)
    y[:,:,0] = np.multiply(z[:,:,0],fact)
    y[:,:,1] = np.multiply(z[:,:,1],fact)
    return y
```

leaving the rest of the function differentiable, therefore we can explicity calculate its gradient:

$$\nabla ROF(u) = \mathbb{K}(\mathbb{K}^* p - f)$$

I will be using a python module that I'm developing called Bilevel Imaging Toolbox [(BIToolbox)](https://github.com/dvillacis/BilevelImagingToolbox), it is still in its early stages, but there you can find an implementation for a projected gradient descent algorithm. The relevant part of the code is detailed:

```python
def projectedGD_ROF(image, clambda, iters=100):
    r""" 2D Dual ROF solver using Projected Gradient Descent Method

    Parameters
    ----------
    image : numpy array
        The noisy image we are processing
    clambda : float
        The non-negative weight in the optimization problem
    iters : int
        Number of iterations allowed

    """
    print("2D Dual ROF solver using Projected Gradient Descent method")

    start_time = timeit.default_timer()
    op = operators.make_finite_differences_operator(image.shape,'fn',1)
    y = op.val(image)
    x = image
    vallog = np.zeros(iters)
    alpha = 0.1 #Line search parameter

    for i in range(iters):
        y -= alpha * op.val(op.conj(y)-image)
        y = operators.prox_project(clambda,y)
        x = image - op.conj(y)
        vallog[i] = ROF_value(image,x,op.val(x),clambda)

    print("Finished Projected Gradient Descent Dual ROF denoising in %d iterations and %f sec"%(iters,timeit.default_timer()-start_time))
    return (x,vallog)
```

In this toolbox there is also a predefined methods for obtaining the discrete gradient matrix *nabla* andd adding different kind of noise to images. If we pass this predefined functions to this algorithm, we will get the following results:

![alt text](https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/gradient_descent_dual_rof.png "ROF Denoised Lena")

And we can take a look to the ROF cost evolution:

![alt text](https://storage.googleapis.com/pagina-personal.appspot.com/img_blog/dual_rof_denoising/ROF_cost_evolution_PGD.png "ROF Cost Evolution")


Therefore, we can make use of traditional optimization techniques to solve this problem. Although, this method requires the gradient to be Lipschitz continuous, in practice this method doesn't work for some values of the step size and $$ \lambda$$. Therefore we will explore other techniques in the next posts. Finally, a couple of words on the descent algorithm itself. In this experiments we are using a fixed step size for the descent step, this type of algorithms are very sensitive to the step length choice. In consequence, if we want to improve the results presented it is important to implement a smarter line search strategy such as *Armijo* or *Wolfe*. Furthermore, various acceleration mechanisms take into account the variation of this parameter to accelerate convergence.
