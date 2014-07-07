% Fun with Filters
% Dominic Steinitz
% 3rd July 2014

Noisy Observation
=================

Suppose we wish to estimate the mean of a sample drawn from a normal
distribution. In the Bayesian approach, we know the prior distribution
for the mean (it could be a non-informative prior) and then we update
this with our observations to create the posterior, the latter giving
us improved information about the distribution of the mean. In symbols

$$
p(\theta \,\vert\, x) \propto p(x \,\vert\, \theta)p(\theta)
$$

Typically, the samples are chosen to be independent, and all of the
data is used to perform the update but, given independence, there is
no particular reason to do that, updates can performed one at a time
and the result is the same; nor is the order of update
important. Being a bit imprecise, we have

$$
p(z \,\vert\, x, y) = p(z, x, y)p(x, y) = p(z, x, y)p(x)p(y) =
p((z \,\vert\, x) \,\vert\, y) =
p((z \,\vert\, y) \,\vert\, x)
$$

The standard notation in Bayesian statistics is to denote the
parameters of interest as $\theta \in \mathbb{R}^p$ and the
observations as $x \in \mathbb{R}^n$. For reasons that will become
apparent, let us change notation and label the parameters as $x$ and
the observations as $y$.


Let us take a very simple example of a prior $X \sim {\cal{N}}(0,
\sigma^2)$ where $\sigma^2$ is known and then sample from a normal
distribution with mean $x$ and variance for the $i$-th sample $c_i^2$
where $c_i$ is known (normally we would not know the variance but
adding this generality would only clutter the exposition
unnecessarily).

$$
p(y_i \,\vert\, x) = \frac{1}{\sqrt{2\pi c_i^2}}\exp\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg)
$$

The likelihood is then

$$
p(\boldsymbol{y} \,\vert\, x) = \prod_{i=1}^n \frac{1}{\sqrt{2\pi c_i^2}}\exp\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg)
$$

As we have already noted, instead of using this with the prior to
calculate the posterior, we can update the prior with each observation
separately. Suppose that we have obtained the posterior given $i - 1$
samples (we do not know this is normally distributed yet but we soon
will):

$$
p(x \,\vert\, y_1,\ldots,y_{i-1}) = {\cal{N}}(\hat{x}_{i-1}, \hat{\sigma}^2_{i-1})
$$

Then we have

$$
\begin{aligned}
p(x \,\vert\, y_1,\ldots,y_{i}) &\propto p(y_i \,\vert\, x)p(x \,\vert\, y_1,\ldots,y_{i-1}) \\
&\propto \exp-\bigg(\frac{(y_i - x)^2}{2c_i^2}\bigg) \exp-\bigg(\frac{(x - \hat{x}_{i-1})^2}{2\hat{\sigma}_{i-1}^2}\bigg) \\
&\propto \exp-\Bigg(\frac{x^2}{c_i^2} - \frac{2xy_i}{c_i^2} + \frac{x^2}{\hat{\sigma}_{i-1}^2} - \frac{2x\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}\Bigg) \\
&\propto \exp-\Bigg( x^2\Bigg(\frac{1}{c_i^2} + \frac{1}{\hat{\sigma}_{i-1}^2}\Bigg) - 2x\Bigg(\frac{y_i}{c_i^2} + \frac{\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}\Bigg)\Bigg)
\end{aligned}
$$

Writing

$$
\frac{1}{\hat{\sigma}_{i}^2} \triangleq \frac{1}{c_i^2} + \frac{1}{\hat{\sigma}_{i-1}^2}
$$

and then completing the square we also obtain

$$
\frac{\hat{x}_{i}}{\hat{\sigma}_{i}^2} \triangleq \frac{y_i}{c_i^2} + \frac{\hat{x}_{i-1}}{\hat{\sigma}_{i-1}^2}
$$

More Formally
-------------

Now let's be a bit more formal about conditional probability and use
the notation of $\sigma$-algebras to define ${\cal{F}}_i =
\sigma\{Y_1,\ldots, Y_i\}$ and $M_i \triangleq \mathbb{E}(X \,\vert\,
{\cal{F}}_i)$ where $Y_i = X + \epsilon_i$, $X$ is as before and
$\epsilon_i \sim {\cal{N}}(0, c_k^2)$. We have previously calculated
that $M_i = \hat{x}_i$ and that ${\cal{E}}((X - M_i)^2 \,\vert\, Y_1,
\ldots Y_i) = \hat{\sigma}_{i}^2$ and the tower law for conditional
probabilities then allows us to conclude ${\cal{E}}((X - M_i)^2) =
\hat{\sigma}_{i}^2$. By [Jensen's
inequality](http://en.wikipedia.org/wiki/Jensen%27s_inequality), we have

$$
{\cal{E}}(M_i^2) = {\cal{E}}({\cal{E}}(X \,\vert\, {\cal{F}}_i)^2)) \leq
{\cal{E}}({\cal{E}}(X^2 \,\vert\, {\cal{F}}_i))) =
{\cal{E}}(X^2) = \sigma^2
$$

Hence $M$ is bounded in $L^2$ and therefore converges in $L^2$ and
almost surely to $M_\infty \triangleq {\cal{E}}(X \,\vert\,
{\cal{F}}_\infty)$.

Importance Sampling
===================

Suppose we have an random variable $X$ with pdf $1/2\exp{-\lvert
x\rvert}$ and we wish to find its second moment numerically. However,
the [random-fu](https://hackage.haskell.org/package/random-fu) package
does not support sampling from such as distribution. We notice that

$$
\int_{-\infty}^\infty x^2 \frac{1}{2} \exp{-\lvert x\rvert} \mathrm{d}x =
\int_{-\infty}^\infty x^2 \frac{\frac{1}{2} \exp{-\lvert x\rvert}}
                               {\frac{1}{\sqrt{8\pi}}{\exp{-x^2/8}}}
                      \frac{1}{\sqrt{8\pi}}{\exp{-x^2/8}}
\,\mathrm{d}x
$$

So we can sample from ${\cal{N}}(0, 4)$ and evaluate

$$
x^2 \frac{\frac{1}{2} \exp{-\lvert x\rvert}}
         {\frac{1}{\sqrt{8\pi}}{\exp{-x^2/8}}}
$$

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> import Control.Monad
> import Data.Random.Source.PureMT
> import Data.Random
> import Control.Monad.State
> import qualified Control.Monad.Writer as W

> import Formatting

> sampleImportance :: RVarT (W.Writer [Double]) ()
> sampleImportance = do
>   x <- rvarT $ Normal 0.0 2.0
>   let x2 = x^2
>       u = x2 * 0.5 * exp (-(abs x))
>       v = (exp ((-x2)/8)) * (recip (sqrt (8*pi)))
>       w = u / v
>   lift $ W.tell [w]
>   return ()

> runImportance :: Int -> [Double]
> runImportance n =
>   snd $
>   W.runWriter $
>   evalStateT (sample (replicateM n sampleImportance))
>              (pureMT 2)

We can run this 10,000 times to get an estimate.

    [ghci]
    format (fixed 2) (sum (runImportance 10000) / 10000)

Since we know that the $n$-th moment of the exponential distribution
is $n! / \lambda^n$ where $\lambda$ is the rate (1 in this example),
the exact answer is 2 which is not too far from our estimate using
importance sampling.

The value of

$$
w(x) = \frac{1}{N}\frac{\frac{1}{2} \exp{-\lvert x\rvert}}
                       {\frac{1}{\sqrt{8\pi}}{\exp{-x^2/8}}}
     = \frac{p(x)}{\pi(x)}
$$

is called the weight, $p$ is the pdf from which we wish to sample and
$\pi$ is the pdf of the importance distribution.

Importance Sampling Approximation of the Posterior
==================================================

Suppose that the posterior distribution of a model in which we are
interested has a complicated functional form and that we therefore
wish to approximate it in some way. First assume that we wish to
calculate the expectation of some arbitrary function $f$ of the
parameters.

$$
{\cal{E}}(f(\boldsymbol{\theta}) \,\vert\, x_1, \ldots x_T) =
\int_\Omega f(\boldsymbol{\theta}) p(\boldsymbol{\theta} \, \vert \, x_1, \ldots x_T) \,\mathrm{d}\boldsymbol{\theta}
$$

As before we can re-write this

$$
\int_\Omega f(\boldsymbol{\theta}) p(\boldsymbol{\theta} \, \vert \, x_1, \ldots x_T) \,\mathrm{d}\boldsymbol{\theta} =
\int_\Omega f(\boldsymbol{\theta}) \frac{p(\boldsymbol{\theta} \, \vert \, x_1, \ldots x_T)}
                                 {\pi(\boldsymbol{\theta} \, \vert \, x_1, \ldots x_T)}
     \pi(\boldsymbol{\theta} \, \vert \, x_1, \ldots x_T) \,\mathrm{d}\boldsymbol{\theta}
$$

We can now sample $\Theta \sim \pi(\boldsymbol{\theta} \, \vert \,
x_1, \ldots x_T)$ repeatedly to obtain

$$
{\cal{E}}(f(\boldsymbol{\theta}) \,\vert\, x_1, \ldots x_T) \approx \frac{1}{N}\sum_1^N
f(\boldsymbol{\Theta_i}) \frac{p(\boldsymbol{\Theta_i} \, \vert \, x_1, \ldots x_T)}
                            {\pi(\boldsymbol{\Theta_i} \, \vert \, x_1, \ldots x_T)} =
\sum_1^N w_if(\boldsymbol{\Theta_i})
$$

where the weights $w_i$ are defined as before by

$$
w_i = \frac{1}{N} \frac{p(\boldsymbol{\Theta_i} \, \vert \, x_1, \ldots x_T)}
                       {\pi(\boldsymbol{\Theta_i} \, \vert \, x_1, \ldots x_T)}
$$
