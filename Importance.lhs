% Importance Sampling
% Dominic Steinitz
% 14th August 2014

---
bibliography: Kalman.bib
---

\newcommand{\condprob} [3] {#1\left(#2 \,\vert\, #3\right)}

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

> module Importance where

> import Control.Monad
> import Data.Random.Source.PureMT
> import Data.Random
> import Data.Random.Distribution.Binomial
> import Data.Random.Distribution.Beta
> import Control.Monad.State
> import qualified Control.Monad.Writer as W


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
    import Formatting
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
{\mathbb{E}}(f({x}) \,\vert\, y_1, \ldots y_T) =
\int_\Omega f({x}) p({x} \, \vert \, y_1, \ldots y_T) \,\mathrm{d}{x}
$$

Using Bayes

$$
\int_\Omega f({x}) \condprob{p}{x}{y_1, \ldots y_T} \,\mathrm{d}{x} =
\frac{1}{Z}\int_\Omega f({x}) \condprob{p}{y_1, \ldots y_T}{x}p(x) \,\mathrm{d}{x}
$$

where $Z$ is some normalizing constant.

As before we can re-write this using a proposal distribution $\pi(x)$

$$
\frac{1}{Z}\int_\Omega f({x}) \condprob{p}{y_1, \ldots y_T}{x}p(x) \,\mathrm{d}{x} =
\frac{1}{Z}\int_\Omega \frac{f({x}) \condprob{p}{y_1, \ldots y_T}{x}p(x)}{\pi(x)}\pi(x) \,\mathrm{d}{x}
$$

We can now sample $X^{(i)} \sim \pi({x})$ repeatedly to obtain

$$
{\mathbb{E}}(f({x}) \,\vert\, y_1, \ldots y_T) \approx \frac{1}{ZN}\sum_1^N
f({X^{(i)}}) \frac{p(y_1, \ldots y_T \, \vert \, {X^{(i)}})p({X^{(i)}})}
                            {\pi({X^{(i)}})} =
\sum_1^N w_if({X^{(i)}})
$$

where the weights $w_i$ are defined as before by

$$
w_i = \frac{1}{ZN} \frac{p(y_1, \ldots y_T \, \vert \, {X^{(i)}})p({X^{(i)}})}
                        {\pi({X^{(i)}})}
$$

We follow [Alex
Cook](http://blog.nus.edu.sg/alexcook/teaching/sph6004/) and use the
example from [@citeulike:5986027]. We take the prior as $\sim
{\cal{Be}}(1,1)$ and use ${\cal{U}}(0.0,1.0)$ as the proposal
distribution. In this case the proposal and the prior are identical
just expressed differently and therefore cancel.

Note that we use the log of the pdf in our calculations otherwise we
suffer from (silent) underflow, e.g.,

    [ghci]
    pdf (Binomial nv (0.4 :: Double)) xv

On the other hand if we use the log pdf form

    [ghci]
    logPdf (Binomial nv (0.4 :: Double)) xv

> xv, nv :: Int
> xv = 51
> nv = 8197

> sampleUniform :: RVarT (W.Writer [Double]) ()
> sampleUniform = do
>   x <- rvarT StdUniform
>   lift $ W.tell [x]
>   return ()

> runSampler :: RVarT (W.Writer [Double]) () ->
>               Int -> Int -> [Double]
> runSampler sampler seed n =
>   snd $
>   W.runWriter $
>   evalStateT (sample (replicateM n sampler))
>              (pureMT (fromIntegral seed))

> sampleSize :: Int
> sampleSize = 1000

> pv :: [Double]
> pv = runSampler sampleUniform 2 sampleSize

> logWeightsRaw :: [Double]
> logWeightsRaw = map (\p -> logPdf (Beta 1.0 1.0) p +
>                            logPdf (Binomial nv p) xv -
>                            logPdf StdUniform p) pv

> logWeightsMax :: Double
> logWeightsMax = maximum logWeightsRaw
>
> weightsRaw :: [Double]
> weightsRaw = map (\w -> exp (w - logWeightsMax)) logWeightsRaw

> weightsSum :: Double
> weightsSum = sum weightsRaw

> weights :: [Double]
> weights = map (/ weightsSum) weightsRaw

> meanPv :: Double
> meanPv = sum $ zipWith (*) pv weights
>
> meanPv2 :: Double
> meanPv2 = sum $ zipWith (\p w -> p * p * w) pv weights
>
> varPv :: Double
> varPv = meanPv2 - meanPv * meanPv

We get the answer

    [ghci]
    meanPv

But if we look at the size of the weights and the effective sample size

    [ghci]
    length $ filter (>= 1e-6) weights
    (sum weights)^2 / (sum $ map (^2) weights)

so we may not be getting a very good estimate. Let's try

> sampleNormal :: RVarT (W.Writer [Double]) ()
> sampleNormal = do
>   x <- rvarT $ Normal meanPv (sqrt varPv)
>   lift $ W.tell [x]
>   return ()

> pvC :: [Double]
> pvC = runSampler sampleNormal 3 sampleSize

> logWeightsRawC :: [Double]
> logWeightsRawC = map (\p -> logPdf (Beta 1.0 1.0) p +
>                             logPdf (Binomial nv p) xv -
>                             logPdf (Normal meanPv (sqrt varPv)) p) pvC

> logWeightsMaxC :: Double
> logWeightsMaxC = maximum logWeightsRawC
>
> weightsRawC :: [Double]
> weightsRawC = map (\w -> exp (w - logWeightsMaxC)) logWeightsRawC

> weightsSumC :: Double
> weightsSumC = sum weightsRawC

> weightsC :: [Double]
> weightsC = map (/ weightsSumC) weightsRawC

> meanPvC :: Double
> meanPvC = sum $ zipWith (*) pvC weightsC

> meanPvC2 :: Double
> meanPvC2 = sum $ zipWith (\p w -> p * p * w) pvC weightsC
>
> varPvC :: Double
> varPvC = meanPvC2 - meanPvC * meanPvC

Now the weights and the effective size are more re-assuring

    [ghci]
    length $ filter (>= 1e-6) weightsC
    (sum weightsC)^2 / (sum $ map (^2) weightsC)

And we can take more confidence in the estimate

    [ghci]
    meanPvC

Bibliography
============