% Fun with Filters
% Dominic Steinitz
% 3rd July 2014

---
bibliography: Kalman.bib
---

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}

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

> {-# LANGUAGE FlexibleContexts             #-}

> import Control.Monad
> import Data.Random.Source.PureMT
> import Data.Random
> import Data.Random.Distribution.Binomial
> import Data.Random.Distribution.Multinomial
> import Control.Monad.State
> import qualified Control.Monad.Writer as W
> import Control.Monad.Loops

> import Formatting


> import Debug.Trace

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
{\cal{E}}(f({x}) \,\vert\, y_1, \ldots y_T) =
\int_\Omega f({x}) p({x} \, \vert \, y_1, \ldots y_T) \,\mathrm{d}{x}
$$

As before we can re-write this

$$
\int_\Omega f({x}) p({x} \, \vert \, y_1, \ldots y_T) \,\mathrm{d}{x} =
\int_\Omega f({x}) \frac{p({x} \, \vert \, y_1, \ldots y_T)}
                                 {\pi({x} \, \vert \, y_1, \ldots y_T)}
     \pi({x} \, \vert \, y_1, \ldots y_T) \,\mathrm{d}{x}
$$

We can now sample $X^{(i)} \sim \pi({x} \, \vert \,
y_1, \ldots y_T)$ repeatedly to obtain

$$
{\cal{E}}(f({x}) \,\vert\, y_1, \ldots y_T) \approx \frac{1}{N}\sum_1^N
f({X^{(i)}}) \frac{p({X^{(i)}} \, \vert \, y_1, \ldots y_T)}
                            {\pi({X^{(i)}} \, \vert \, y_1, \ldots y_T)} =
\sum_1^N w_if({X^{(i)}})
$$

where the weights $w_i$ are defined as before by

$$
w_i = \frac{1}{N} \frac{p({X^{(i)}} \, \vert \, y_1, \ldots y_T)}
                       {\pi({X^{(i)}} \, \vert \, y_1, \ldots y_T)}
$$

We follow [Alex
Cook](http://blog.nus.edu.sg/alexcook/teaching/sph6004/) and use the
example from [@citeulike:5986027].

> xv, nv :: Int
> xv = 51
> nv = 8197

> sampleUniform :: RVarT (W.Writer [Double]) ()
> sampleUniform = do
>   x <- rvarT StdUniform
>   lift $ W.tell [x]
>   return ()

> sampleSize = 10000

> pv = runSampler sampleUniform 2 sampleSize

> weightsRaw = map (\p -> pdf (Binomial nv p) xv) pv
> weightsSum = sum weightsRaw
> weights = map (/ weightsSum) weightsRaw

> meanPv = sum $ zipWith (*) pv weights

But

    [ghci]
    length $ filter (>= 1e-6) weights

so we may not be getting a very good estimate.

> sampleNormal :: RVarT (W.Writer [Double]) ()
> sampleNormal = do
>   let xvd   = fromIntegral xv
>       nvd   = fromIntegral nv
>       mu    = xvd / nvd
>       sigma = (sqrt (xvd / nvd) * (1 - xvd / nvd)) / nvd
>   x <- rvarT $ Normal mu sigma
>   lift $ W.tell [x]
>   return ()

> runSampler :: RVarT (W.Writer [Double]) () -> Int -> Int -> [Double]
> runSampler sampler seed n =
>   snd $
>   W.runWriter $
>   evalStateT (sample (replicateM n sampler))
>              (pureMT (fromIntegral seed))

> pvC = runSampler sampleNormal 2 sampleSize

> weightsRawC = map (\p -> pdf (Binomial nv p) xv) pvC
> weightsSumC = sum weightsRawC
> weightsC = map (/ weightsSumC) weightsRawC

> meanPvC = sum $ zipWith (*) pvC weightsC

    [ghci]
    length $ filter (>= 1e-6) weightsC

Sequential Importance Sampling
==============================

Now let us generalize our Bayesian model and allow
${X}_i$ to depend on ${X}_{i-1}$ and
${Y}_i$ depend on ${X}_i$, that is

$$
\begin{aligned}
{X}_i &\sim \condprob{p}{{x}_i}{{x}_{i-1}} \\
{Y}_i &\sim \condprob{p}{{y}_i}{{x}_{i}} \\
\end{aligned}
$$

Recall in our original model we had

$$
\begin{aligned}
{X}_i &\sim p({{x}}) \\
{Y}_i &\sim \condprob{p}{{y}_i}{{x}} \\
\end{aligned}
$$

And, as before, by Bayes and conditional independence we also had

$$
\begin{aligned}
\condprob{p}{{x}}{{y}_1, \ldots, {y}_T} & \propto
p({x})\condprob{p}{{y}_1, \ldots, {y}_T}{{x}} \\
&  =p({x})\prod_{i=1}^T\condprob{p}{{y}_i}{{x}} \\
& = p({x})\condprob{p}{{y}_T}{{x}}\prod_{i=1}^{T-1}\condprob{p}{{y}_i}{{x}} \\
& \propto \condprob{p}{{y}_T}{{x}}\condprob{p}{{x}}{{y}_1, \ldots, {y}_{T-1}}
\end{aligned}
$$

In our new and generalized model we have a whole sequence of
parameters for which we wish to find the posterior distribution.

$$
\condprob{p}{{x}_0,\ldots,{x}_T}{{y}_1,\ldots,{y}_T}
$$

Applying Bayes and the restrictions of our model that the parameters are Markovian

$$
\begin{aligned}
\condprob{p}{{x}_0,\ldots,{x}_T}{{y}_1,\ldots,{y}_T}
& \propto
\condprob{p}{{y}_T}{{x}_0,\ldots,{x}_T, {y}_1,\ldots,{y}_{T-1}}
\condprob{p}{{x}_0,\ldots,{x}_T}{{y}_1,\ldots,{y}_{T-1}} \\
& =
\condprob{p}{{y}_T}{{x}_T}
\condprob{p}{{x}_T}{{x}_0,\ldots,{x}_{T-1}, {y}_1,\ldots,{y}_{T-1}}\condprob{p}{{x}_0,\ldots,{x}_{T-1}}{{y}_1,\ldots,{y}_{T-1}} \\
& =
\condprob{p}{{y}_T}{{x}_T}
\condprob{p}{{x}_T}{{x}_{T-1}}\condprob{p}{{x}_0,\ldots,{x}_{T-1}}{{y}_1,\ldots,{y}_{T-1}} \\
\end{aligned}
$$

Now let us take an importance distribution
$\condprob{\pi}{x_0,\ldots,x_T}{y_1,\ldots,y_T}$ and as before compute
the un-normalized weights by substituting in the posterior we have
just calculated.

$$
w_T^{(i)} = \frac
{\condprob{p}{{y}_T}{{X_T^{(i)}}}
 \condprob{p}{{X_T^{(i)}}}{{X_{T-1}^{(i)}}}\condprob{p}{{X_0^{(i)}},\ldots,{X_{T-1}^{(i)}}}{{y}_1,\ldots,{y}_{T-1}}
}
{\condprob{\pi}{X^{(i)}_0,\ldots,X^{(i)}_T}{y_1,\ldots,y_T}}
$$

The problem now is how to select an appropriate importance distribution.
Suppose we can do this recursively.

$$
{\condprob{\pi}{x_0,\ldots,x_T}{y_1,\ldots,y_T}} =
{\condprob{\pi}{x_T}{x_0,\ldots,x_{T-1}, y_1,\ldots,y_T}}
{\condprob{\pi}{x_0,\ldots,x_{T-1}}{y_1,\ldots,y_{T-1}}}
$$

then

$$
\begin{aligned}
w_T^{(i)} & = \frac
{\condprob{p}{{y}_T}{{X_T^{(i)}}}
 \condprob{p}{{X_T^{(i)}}}{{X_{T-1}^{(i)}}}
}
{
{\condprob{\pi}{X_T^{(i)}}{X_0^{(i)},\ldots,X_{T-1}^{(i)}, y_1,\ldots,y_T}}
}
\frac
\condprob{p}{{X_0^{(i)}},\ldots,{X_{T-1}^{(i)}}}{{y}_1,\ldots,{y}_{T-1}}
{\condprob{\pi}{X_0^{(i)},\ldots,X_{T-1}^{(i)}}{y_1,\ldots,y_{T-1}}} \\
& =
\frac
{\condprob{p}{{y}_T}{{X_T^{(i)}}}
 \condprob{p}{{X_T^{(i)}}}{{X_{T-1}^{(i)}}}
}
{
{\condprob{\pi}{X_T^{(i)}}{X_0^{(i)},\ldots,X_{T-1}^{(i)}, y_1,\ldots,y_T}}
}w^{(i)}_{T-1}
\end{aligned}
$$

We can start the whole sampling process by sampling $X_0^{(i)} \sim
p(x_0)$ from the prior and setting $w_0^{(i)} = 1/N$ where $N$ is the
total number of samples.

Further, if we choose ${\condprob{\pi}{X_T}{X_0,\ldots,X_{T-1},
y_1,\ldots,y_T}} = {\condprob{\pi}{X_T}{X_{T-1}, y_1,\ldots,y_T}}$ then we
do not need to keep the whole history of parameters
$X_0^{(i)},\ldots,X_{T-1}^{(i)}$ only the previous parameters.

This simplifies the recursive weight equations

$$
\begin{aligned}
w_T^{(i)}
& =
\frac
{\condprob{p}{{y}_T}{{X_T^{(i)}}}
 \condprob{p}{{X_T^{(i)}}}{{X_{T-1}^{(i)}}}
}
{
{\condprob{\pi}{X_T^{(i)}}{X_{T-1}^{(i)}, y_1,\ldots,y_T}}
}w^{(i)}_{T-1}
\end{aligned}
$$

If we further take $\condprob{\pi}{X_T^{(i)}}{X_{T-1}^{(i)},
y_1,\ldots,y_T}$ to be $\condprob{p}{{X_T^{(i)}}}{{X_{T-1}^{(i)}}}$
weighted by $w_{T-1}^{(i)}$ then

$$
\begin{aligned}
w_T^{(i)}
& =
\condprob{p}{{y}_T}{{\tilde{X}_T^{(i)}}}
\end{aligned}
$$

where

$$
\tilde{X}_T \sim \sum_{i=1}^N w_{T-1}^{(i)} \condprob{p}{X_T}{X_{T-1}^{(i)}}
$$

This second sampling according to the weights is known as **resampling**.

> muPrior = 0.0
> sigmaPrior = 1.0
> muLikelihood = 0.0
> cs = repeat 1.0
>
> bigN = 400

> normalPdf :: Double -> Double -> Double -> Double
> normalPdf mu sigma x =
>   (recip (sqrt (2 * pi * sigma2))) * (exp ((-(x - mu)^2) / (2 * sigma2)))
>   where
>     sigma2 = sigma^2

> sir :: [(Double, Double)] -> Double -> Double ->
>        RVarT (W.Writer [([Double], [Int])]) [(Double, Double)]
> sir weightsMusPrev y sigma = do
>   let n           = length weightsMusPrev
>       weightsPrev = map fst weightsMusPrev
>
>   nParticless <- rvarT $ Multinomial weightsPrev n
>
>   let musPrev  = map snd weightsMusPrev
>   let _musTilde = concatMap (\i -> replicate (nParticless!!i) (musPrev!!i)) [0..n - 1]
>
>   musNew <- return musPrev
>
>   let weightsNew = map (\i -> normalPdf (musNew!!i) sigma y) [0..n - 1]
>
>   lift $ W.tell [(musNew, nParticless)]
>   return (zip weightsNew musNew)

> initSir :: Int -> RVarT (W.Writer [([Double], [Int])]) [(Double, Double)]
> initSir n = do
>   mus <- replicateM n (rvarT $ Normal muPrior sigmaPrior)
>   lift $ W.tell [(mus, replicate n 1)]
>   return (zip (replicate n (recip (fromIntegral n))) mus)

> createObs :: Int -> RVar [Double]
> createObs n = do
>   x <- rvarT (Normal muPrior sigmaPrior)
>   trace (show x) $ return ()
>   ys <- mapM (\c -> rvarT (Normal x c)) (take n cs)
>   return ys

> obss = evalState (sample (createObs 10)) (pureMT 2)

> testPF :: RVarT (W.Writer [([Double], [Int])]) [(Double, Double)]
> testPF = do
>   wms <- initSir bigN
>   let sirs = zipWith (\obs c -> (\wm -> sir wm obs c)) obss cs
>   foldr (>=>) return sirs wms

> runPF :: [([Double], [Int])]
> runPF = snd (W.runWriter (evalStateT (sample testPF) (pureMT 2)))

Bibliography
============