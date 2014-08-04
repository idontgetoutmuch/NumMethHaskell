% Fun with (Kalman) Filters Part I
% Dominic Steinitz
% 21st July 2014

---
bibliography: Kalman.bib
---

\newcommand{\condprob} [3] {#1 \left( #2 \,\vert\, #3 \right)}
\newcommand{\normal} [2] {{\cal{N}}\big( #1 , #2\big)}

Suppose we have particle moving in 1 dimension, where the initial
velocity is sampled from a distribution. We can observe the position
of the particle at fixed intervals and we wish to estimate its initial
velocity. For generality, let us assume that the positions and the
velocities can be perturbed at each interval and that our measurements
are noisy.

$$
\begin{aligned}
x_i &= x_{i-1} + \Delta T v_{i-1} + \psi^{(x)}_i \\
v_i &= v_{i-1} + \psi^{(v)}_i \\
y_i &= a_i x_i + \upsilon_i
\end{aligned}
$$

where $\psi^{(x)}_i, \psi^{(v)}_i$ and $\upsilon_i$ are all IID
normal with means of 0 and variances of $\sigma^2_x,
\sigma^2_v$ and $\sigma^2_y$

We can re-write this as

$$
\begin{aligned}
\boldsymbol{x}_i &= \boldsymbol{A}_{i-1}\boldsymbol{x}_{i-1} + \boldsymbol{\psi}_{i-1} \\
\boldsymbol{y}_i &= \boldsymbol{H}_i\boldsymbol{x}_i + \boldsymbol{\upsilon}_i
\end{aligned}
$$

where

$$
\boldsymbol{A}_i =
  \begin{bmatrix}
    1 & \Delta T\\
    0 & 1\\
  \end{bmatrix}
,\quad
\boldsymbol{H}_i =
  \begin{bmatrix}
    a_i \\
  \end{bmatrix}
,\quad
\boldsymbol{\psi}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(x)}_i\big)
,\quad
\boldsymbol{\Sigma}^{(x)}_i =
  \begin{bmatrix}
    \sigma^2_{x} & 0\\
    0 & \sigma^2_{v} \\
  \end{bmatrix}
,\quad
\boldsymbol{\upsilon}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(y)}_i\big)
,\quad
\boldsymbol{\Sigma}^{(y)}_i =
  \begin{bmatrix}
    \sigma^2_{z} \\
  \end{bmatrix}
$$

Let us denote the mean and variance of
$\boldsymbol{X}_i\,\vert\,\boldsymbol{Y}_{i-1}$ as
$\hat{\boldsymbol{x}}^\flat_i$ and $\hat{\boldsymbol{\Sigma}}^\flat_i$
respectively and note that

$$
\begin{aligned}
{\boldsymbol{Y}_i}\,\vert\,{\boldsymbol{Y}_{i-1}} =
{\boldsymbol{H}_i\boldsymbol{X}_i\,\vert\,{\boldsymbol{Y}_{i-1}} + \boldsymbol{\Upsilon}_i}\,\vert\,{\boldsymbol{Y}_{i-1}} =
{\boldsymbol{H}_i\boldsymbol{X}_i\,\vert\,{\boldsymbol{Y}_{i-1}} + \boldsymbol{\Upsilon}_i}
\end{aligned}
$$

Since ${\boldsymbol{X}_i}\,\vert\,{\boldsymbol{Y}_{i-1}}$ and
${\boldsymbol{Y}_i}\,\vert\,{\boldsymbol{Y}_{i-1}}$ are jointly
Gaussian and recalling that ${\hat{\boldsymbol{\Sigma}}^\flat_i}^\top =
\hat{\boldsymbol{\Sigma}}^\flat_i$ as covariance matrices are
symmetric, we can calculate their mean and covariance matrix as

$$
\begin{bmatrix}
    \hat{\boldsymbol{x}}^\flat_i \\
    \boldsymbol{H}_i\hat{\boldsymbol{x}}^\flat_i
\end{bmatrix}
,\quad
\begin{bmatrix}
    \hat{\boldsymbol{\Sigma}}^\flat_i & \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top \\
     \boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i & \boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i \\
\end{bmatrix}
$$

We can now use [standard formul&aelig;](http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Conditional_distributions) which says if

$$
\begin{bmatrix}
    \boldsymbol{X} \\
    \boldsymbol{Y}
\end{bmatrix}
\sim
{\cal{N}}
\begin{bmatrix}
\begin{bmatrix}
    \boldsymbol{\mu}_x \\
    \boldsymbol{\mu}_y
\end{bmatrix}
&
,
&
\begin{bmatrix}
    \boldsymbol{\Sigma}_x & \boldsymbol{\Sigma}_{xy} \\
    \boldsymbol{\Sigma}^\top_{xy} & \boldsymbol{\Sigma}_y
\end{bmatrix}
\end{bmatrix}
$$

then

$$
\boldsymbol{X}\,\vert\,\boldsymbol{Y}=\boldsymbol{y} \sim \normal{\boldsymbol{\mu}_x + \boldsymbol{\Sigma}_{xy}\boldsymbol{\Sigma}^{-1}_y(\boldsymbol{y} - \boldsymbol{\mu}_y)}{\boldsymbol{\Sigma}_x - \boldsymbol{\Sigma}_{xy}\boldsymbol{\Sigma}^{-1}_y\boldsymbol{\Sigma}^\top_{xy}}
$$

and apply this to

$$
(\boldsymbol{X}_i\,\vert\, \boldsymbol{Y}_{i-1})\,\vert\,(\boldsymbol{Y}_i\,\vert\, \boldsymbol{Y}_{i-1})
$$

to give

$$
\boldsymbol{X}_i\,\vert\, \boldsymbol{Y}_{i} = \boldsymbol{y}_i
\sim
\normal{\hat{\boldsymbol{x}}^\flat_i + \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top
\big(\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i\big)^{-1}
(\boldsymbol{y}_i - \boldsymbol{H}_i\hat{\boldsymbol{x}}^\flat_i)}
{\hat{\boldsymbol{\Sigma}}^\flat_i - \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top(\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i)^{-1}\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i}
$$

This is called the measurement update; more explicitly

$$
\begin{aligned}
\hat{\boldsymbol{x}}^i &\triangleq
\hat{\boldsymbol{x}}^\flat_i +
\hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top
\big(\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i\big)^{-1}
(\boldsymbol{y}_i - \boldsymbol{H}_i\hat{\boldsymbol{x}}^\flat_i) \\
\hat{\boldsymbol{\Sigma}}_i &\triangleq
{\hat{\boldsymbol{\Sigma}}^\flat_i - \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top(\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i \boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i)^{-1}\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i}
\end{aligned}
$$

Sometimes the measurement residual $\boldsymbol{v}_i$, the measurement
prediction covariance $\boldsymbol{S}_i$ and the filter gain
$\boldsymbol{K}_i$ are defined and the measurement update is written
as

$$
\begin{aligned}
\boldsymbol{v}_i & \triangleq
\boldsymbol{y}_i - \boldsymbol{H}_i\hat{\boldsymbol{x}}^\flat_i \\
\boldsymbol{S}_i & \triangleq
\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i \\
\boldsymbol{K}_i & \triangleq \hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top\boldsymbol{S}^{-1}_i \\
\hat{\boldsymbol{x}}^i &\triangleq \hat{\boldsymbol{x}}^\flat_i + \boldsymbol{K}_i\boldsymbol{v}_i \\
\hat{\boldsymbol{\Sigma}}_i &\triangleq \hat{\boldsymbol{\Sigma}}^\flat_i - \boldsymbol{K}_i\boldsymbol{S}_i\boldsymbol{K}^\top_i
\end{aligned}
$$

We further have that

$$
\begin{aligned}
{\boldsymbol{X}_i}\,\vert\,{\boldsymbol{Y}_{i-1}} =
{\boldsymbol{A}_i\boldsymbol{X}_{i-1}\,\vert\,{\boldsymbol{Y}_{i-1}} + \boldsymbol{\Psi}_{i-1}}\,\vert\,{\boldsymbol{Y}_{i-1}} =
{\boldsymbol{A}_i\boldsymbol{X}_{i-1}\,\vert\,{\boldsymbol{Y}_{i-1}} + \boldsymbol{\Psi}_i}
\end{aligned}
$$

We thus obtain the Kalman filter prediction step:

$$
\begin{aligned}
\hat{\boldsymbol{x}}^\flat_i &=
\boldsymbol{A}_{i-1}\hat{\boldsymbol{x}}_{i-1} \\
\hat{\boldsymbol{\Sigma}}^\flat_i &= \boldsymbol{A}_{i-1}
                                     \hat{\boldsymbol{\Sigma}}_{i-1}
                                     \boldsymbol{A}_{i-1}^\top
                                   + \boldsymbol{\Sigma}^{(x)}_{i-1}
\end{aligned}
$$

Further information can be found in [@Boyd:EE363:Online], [@kleeman1996understanding] and [@sarkka2013bayesian].

A Haskell Implementation
========================

Let us start with a prior normal distribution with a mean position and
velocity of 0 with large variances and no correlation.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE DataKinds                    #-}
> {-# LANGUAGE ScopedTypeVariables          #-}
> {-# LANGUAGE RankNTypes                   #-}

> import GHC.TypeLits
> import Numeric.LinearAlgebra.Static

> import Data.Maybe ( fromJust )

> import Data.Random.Source.PureMT
> import Data.Random hiding ( gamma )
> import Data.Random.Distribution.Multinomial
> import Control.Monad.State
> import qualified Control.Monad.Writer as W

> muPrior :: R 2
> muPrior = vector [0.0, 0.0]

> sigmaPrior :: Sq 2
> sigmaPrior = matrix [ 1000.0,    0.0
>                     ,    0.0, 1000.0
>                     ]

> deltaT :: Double
> deltaT = 0.1

> bigA :: Sq 2
> bigA = matrix [ 1, deltaT
>               , 0,      1
>               ]

> a :: Double
> a = 2.0

> bigH :: L 1 2
> bigH = matrix [ a, 0
>               ]

> bigSigmaY :: Sq 1
> bigSigmaY = matrix [ 0.1 ]

> bigSigmaX :: Sq 2
> bigSigmaX = matrix [ 0.1, 0.0
>                    , 0.0, 0.1
>                    ]

> inv :: KnownNat n => Sq n -> Maybe (Sq n)
> inv m = linSolve m eye

> outer ::  forall m n . (KnownNat m, KnownNat n) =>
>           R n -> Sq n -> L m n -> Sq m -> Sq n -> Sq n -> [R m] -> [(R n, Sq n)]
> outer muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = result
>   where
>     result = scanl update (muPrior, sigmaPrior) ys
>
>     update :: (R n, Sq n) -> R m -> (R n, Sq n)
>     update (xHatFlat, bigSigmaHatFlat) y = (xHatFlatNew, bigSigmaHatFlatNew)
>       where
>         v = y - bigH #> xHatFlat
>         bigS = bigH <> bigSigmaHatFlat <> (tr bigH) + bigSigmaY
>         bigK = bigSigmaHatFlat <> (tr bigH) <> (fromJust (inv bigS))
>         xHat = xHatFlat + bigK #> v
>         bigSigmaHat = bigSigmaHatFlat - bigK <> bigS <> (tr bigK)
>         xHatFlatNew = bigA #> xHat
>         bigSigmaHatFlatNew = bigA <> bigSigmaHat <> (tr bigA) + bigSigmaX


> createObs :: RVar Double
> createObs = do
>   x <- rvarT (Normal 0.0 0.1)
>   return x

> foo seed nSamples =
>   evalState (replicateM nSamples (sample (Normal 0.0 0.1)))
>   (pureMT $ fromIntegral seed)

Bibliography
============