% The Existence of Brownian Motion
% Dominic Steinitz
% 24th September 2014

---
bibliography: Kalman.bib
---

In 1905, Einstein published five remarkable papers including, "Über
die von der molekularkinetischen Theorie der Wärme geforderte Bewegung
von in ruhenden Flüssigkeiten suspendierten Teilchen." which roughly
translates as "On the movement of small particles suspended in a
stationary liquid demanded by the molecular-kinetic theory of heat."
giving the first explanation of the phenomenon observed by Robert
Brown in 1827 of small particles of pollen moving at random when
suspended in water.

The eponymously named Brownian motion is a stochastic process
$\big(W_t\big)_{0 \le t \le 1}$ (presumably $W$ for Wiener) such that

1. $W_0(\omega) = 0$ for all $\omega$

2. For all $0 \le t_1 \le t \le t_2 \le, \ldots, t_n \le 1$, $W_{t_1},
W_{t_2} - W_{t_3}, \ldots W_{t_{n-1}} - W_{t_n}$ are independent.

3. $W_{t+h} - W_t \sim {\cal{N}}(0,h)$.

4. The map $t \mapsto W_t(\omega)$ is a continuous function of $t$ for
all $\omega$.

The
[Kolmogorov-Daniell](http://www.hss.caltech.edu/~kcb/Notes/Kolmogorov.pdf)
theorem guarantees that a stochastic process satisfying the first two
conditions exists but does not tell us that the paths are
continuous.

For example, suppose that we have constructed Brownian motion as
defined by the above conditions then take a random variable $U \sim
{\cal{U}}[0,1]$. Define a new process

$$
 \tilde{W}_t =
  \begin{cases}
   W & \text{if } t \neq U \\
   0 & \text{if } t = U
  \end{cases}
$$

This has the same finite dimensional distributions as $W_t$ as
$(W_{t_1}, W_{t_2}, W_{t_3}, \ldots W_{t_{n-1}}, W_{t_n})$ and
$(\tilde{W}_{t_1}, \tilde{W}_{t_2}, \tilde{W}_{t_3}, \ldots
\tilde{W}_{t_{n-1}}, \tilde{W}_{t_n})$ are equal unless $U \in \{t_1,
\ldots, t_n\}$ and this set has measure 0. This process satisifes
conditions 1--3 but is not continuous

$$
\lim_{t \uparrow U} \tilde{W}_t = \lim_{t \uparrow U} W_t = W_U
$$

and

$$
\mathbb{P}(W_U = 0) = \int_0^1
$$

Further this theorem is not constructive relying on the
axion of choice. Instead let us follow [@Ciesielski61].

It is tempting to think of Brownian motion as the limit in some sense
of a random walk.

$$
x_t^{(N)} = \sum_{i=1}^\floor{Nt}\frac{\xi_i}{N}
$$

where $0 \le t \le 1$ However, the processes, $x_t^{(1)}, x_t^{(2)},
\ldots$ are discontinuous and it is not clear how one would prove that
the limit of discontinuous processes is in fact continuous.

Let $\{\phi_i\}$ be a complete orthonormal basis for $L^2[0,1]$. That
is any $f$ for which $\int_0^1 f^2 \mathrm{d}\mu$ exists then

$$
f = \sum_{i=1}^\infty \langle f, \phi_i\rangle
$$

where $\mu$ is Lesbegue measure and $\langle\ldots,\ldots\rangle$ is
the inner product for $L^2$ defined as usual by

$$
\langle f, g\rangle = \int_0^1 fg\mathrm{d}\mu$
$$

We know such bases exist, for example, the [Fourier
expansion](http://en.wikipedia.org/wiki/Fourier_series) and [Legendre
polynomials](http://en.wikipedia.org/wiki/Legendre_polynomials).

We defined the so-called Haar wavelets for $n = 0, 1, \ldots$ and $k = 1, 3, 5, \ldots, 2^n - 1$.

$$
H_{n,k}(t) = +2^{(n-1)/2)} (k - 1)2^{-n} < t \le k2^{-n}
           = -2^{(n-1)/2)} k2^{-n}       < t \le (k + 1)2^{-n}
$$

Using Haskell's capabilities for dependently typed programming we can
express these as shown below.

> {-# OPTIONS_GHC -Wall                     #-}
> {-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
> {-# OPTIONS_GHC -fno-warn-type-defaults   #-}
> {-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
> {-# OPTIONS_GHC -fno-warn-missing-methods #-}
> {-# OPTIONS_GHC -fno-warn-orphans         #-}

> {-# LANGUAGE DataKinds           #-}
> {-# LANGUAGE TypeOperators       #-}
> {-# LANGUAGE KindSignatures      #-}
> {-# LANGUAGE TypeFamilies        #-}
> {-# LANGUAGE TypeOperators       #-}
> {-# LANGUAGE Rank2Types          #-}
> {-# LANGUAGE ScopedTypeVariables #-}
> {-# LANGUAGE PolyKinds           #-}

> module Brownian where

> import GHC.TypeLits
> import Data.Proxy

> import Numeric.Integration.TanhSinh

> data Haar (a :: Nat) (b :: Nat) = Haar { unHaar :: Double -> Double }

> haar :: forall n k .
>         (KnownNat n, KnownNat k, (2 * k - 1 <=? 2^n - 1) ~ 'True) =>
>         Haar n (2 * k - 1)
> haar = Haar g where
>   g t | (k' - 1) * 2 ** (-n') < t && t <= k'       * 2 ** (-n') =  2 ** ((n' - 1) / 2)
>       | k'       * 2 ** (-n') < t && t <= (k' + 1) * 2 ** (-n') = -2 ** ((n' - 1) / 2)
>       | otherwise                                               =  0
>     where
>         n' = fromIntegral (natVal (Proxy :: Proxy n))
>         k' = 2 * (fromIntegral (natVal (Proxy :: Proxy k))) - 1

Now for example we can evaluate

> haar11 :: Double -> Double
> haar11 = unHaar (haar :: Haar 1 1)

    [ghci]
    haar11 0.75

but we if we try to evaluate *haar :: Haar 1 2* we get a type error.

> type family ZipWith (f :: a -> b -> c) (as :: [a]) (bs :: [b]) :: [c] where
>   ZipWith f (a ': as) (b ': bs) = (f a b) ': (ZipWith f as bs)
>   ZipWith f as        bs        = '[]

Rather than go too far into dependently typed programming which would
distract us from the existence proof, let us re-express this function
in a more traditional way (which will only give us errors at runtime).

> haarEtouffee :: Int -> Int -> Double -> Double
> haarEtouffee n k t
>   | n <= 0               = error "n must be >= 1"
>   | k `mod`2 == 0        = error "k must be odd"
>   | k < 0 || k > 2^n - 1 = error "k must be >=0 and <= 2^n -1"
>   | (k' - 1) * 2 ** (-n') < t && t <= k'       * 2 ** (-n') =  2 ** ((n' - 1) / 2)
>   | k'       * 2 ** (-n') < t && t <= (k' + 1) * 2 ** (-n') = -2 ** ((n' - 1) / 2)
>   | otherwise                                               =  0
>   where
>     k' = fromIntegral k
>     n' = fromIntegral n

Here are the first few Haar wavelets.

```{.dia height='300'}
import Brownian
import BrownianChart

dia = diag 2 "n = 1, k = 1" (xss!!0)
````

```{.dia height='300'}
import Brownian
import BrownianChart

dia = (diag 2 "n = 2, k = 1" (xss!!1) |||
       diag 2 "n = 2, k = 3" (xss!!2))
````

```{.dia height='300'}
import Brownian
import BrownianChart

dia = ((diag 2 "n = 3, k = 1" (xss!!3) |||
        diag 2 "n = 3, k = 3" (xss!!4) |||
        diag 2 "n = 3, k = 5" (xss!!5) |||
        diag 2 "n = 3, k = 7" (xss!!6)))
````

> schauderEtouffee :: Int -> Int -> Double -> Double
> schauderEtouffee n k t = result (absolute 1e-6 (parSimpson (haarEtouffee n k) 0 t))

> n :: Int
> n = 1000

> xss :: [[(Double, Double)]]
> xss = map (\(m, k) -> map (\i -> let x = fromIntegral i / fromIntegral n in (x, haarEtouffee m k x)) [0..n - 1]) [(1,1), (2,1), (2,3), (3,1), (3,3), (3,5), (3,7)]

> yss :: [[(Double, Double)]]
> yss = map (\(m, k) -> map (\i -> let x = fromIntegral i / fromIntegral n in (x, schauderEtouffee m k x)) [0..n - 1]) [(1,1), (2,1), (2,3), (3,1), (3,3), (3,5), (3,7)]

```{.dia height='300'}
import Brownian
import BrownianChart

dia = diag 0.5 "Foo" (yss!!0)
````

```{.dia height='300'}
import Brownian
import BrownianChart

dia = (diag 0.5 "Foo" (yss!!1) |||
       diag 0.5 "Baz" (yss!!2))
````

```{.dia height='300'}
import Brownian
import BrownianChart

dia = ((diag 0.5 "Foo" (yss!!3) |||
        diag 0.5 "Baz" (yss!!4) |||
        diag 0.5 "Baz" (yss!!5) |||
        diag 0.5 "Baz" (yss!!6)))
````

Bibliography and Resources
--------------------------
