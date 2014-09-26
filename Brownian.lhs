% The Existence of Brownian Motion
% Dominic Steinitz
% 24th September 2014

---
bibliography: Kalman.bib
---

In 1905, Einstein published five remarkable papers including:

Über die von der molekularkinetischen Theorie der Wärme geforderte
Bewegung von in ruhenden Flüssigkeiten suspendierten Teilchen.

which roughly translates as

On the movement of small particles suspended in a stationary liquid
demanded by the molecular-kinetic theory of heat.

which explained why the phenomenon observed by

Brownian motion is a stochastic process $\big(W_t\big)_{0 \le t \le 1}$ (presumably
$W$ for Wiener) such that

* $W_0(\omega) = 0$ for all $\omega$

* For every $0 \le s \le t \le 1$, $W_t - W_s$ is independent of $\{
W_u : u \le s\}$ and $~ {\cal{N}}(0,t-s)$.

* The map $t \mapsto W_t(\omega)$ is a continuous function of $t$ for
all $\omega$.

The
[Kolmogorov-Daniell](http://www.hss.caltech.edu/~kcb/Notes/Kolmogorov.pdf)
theorem guarantees that a stochastic process satisfying the first two
conditions exists but does not tell us that the paths are continuous. Further this theorem is not constructive relying on the axion of choice. Instead let us follow [@Ciesielski61].

> {-# LANGUAGE DataKinds           #-}
> {-# LANGUAGE TypeOperators       #-}
> {-# LANGUAGE KindSignatures      #-}
> {-# LANGUAGE TypeFamilies        #-}
> {-# LANGUAGE TypeOperators       #-}
> {-# LANGUAGE Rank2Types          #-}
> {-# LANGUAGE ScopedTypeVariables #-}

> import GHC.TypeLits
> import Data.Proxy

> import Diagrams.Backend.CmdLine
> import Diagrams.Backend.Cairo.CmdLine

> import Diagrams.Prelude hiding ( sample, render, Renderable )
> import Graphics.Rendering.Chart
> import Graphics.Rendering.Chart.Backend.Diagrams
> import Data.Default.Class
> import System.IO.Unsafe

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

> haarEtouffe :: Int -> Int -> Double -> Double
> haarEtouffe n k t
>   | n <= 0               = error "n must be >= 1"
>   | k `mod`2 == 0        = error "k must be odd"
>   | k < 0 || k > 2^n - 1 = error "k must be >=0 and <= 2^n -1"
>   | (k' - 1) * 2 ** (-n') < t && t <= k'       * 2 ** (-n') =  2 ** ((n' - 1) / 2)
>   | k'       * 2 ** (-n') < t && t <= (k' + 1) * 2 ** (-n') = -2 ** ((n' - 1) / 2)
>   | otherwise                                               =  0
>   where
>     k' = fromIntegral k
>     n' = fromIntegral n

> n :: Int
> n = 100

> xss :: [[(Double, Double)]]
> xss = map (\(m, k) -> map (\i -> let x = fromIntegral i / fromIntegral n in (x, haarEtouffe m k x)) [0..n - 1]) [(1,1), (2,1), (2,3), (3,1), (3,3), (3,5), (3,7)]

Bibliography and Resources
--------------------------
