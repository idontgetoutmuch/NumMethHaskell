% A Dynamical System Example
% Dominic Steinitz
% 12th December 2015

---
bibliography: Kalman.bib
---

Introduction
============

Apparently and under special conditions, the concentrations of the
metabolites involved in yeast glycolysis can oscillate:
@Bier2000. They give the following dynamical system.

$$
\begin{aligned}
\frac{\mathrm{d}[G]}{\mathrm{d}t} &= V_{\mathrm{in}} - k_1 [G] [ATP] \\
\frac{\mathrm{d}[ATP]}{\mathrm{d}t} &= \pi
\end{aligned}
$$

where $V_{\mathrm{in}}$ represents the rate at which glucose enters
the system and $k_1$ represents the rate at which glucose is converted
into [ATP](https://en.wikipedia.org/wiki/Adenosine_triphosphate).

> module BierEtAl where

> import Numeric.GSL.ODE
> import Numeric.LinearAlgebra

> xdot t [x,v] = [v, -0.95*x - 0.1*v]

> ts = linspace 100 (0,20 :: Double)

> sol = odeSolve xdot [10,0] ts

main = mplot (ts : toColumns sol)

Bier, Bakker, & Westerhoff published a very simple one
(Biophys. J. 78:1087-1093, 2000)

> vIn' = 0.10
> k1' = 0.02
> kp' = 6
> bigKm' = 12.0
> atpInit = 4.0
> bigGInit = 3.0

> deltaT = 0.1
> totTime = 1000.0
> bigN = floor $ totTime / deltaT

> data Bier = Bier { vin :: Double
>                  , k1 :: Double
>                  , kp :: Double
>                  , km :: Double
>                  }
>   deriving Show

> bier vIn k1 kp bigKm = ydot
>   where
>     ydot t [atp, g] = [ 2 * k1 * g * atp - (kp * atp) / (atp + bigKm)
>                       , vIn - k1 * g * atp]


> us = linspace bigN (0, 500 :: Double)
> us' = toList us

> tol :: Bier -> [[Double]]
> tol b =
>   map toList $
>   toColumns $
>   odeSolve (bier (vin b) (k1 b) (kp b) (km b)) [atpInit, bigGInit] us

> bier11 vin = tol (Bier {vin = vin, k1 = k1', kp = kp', km = bigKm'})
> 

*Main> let mins = map (map maximum . map (drop 9000) . bier11) [0.1,0.2..1.6]
*Main> let mins = map (map minimum . map (drop 9000) . bier11) [0.1,0.2..1.6]
*Main> let maxs = map (map maximum . map (drop 9000) . bier11) [0.1,0.2..1.6]
*Main> zipWith (zipWith (-)) maxs mins
[[10.203632755679338,4.317902636079657],
[17.41769751835152,20.356676523447213],
[0.3447276080418492,14.50788840287894],
[16.544243323228525,17.543667352329248],
[15.737773849466983,15.272180092281396],
[14.480447001826633,14.1415645039917],
[12.600308792865956,11.916503279523486],
[9.749314488605487,8.985392272576505],
[5.300678173001263,4.7422363608151565],
[1.06848710408807,0.8628621763411157],
[7.199939417051748e-2,5.381937464709985e-2],
[1.6850502574889958e-3,1.3224805287066488e-3],
[1.9259713305075365e-5,1.2806362158279683e-5],
[9.038267378969067e-8,6.782357075962864e-8],
[1.7469403701397823e-10,6.895017889974042e-11],
[0.0,0.0]]