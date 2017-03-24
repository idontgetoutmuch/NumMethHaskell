% Trouble with Tribbles
% Dominic Steinitz
% 14th March 2017

---
bibliography: DynSys.bib
---

Introduction
============

McKendrick / von Foerster
-------------------------

[McKendrick](https://en.wikipedia.org/wiki/Anderson_Gray_McKendrick)
and [von Foerster](https://en.wikipedia.org/wiki/Heinz_von_Foerster)
independently derived a model of age-dependent population
growth.

Let $n(a,t)$ be the density of females of age $a$ at time $t$. The
number of females between ages $a$ and $a + \delta a$ are thus $n(a,
t)\delta a$. Assuming individuals are born at age $0$, we have

$$
\frac{\partial}{\partial t}(n(a, t)\delta a) =
J(a, t) - J(a + \delta a, t) - \mu(a, t)n(a, t)\delta a
$$

where $\mu(a, t)$ is the death rate density and $J(a, t)$ denotes the
rate of entry to the cohort of age $a$. Dividing by $\delta a$ we obtain

$$
\frac{\partial}{\partial t}n(a, t) =
 - \frac{J(a + \delta a, t) - J(a, t)}{\delta a} - \mu(a, t)n(a, t)
$$

which in the limit becomes

$$
\frac{\partial}{\partial t}n(a, t) =
 - \frac{\partial J(a, t)}{\partial a} - \mu(a, t)n(a, t)
$$

We can further assume that the rate of entry to a cohort is
proportional to the density of individuals times a velocity of aging
$v(a, t)$.

$$
J(a, t) = n(a, t)v(a, t)
$$

Occasionally there is some reason to assume that aging one year is
different to experiencing one year but we further assume $v = 1$.

We thus obtain

$$
\frac{\partial}{\partial t}n(a, t) + \frac{\partial n(a, t)}{\partial a} =
- \mu(a, t)n(a, t)
$$

Erythropoiesis
==============

The production of red blood cells which contain haemoglobin (aka
 hemoglobin),
 [Erythropoiesis](https://en.wikipedia.org/wiki/Erythropoiesis), is
 modulated in a feedback loop by
 [erythropoietin](https://en.wikipedia.org/wiki/Erythropoietin). As
 can be seen in the overview by @Torbett2009, the full feedback loop
 is complex. So as not to lose ourselves in the details and following
 @Thibodeaux2011 and @BELAIR1995317, we consider a model with two
 compartments.

* Precursors: prototype erythrocytes developing in the bone marrow
  with $p(\mu, t)$ being the density of such cells of age $\mu$ at
  time $t$.

* Erythrocytes: mature red blood cells circulating in the blood with
  $m(\mu, t)$ being the density of such cells of age $\nu$ at time
  $t$.

$$
\begin{aligned}
\frac{\partial p(\mu, t)}{\partial t} + g(E(t))\frac{\partial p(\mu, t)}{\partial \mu} &=
\sigma(\mu, t, E(t))p(\mu, t) & 0 < \mu < \mu_F & & 0 < t < T \\
\frac{\partial m(\nu, t)}{\partial t} + \phantom{g(E(t))}\frac{\partial m(\nu, t)}{\partial \nu} &=
-\gamma(\nu, t, M(t))m(\nu, t) & 0 < \nu < \nu_F & &  0 < t < T
\end{aligned}
$$

where $\sigma(\mu, t, E(t))$ is the birth rate of precursors and
$\gamma(\nu, t, M(t))$ is the death rate of erythrocytes, $g(E(t))$ is
the maturation rate of precursors and where

$$
M(t) = \int_0^{\nu_F} p(\nu, t) \,\mathrm{d}\nu
$$

As boundary conditions, we have that the number of precursors maturing
must equal the production of number of erythrocytes

$$
m(0, t) = g(E(t))p(\mu_F, t)
$$

and the production of the of the number of precursors depends on the
level of erythropoietin

$$
g(E(t))p(0, t) = \phi(t)E(t)
$$

where $\phi(t)$ is some proportionality function.

As initial conditions, we have

$$
\begin{aligned}
p(\mu, 0) &= p_0(\mu) \\
m(\nu, 0) &= m_0(\nu)
\end{aligned}
$$

We can further model the erythropoietin dynamics as

$$
\frac{\mathrm{d}E(t)}{\mathrm{d}t} = f(M(t), t) - a_E(P(t))E(t)
$$

where $f$ is the feedback function from the kidneys and the decay
rate, $a_E$ depends on the total precursor population $P(t)$
(@sawyer1987binding) although this often is taken to be a constant and
I would feel more comfortable with a more recent citation and where

$$
P(t) = \int_0^{\mu_F} p(\mu, t) \,\mathrm{d}\mu
$$

As initial condition we have

$$
E(0) = E_0
$$

A Finite Difference Attempt
===========================

Let us try solving the above model using a finite difference scheme
observing that we currently have no basis for whether it has a
solution and whether the finite difference scheme approximates such a
solution!

Divide up the age and time ranges, $[0, \mu_F]$, $[0, \nu_F]$ and $[0, T]$
into equal sub-intervals,
$[\mu_i, \mu_{i+1}]$, $[\nu_j, \nu_{j+1}]$ and $[t_k, t_{k+1}]$
where

$$
\begin{aligned}
\mu_i &= i\Delta\mu & & \mathrm{for} & i = 1 \ldots n_1 \\
\nu_j &= j\Delta\nu & & \mathrm{for} & j = 1 \ldots n_2 \\
t_k   &= k\Delta t  & & \mathrm{for} & k = 1 \ldots K
\end{aligned}
$$

where $\Delta\mu = \mu_F / n_1$, $\Delta\nu = \nu_F / n_2$ and $\Delta t = T / K$.

Denoting $p(\mu_i, t_k) = p_i^k$ and similarly we obtain

$$
\begin{aligned}
\frac{p_i^{k+1} - p_i^k}{\Delta t} + g^k\frac{p_i^{k+1} - p_{i-1}^{k+1}}{\Delta\mu} &= \sigma_i^k p_i^{k+1} \\
\frac{m_j^{k+1} - m_j^k}{\Delta t} + \phantom{g^k}\frac{m_j^{k+1} - m_{j-1}^{k+1}}{\Delta\mu} &= -\gamma_j^k m_j^{k+1}
\end{aligned}
$$

and

$$
\begin{aligned}
\frac{E^{k+1} - E^k}{\Delta t} &= f^k - a_E^k E^{k+1} \\
g^k p_0^{k+1} &= \phi^k E^k \\
m_0^{k+1}     &= g^k m_{n_1}^{k+1}
\end{aligned}
$$

Re-arranging we get

$$
\begin{aligned}
-g^k\frac{\Delta t}{\Delta \mu}p_{i-1}^{k+1} +
\bigg(1 + g^k\frac{\Delta t}{\Delta \mu} - \Delta t \sigma_i^k\bigg)p_i^{k+1} &=
p_i^k \\
\frac{\Delta t}{\Delta \mu}m_{j-1}^{k+1} +
\bigg(1 + \frac{\Delta t}{\Delta \mu} + \Delta t \gamma_j^k\bigg)m_j^{k+1} &=
m_j^k
\end{aligned}
$$

Writing

$$
d _{1,i}^k = 1 + g^k\frac{\Delta t}{\Delta \mu} - \Delta t \sigma_i^k
$$

$$
\begin{bmatrix}
g^k & 0 & 0 & \ldots & 0 & 0 \\
-g^k\frac{\Delta t}{\Delta \mu} & d_{1,1}^k & 0 & \ldots & 0 & 0\\
0 & -g^k\frac{\Delta t}{\Delta \mu} & d_{1,2}^k & \ldots & 0 & 0 \\
\ldots & \ldots & \ldots & \ldots & \ldots & \ldots \\
0 & 0 & 0 & \ldots &\ -g^k\frac{\Delta t}{\Delta \mu} & d_{1,n_1}^k \\
\end{bmatrix}
\begin{bmatrix}
p_0^{k+1} \\
p_1^{k+1} \\
p_2^{k+1} \\
\ldots \\
p_{n_1}^{k+1}
\end{bmatrix}
=
\begin{bmatrix}
\phi^k E^k \\
p_1^k \\
p_2^k \\
\ldots \\
p_{n_1}^k \\
\end{bmatrix}
$$

\begin{code}
{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE TypeOperators #-}

module Main where

-- import Data.Metrology
-- import Data.Units.SI
-- import Data.Units.SI.Prefixes
-- import Data.Metrology.Poly

import qualified Prelude

import Numeric.Units.Dimensional.Prelude hiding (Unit)
import Numeric.Units.Dimensional.NonSI (mile)
import Numeric.Units.Dimensional
import Numeric.Units.Dimensional.UnitNames (UnitName, atom)

-- data EpoInternationUnit = EpoInternationUnit

-- instance Unit EpoInternationUnit where
--   type BaseUnit EpoInternationUnit = Gram
--   conversionRatio _ = error "Epo should never be converted from International Units"

-- a_E :: Double -> Double
a_E p = (13.8 Prelude.* p Prelude.+ 0.04) Prelude./ (0.08 Prelude.+ p)

mu_F, nu_F, gamma :: Double

mu_F = 5.9

nu_F = 120.0

muPerMl :: (Fractional a, Num a) => Unit NonMetric DConcentration a
muPerMl = (milli mole) / (milli litre)

bigE_0 :: Concentration Double
bigE_0 = 15.0 *~ muPerMl

gamma = 0.0083

\end{code}

@Ackleh200621 gives $f$ and $g$ as

\begin{code}
fAckleh t m = a Prelude./ (1 Prelude.+ k Prelude.* (m Prelude.** r))
  where
    a = 15600
    k = 0.0382
    r = 6.96

gAckleh e = 1.0
\end{code}

@BELAIR1995317 gives $f$ as

\begin{code}
fBelair t m = a Prelude./ (1 + k Prelude.* (m Prelude.** r))
  where
    a = 6570
    k = 0.0382
    r = 6.96
\end{code}

@Thibodeaux2011 give $g$ as

\begin{code}
gThibodeaux :: Concentration Double  -> Dimensionless Double
gThibodeaux e = d / n
  where
    n = ((3.02 *~ one) * e + (0.31 *~ muPerMl))
    d = (30.61 *~ muPerMl) + e
\end{code}

From@Ackleh200621 but note that @Thibodeaux2011 seem to have $T = 20$.

\begin{code}
deltaT, deltaMu, deltaNu :: Double

deltaT = 0.05
deltaMu = 0.01
deltaNu = 0.05

bigT :: Double
bigT = 100.0

-- ts :: [Double]
ts = [0.0,deltaT..bigT Prelude./ deltaT]
\end{code}

\begin{code}
g'0 = gThibodeaux bigE_0

betaAckleh :: Time Double -> Frequency Double
betaAckleh mu
  | mu < (0 *~ day) = error "betaAckleh: negative age"
  | mu < (3 *~ day) = 2.773 *~ (one / day)
  | otherwise       = 0.000 *~ (one / day)

-- sigmaThibodeaux mu t e
--   | mu < 0    = error "negative mu"
--   | mu <= 3   = 2.773 - 0.5 / (1 + e)
--   | otherwise = - 0.5 / (1 + e)

-- d_1'0 :: Prelude.Int -> Prelude.Double
-- d_1'0 i = 1 Prelude.+ g'0 Prelude.* deltaT Prelude./ deltaMu Prelude.- deltaT Prelude.* undefined

-- bigAd_1'1 = gThibodeaux bigE_0 : []
\end{code}

\begin{code}
main :: IO ()
main = undefined
\end{code}

References
==========
