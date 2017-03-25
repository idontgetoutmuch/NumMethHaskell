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

import qualified Prelude as P

import Numeric.Units.Dimensional.Prelude hiding (Unit)
import Numeric.Units.Dimensional

-- a_E :: Double -> Double
a_E p = (13.8 P.* p P.+ 0.04) P./ (0.08 P.+ p)

mu_F, nu_F, gamma :: Double

mu_F = 5.9

nu_F = 120.0

muPerMl :: (Fractional a, Num a) => Unit 'NonMetric DConcentration a
muPerMl = (milli mole) / (milli litre)

bigE_0 :: Concentration Double
bigE_0 = 15.0 *~ muPerMl

gamma = 0.0083

\end{code}

@Ackleh200621 gives $f$ and $g$ as

\begin{code}
fAckleh _t m = a P./ (1 P.+ k P.* (m P.** r))
  where
    a = 15600
    k = 0.0382
    r = 6.96
\end{code}

@BELAIR1995317 gives $f$ as

\begin{code}
fBelair _t m = a P./ (1 + k P.* (m P.** r))
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
deltaT, deltaMu, deltaNu :: Time Double

deltaT = 0.05 *~ day
deltaMu = 0.01 *~ day
deltaNu = 0.05 *~ day

bigT :: Time Double
bigT = 100.0 *~ day

muF, nuF :: Time Double
muF = 5.9 *~ day
nuF = 120.0 *~ day

bigK :: Int
bigK = floor (bigT / deltaT /~ one)

n1 :: Int
n1 = floor (muF / deltaT /~ one)

n2 :: Int
n2 = floor (nuF / deltaT /~ one)

ts :: [Time Double]
ts = take bigK $ 0.0 *~ day : (map (+ deltaT) ts)
\end{code}

\begin{code}
g'0 :: Dimensionless Double
g'0 = gThibodeaux bigE_0

betaAckleh :: Time Double -> Frequency Double
betaAckleh mu
  | mu < (0 *~ day) = error "betaAckleh: negative age"
  | mu < (3 *~ day) = 2.773 *~ (one / day)
  | otherwise       = 0.000 *~ (one / day)

gAckleh :: Concentration Double -> Dimensionless Double
gAckleh _e = 1.0 *~ one

sigmaAckleh :: Time Double ->
               Time Double ->
               Concentration Double ->
               Frequency Double
sigmaAckleh mu _t e = betaAckleh mu * gAckleh e

sigmaThibodeaux :: Time Double ->
                   Time Double ->
                   Concentration Double ->
                   Frequency Double
sigmaThibodeaux mu _t e
  | mu < (0 *~ day) = error "sigmaThibodeaux: negative age"
  | mu < (3 *~ day) = (2.773 *~ (one / day))
                      - (0.5 *~ (muPerMl / day)) / ((1 *~ muPerMl) + e)
  | otherwise       = (0.0 *~ (one /day))
                      - (0.5 *~ (muPerMl / day)) / ((1 *~ muPerMl) + e)

d_1'0 :: Int -> Dimensionless Double
d_1'0 i = (1 *~ one) + (g'0 * deltaT / deltaMu)
          - deltaT * sigmaThibodeaux ((fromIntegral i *~ one) * deltaMu) undefined bigE_0

s_0 :: Time Double -> Quantity (DAmountOfSubstance / DConcentration) Double
s_0 = const (4.45e7 *~ (mole / muPerMl))

lowers :: [Dimensionless Double]
lowers = replicate n1 (g'0 * deltaT / deltaMu)

diags :: [Dimensionless Double]
diags = g'0 : map d_1'0 [1..n1]

uppers :: [Dimensionless Double]
uppers = replicate n1 (0.0 *~ one)
\end{code}

As in @Thibodeaux2011 we give quantities in terms of cells per
kilogram of body weight.

\begin{code}
p_0 :: Time Double -> Quantity (DAmountOfSubstance / DTime / DMass) Double
p_0 mu = (1e11 *~ one) * pAux mu
  where
    pAux mu
      | mu < (0 *~ day) = error "P_0: negative age"
      | mu < (3 *~ day) = 8.55e-6 *~ (mole / day / kilo gram) *
                          exp ((2.7519 *~ (one / day)) * mu)
      | otherwise       = 8.55e-6 *~ (mole / day / kilo gram) *
                          exp (8.319 *~ one - (0.0211 *~ (one /day)) * mu)
\end{code}

\begin{code}
main :: IO ()
main = undefined
\end{code}

References
==========
