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

> {-# OPTIONS_GHC -Wall #-}
> {-# LANGUAGE TypeFamilies #-}
> {-# LANGUAGE NoImplicitPrelude #-}
> {-# LANGUAGE FlexibleContexts #-}
> {-# LANGUAGE DataKinds #-}
> {-# LANGUAGE TypeOperators #-}

> module Tribbles where

> import qualified Prelude as P

> import Numeric.Units.Dimensional.Prelude hiding (Unit)
> import Numeric.Units.Dimensional

> import Numeric.LinearAlgebra
> import Numeric.Integration.TanhSinh

> import Control.Monad.Writer
> import Control.Monad.Loops


> muPerMl :: (Fractional a, Num a) => Unit 'NonMetric DConcentration a
> muPerMl = (milli mole) / (milli litre)

> bigE'0 :: Concentration Double
> bigE'0 = 15.0 *~ muPerMl

@Thibodeaux2011 give $g$ as

> gThibodeaux :: Concentration Double  -> Dimensionless Double
> gThibodeaux e = d / n
>   where
>     n = ((3.02 *~ one) * e + (0.31 *~ muPerMl))
>     d = (30.61 *~ muPerMl) + e

From@Ackleh200621 but note that @Thibodeaux2011 seem to have $T = 20$.

> deltaT, deltaMu, deltaNu :: Time Double

> deltaT = 0.05 *~ day
> deltaMu = 0.01 *~ day
> deltaNu = 0.05 *~ day

> bigT :: Time Double
> bigT = 100.0 *~ day

> muF, nuF :: Time Double
> muF = 5.9 *~ day
> nuF = 120.0 *~ day

> bigK :: Int
> bigK = floor (bigT / deltaT /~ one)

> n1 :: Int
> n1 = floor (muF / deltaMu /~ one)

> n2 :: Int
> n2 = floor (nuF / deltaNu /~ one)

> ts :: [Time Double]
> ts = take bigK $ 0.0 *~ day : (map (+ deltaT) ts)

> g'0 :: Dimensionless Double
> g'0 = gThibodeaux bigE'0

> betaAckleh :: Time Double -> Frequency Double
> betaAckleh mu
>   | mu < (0 *~ day) = error "betaAckleh: negative age"
>   | mu < (3 *~ day) = 2.773 *~ (one / day)
>   | otherwise       = 0.000 *~ (one / day)

> gAckleh :: Concentration Double -> Dimensionless Double
> gAckleh _e = 1.0 *~ one

> sigmaAckleh :: Time Double ->
>                Time Double ->
>                Concentration Double ->
>                Frequency Double
> sigmaAckleh mu _t e = betaAckleh mu * gAckleh e

> betaThibodeaux :: Time Double ->
>                   Frequency Double
> betaThibodeaux mu
>   | mu < (0 *~ day) = error "betaThibodeaux: negative age"
>   | mu < (3 *~ day) = (2.773 *~ (one / day))
>   | otherwise       = (0.0 *~ (one /day))

> alphaThibodeaux :: Concentration Double ->
>                    Frequency Double
> alphaThibodeaux e = (0.5 *~ (muPerMl / day)) / ((1 *~ muPerMl) + e)

> sigmaThibodeaux :: Time Double ->
>                    Time Double ->
>                    Concentration Double ->
>                    Frequency Double
> sigmaThibodeaux mu _t e = gThibodeaux e * (betaThibodeaux mu - alphaThibodeaux e)

> d_1'0 :: Int -> Dimensionless Double
> d_1'0 i = (1 *~ one) + (g'0 * deltaT / deltaMu)
>           - deltaT * sigmaThibodeaux ((fromIntegral i *~ one) * deltaMu) undefined bigE'0

> lowers :: [Dimensionless Double]
> lowers = replicate n1 (negate $ g'0 * deltaT / deltaMu)

> diags :: [Dimensionless Double]
> diags = g'0 : map d_1'0 [1..n1]

> uppers :: [Dimensionless Double]
> uppers = replicate n1 (0.0 *~ one)

As in @Thibodeaux2011 we give quantities in terms of cells per
kilogram of body weight.

> type CellDensity = Quantity (DAmountOfSubstance / DTime / DMass)

> p'0 :: Time Double -> CellDensity Double
> p'0 mu' = (1e11 *~ one) * pAux mu'
>   where
>     pAux mu
>       | mu < (0 *~ day) = error "p'0: negative age"
>       | mu < (3 *~ day) = 8.55e-6 *~ (mole / day / kilo gram) *
>                           exp ((2.7519 *~ (one / day)) * mu)
>       | otherwise       = 8.55e-6 *~ (mole / day / kilo gram) *
>                           exp (8.319 *~ one - (0.0211 *~ (one / day)) * mu)

> m_0 :: Time Double -> CellDensity Double
> m_0 nu' = (1e11 *~ one) * mAux nu'
>   where
>     mAux nu
>       | nu < (0 *~ day) = error "m_0: age less than zero"
>       | otherwise       = 0.039827  *~ (mole / day / kilo gram) *
>                           exp (((-0.0083) *~ (one / day)) * nu)

Let's check that these give plausible results.

> m_0Untyped :: Double -> Double
> m_0Untyped nu = m_0 (nu *~ day) /~ (mole / day / kilo gram)

> p'0Untyped :: Double -> Double
> p'0Untyped mu = p'0 (mu *~ day) /~ (mole / day / kilo gram)


    [ghci]
    import Numeric.Integration.TanhSinh
    result $ relative 1e-6 $ parTrap m_0Untyped 0.001 (nuF /~ day)
    result $ relative 1e-6 $ parTrap p'0Untyped 0.001 (muF /~ day)

@Thibodeaux2011 does not give a definition for $\phi$ so we use the
equivalent $s_0$ from @Ackleh200621 which references @Banks2003:
"$\times 10^{11}$ erythrocytes/kg body weight $\times$ mL plasma/mU Epo/day"

> s_0 :: Time Double ->
>        Quantity (DAmountOfSubstance / DTime / DMass / DConcentration) Double
> s_0 = const ((1e11 *~ one) * (4.45e-7 *~ (mole / day / kilo gram / muPerMl)))

> b'0 :: [CellDensity Double]
> b'0 = (s_0 (0.0 *~ day) * bigE'0) : (take n1 $ map p'0 (iterate (+ deltaMu) deltaMu))

> p'1 :: Matrix Double
> p'1 = triDiagSolve (fromList (map (/~ one) lowers))
>                    (fromList (map (/~ one) diags))
>                    (fromList (map (/~ one) uppers))
>                    (((n1 P.+1 )><1) (map (/~ (mole / second / kilo gram)) b'0))

Death rates of mature erythrocytes.

> gammaThibodeaux :: Time Double ->
>                    Time Double ->
>                    Quantity (DAmountOfSubstance / DMass) Double ->
>                    Frequency Double
> gammaThibodeaux _nu _t _bigM = 0.0083 *~ (one / day)

> gammaAckleh :: Time Double ->
>                Time Double ->
>                Quantity (DAmountOfSubstance / DMass) Double ->
>                Frequency Double
> gammaAckleh _nu _t bigM = (0.01 *~ (kilo gram / mole / day)) * bigM + 0.0001 *~ (one / day)

For the intial mature erythrocyte population we can either use the
integral of the initial distribution

> bigM'0 :: Quantity (DAmountOfSubstance / DMass) Double
> bigM'0 = r *~ (mole / kilo gram)
>  where
>    r = result $ relative 1e-6 $ parTrap m_0Untyped 0.001 (nuF /~ day)

    [ghci]
    bigM'0

or we can use the sum of the values used in the finite difference approximation

> bigM'0' :: Quantity (DAmountOfSubstance / DMass) Double
> bigM'0' = (* deltaNu) $ sum $ map m_0 $ take n2 $ iterate (+ deltaNu) (0.0 *~ day)

    [ghci]
    bigM'0'

> d_2'0 :: Int -> Dimensionless Double
> d_2'0 i = (1 *~ one) + (g'0 * deltaT / deltaNu)
>           + deltaT * gammaThibodeaux ((fromIntegral i *~ one) * deltaNu) undefined bigM'0

> lowers2 :: [Dimensionless Double]
> lowers2 = replicate n2 (negate $ deltaT / deltaNu)

> diags2 :: [Dimensionless Double]
> diags2 = (1.0 *~ one) : map d_2'0 [1..n2]

> uppers2 :: [Dimensionless Double]
> uppers2 = replicate n2 (0.0 *~ one)

> b_2'0 :: [CellDensity Double]
> b_2'0 = (g'0 * ((p'1 `atIndex` (n1,0)) *~ (mole / second / kilo gram))) :
>         (take n2 $ map m_0 (iterate (+ deltaNu) deltaNu))

> m'1 :: Matrix Double
> m'1 = triDiagSolve (fromList (map (/~ one) lowers2))
>                    (fromList (map (/~ one) diags2))
>                    (fromList (map (/~ one) uppers2))
>                    (((n2 P.+ 1)><1) (map (/~ (mole / second / kilo gram)) b_2'0))

@Ackleh2013 and @Ackleh200621 give $f$

> fAckleh :: Time Double ->
>            Quantity (DAmountOfSubstance / DMass) Double ->
>            Quantity (DConcentration / DTime) Double
> fAckleh _t bigM = a / ((1.0 *~ one) + k * (bigM' ** r))
>   where
>     a = 15600 *~ (muPerMl / day)
>     k = 0.0382 *~ one
>     r = 6.96 *~ one
>     bigM' = ((bigM /~ (mole / kilo gram)) *~ one) * (1e-11 *~ one)

The much older @BELAIR1995317 gives $f$ as

fBelair _t m = a P./ (1 + k P.* (m P.** r))
  where
    a = 6570
    k = 0.0382
    r = 6.96

For the intial precursor population we can either use the
integral of the initial distribution

    result $ relative 1e-6 $ parTrap p'0Untyped 0.001 (muF /~ day)

> bigP'0 :: Quantity (DAmountOfSubstance / DMass) Double
> bigP'0 = r *~ (mole / kilo gram)
>  where
>    r = result $ relative 1e-6 $ parTrap p'0Untyped 0.001 (muF /~ day)

    [ghci]
    bigP'0

or we can use the sum of the values used in the finite difference approximation

> bigP'0' :: Quantity (DAmountOfSubstance / DMass) Double
> bigP'0' = (* deltaMu) $ sum $ map p'0 $ take n1 $ iterate (+ deltaMu) (0.0 *~ day)

    [ghci]
    bigP'0'

@Thibodeaux2011 give the following for $a_E$

> a_E :: Quantity (DAmountOfSubstance / DMass) Double -> Frequency Double
> a_E bigP = ((n / d) /~ one) *~ (one / day)
>   where
>     n :: Dimensionless Double
>     n = bigP * (13.8 *~ (kilo gram / mole)) + 0.04 *~ one
>     d :: Dimensionless Double
>     d = (bigP /~ (mole / kilo gram)) *~ one + 0.08 *~ one

*but* from @Ackleh200621

 > The only biological basis for the latter is that the decay rate of
 > erythropoietin should be an increasing function of the precursor
 > population and this function remains in the range 0.50–6.65
 > $\mathrm{days}^{-1}$

and, given this is at variance with their given function, it may be
safer to use their alternative of

> a_E' :: Quantity (DAmountOfSubstance / DMass) Double -> Frequency Double
> a_E' _bigP = 6.65 *~ (one / day)

> f'0 :: Quantity (DConcentration / DTime) Double
> f'0 = fAckleh undefined bigM'0

> bigE'1 :: Concentration Double
> bigE'1 = (bigE'0 + deltaT * f'0) / (1.0 *~ one + deltaT * a_E' bigP'0)

> d_1 :: Dimensionless Double ->
>        Concentration Double ->
>        Int ->
>        Dimensionless Double
> d_1 g e i = (1 *~ one) + (g * deltaT / deltaMu)
>           - deltaT * sigmaThibodeaux ((fromIntegral i *~ one) * deltaMu) undefined e

> d_2 :: Quantity (DAmountOfSubstance / DMass) Double ->
>        Int ->
>        Dimensionless Double
> d_2 bigM i = (1 *~ one) + deltaT / deltaNu
>            + deltaT * gammaThibodeaux ((fromIntegral i *~ one) * deltaNu) undefined bigM

> oneStepM :: (Matrix Double, Matrix Double, Concentration Double, Time Double) ->
>             Writer [(Quantity (DAmountOfSubstance / DMass) Double,
>                      Quantity (DAmountOfSubstance / DMass) Double,
>                      Concentration Double)]
>                    (Matrix Double, Matrix Double, Concentration Double, Time Double)
> oneStepM (psPrev, msPrev, ePrev, tPrev) = do
>   let
>     g  = gThibodeaux ePrev
>     ls = replicate n1 (negate $ g * deltaT / deltaMu)
>     ds = g : map (d_1 g ePrev)  [1..n1]
>     us = replicate n1 (0.0 *~ one)
>     b1'0 = (s_0 tPrev * ePrev) /~ (mole / second / kilo gram)
>     b1 = asColumn $ vjoin [scalar b1'0, subVector 1 n1 $ flatten psPrev]
>     psNew :: Matrix Double
>     psNew = triDiagSolve (fromList (map (/~ one) ls))
>                          (fromList (map (/~ one) ds))
>                          (fromList (map (/~ one) us))
>                          b1
>     ls2 = replicate n2 (negate $ deltaT / deltaNu)
>     bigM :: Quantity (DAmountOfSubstance / DMass) Double
>     bigM = (* deltaNu) $ ((sumElements msPrev) *~ (mole / kilo gram / second))
>     ds2 = (1.0 *~ one) : map (d_2 bigM) [1..n2]
>     us2 = replicate n2 (0.0 *~ one)
>     b2'0 = (g * (psNew `atIndex` (n1, 0) *~ (mole / second / kilo gram))) /~
>            (mole / second / kilo gram)
>     b2 = asColumn $ vjoin [scalar b2'0, subVector 1 n2 $ flatten msPrev]
>     msNew :: Matrix Double
>     msNew = triDiagSolve (fromList (map (/~ one) ls2))
>                          (fromList (map (/~ one) ds2))
>                          (fromList (map (/~ one) us2))
>                          b2
>     bigP :: Quantity (DAmountOfSubstance / DMass) Double
>     bigP = (* deltaMu) $ sumElements psPrev *~ (mole / kilo gram / second)
>     f :: Quantity (DConcentration / DTime) Double
>     f = fAckleh undefined bigM
>     eNew :: Concentration Double
>     eNew = (ePrev + deltaT * f) / (1.0 *~ one + deltaT * a_E' bigP)
>     tNew = tPrev + deltaT
>   tell [(bigP, bigM, ePrev)]
>   return (psNew, msNew, eNew, tNew)

> ys :: [(Quantity (DAmountOfSubstance / DMass) Double,
>         Quantity (DAmountOfSubstance / DMass) Double,
>         Concentration Double)]
> ys = take 2000 $
>      snd $
>      runWriter $
>      iterateM_ oneStepM ((((n1 P.+1 )><1) (map (/~ (mole / second / kilo gram)) b'0)),
>                          (((n2 P.+ 1)><1) $ (map (/~ (mole / second / kilo gram)) b_2'0)),
>                          bigE'0,
>                          (0.0 *~ day))

References
==========
