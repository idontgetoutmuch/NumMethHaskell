\documentclass[presentation]{beamer}

%include polycode.fmt
%options ghci -fglasgow-exts

%format muPrior         = "\boldsymbol{\mu_0}"
%format sigmaPrior      = "\boldsymbol{\Sigma_0}"
%format bigH            = "\boldsymbol{H}"
%format bigSigmaY       = "\boldsymbol{\Sigma}^{(y)}"
%format bigA            = "\boldsymbol{A}"
%format bigSigmaX       = "\boldsymbol{\Sigma}^{(x)}"
%format xHat            = "\hat{\boldsymbol{x}}"
%format xHatFlat        = "\hat{\boldsymbol{x}}^\flat"
%format bigS            = "\boldsymbol{S}"
%format bigSigmaHat     = "\hat{\boldsymbol{\Sigma}}"
%format bigSigmaHatFlat = "\hat{\boldsymbol{\Sigma}}^\flat"
%format xHatFlatNew     = "\hat{\boldsymbol{x}}^\flat_{\mathrm{new}}"
%format bigSigmaHatFlatNew = "\hat{\boldsymbol{\Sigma}}^\flat_{\mathrm{new}}"
%format bigK            = "\boldsymbol{K}"
%format vv              = "\hat{\boldsymbol{v}}"
%format yy              = "\hat{\boldsymbol{y}}"


\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{fixltx2e}
\usepackage{graphicx}
\usepackage{grffile}
\usepackage{longtable}
\usepackage{wrapfig}
\usepackage{rotating}
\usepackage[normalem]{ulem}
\usepackage{amsmath}
\usepackage{textcomp}
\usepackage{amssymb}
\usepackage{capt-of}
\usepackage{hyperref}
\usepackage{listings}
\usepackage{color}
\usepackage{verse}
\RequirePackage{fancyvrb}
\DefineVerbatimEnvironment{verbatim}{Verbatim}{fontsize=\scriptsize}
\usepackage[style=alphabetic]{biblatex}
\usetheme{Frankfurt}
\author{Dominic Steinitz}
\date{Thursday 12 October 17}
\title{Making Kalman Filtering Correct with Types}
\hypersetup{
 pdfauthor={Dominic Steinitz},
 pdfkeywords={},
 pdflang={English}}

\usepackage{dramatist}

\renewcommand{\hscodestyle}{\small}

\begin{document}

\maketitle
\begin{frame}{Outline}
\tableofcontents
\end{frame}

\section{Introduction}

\begin{frame}{Apollo 8 launched on December 21, 1968}

\StageDir{03:17:45:17 (Dec. 25, 1968, 6:36 a.m. UTC)}

\begin{drama}
  \Character{Jim Lovell (Commander Module Pilot)}{jim}
  \Character{Ken Mattingly (CAPCOM)}{ken}

  \jimspeaks: Roger. Do you wish me to reinitialize the W-matrix at this time?
\end{drama}

\StageDir{03:17:45:26}

\begin{drama}
  \Character{Jim Lovell (Commander Module Pilot)}{jim}
  \Character{Ken Mattingly (CAPCOM)}{ken}

  \kenspeaks: Affirmative, Apollo 8
\end{drama}

\section{Introducing the Reverend Bayes}
\end{frame}

\begin{frame}{Game}

  \begin{block}{Game}
    \begin{itemize}
    \item I select a number at random from a normal distribution.
    \item At time 1 I give you some information: the number with added noise.
    \item At time 2 I give you more information: the same number but with different added noise.
    \item And so on $\ldots$
    \end{itemize}
  \end{block}

\end{frame}

\begin{frame}{Bayes' Theorem}

  $$
  \mathbb{P}(A \,\vert\, B) \triangleq \frac{\mathbb{P}(A \cap B)}{\mathbb{P}(B)}
  $$

  Also

  $$
  \mathbb{P}(B \,\vert\, A) \triangleq \frac{\mathbb{P}(A \cap B)}{\mathbb{P}(A)}
  $$

  Thus

  $$
  \mathbb{P}(A \,\vert\, B) \propto {\mathbb{P}(B \,\vert\, A)}{\mathbb{P}(A)}
  $$

\end{frame}


\begin{frame}{Prior and Secret}

%if style == newcode
\begin{code}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}

{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes          #-}

{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

import Control.Monad
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State
import Control.Monad.Loops as ML

import qualified Data.Vector.Unboxed as V
import Data.Histogram.Fill
import qualified Data.Histogram as H
import Data.Histogram.Generic ( Histogram )

import qualified Control.Foldl as F

import           Numeric.LinearAlgebra.HMatrix ( (<>), tr, (#>), inv )
import qualified Numeric.LinearAlgebra.HMatrix as M
import qualified Numeric.LinearAlgebra.Static as S

import qualified Data.Random.Distribution.MultivariateNormal as G

import           GHC.TypeLits

numBins :: Int
numBins = 100

runSampler :: RVar Double -> Int -> Int -> [Double]
runSampler sampler seed n =
  evalState (sample (replicateM n sampler))
             (pureMT (fromIntegral seed))

stats :: (F.Foldable f, Fractional a) =>
         f a -> (a, a, a)
stats v = F.fold stats' v
  where
    stats' = f <$> (F.premap (\x -> x * x) F.sum) <*> F.sum <*> F.genericLength
    f x2Sum xSum n = (var, mean, n)
      where
        mean = xSum / n
        mean2 = x2Sum / n
        var = n * (mean2 - mean * mean) / (n - 1)

hb :: F.Foldable f => f Double -> HBuilder Double (Histogram V.Vector BinD Double)
hb xs = forceDouble -<< mkSimple (binD lower numBins upper)
  where
    (varX, xBar, _) = stats xs
    lower = xBar - 4.0 * sqrt varX
    upper = xBar + 4.0 * sqrt varX

hist :: F.Foldable f => f Double -> Histogram V.Vector BinD Double
hist xs = fillBuilder (hb xs) xs
\end{code}
%endif

I give you the prior $\mathbb{P}(A)$

\begin{code}
mu0, sigma0 :: Double
mu0 = 0.0
sigma0 = 1.00

priors :: Histogram V.Vector BinD Double
priors = hist $ runSampler (normal mu0 sigma0) 2 100000
\end{code}

I sample my secret value

\begin{code}
x0 :: Double
x0 = 1.00
\end{code}

\end{frame}

\begin{frame}{Model and Data}

I tell you the model aka the likelihood $\mathbb{P}(A \,\vert\, B)$

\begin{code}
sigma :: Double
sigma = 0.81

likelihood :: Double -> Double -> Double
likelihood a b = n / d
  where
    n = exp (-(a-b)^2 / (2 * sigma^2))
    d = sqrt (2 * sigma^2)
\end{code}

Finally I give you some noisy data

\begin{code}
ds :: [Double]
ds = runSampler (normal x0 sigma) 2 10
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

Now you can use Bayes' to determine my secret value

\begin{code}
posteriorize :: (BinValue bin ~ Double, Bin bin) =>
                H.Histogram bin Double ->
                Double ->
                H.Histogram bin Double
posteriorize q d = H.bmap (\v p -> p * likelihood v d) q

qs :: [H.Histogram BinD Double]
qs = scanl posteriorize priors ds

ss :: [Double]
ss = map H.sum qs

ns :: [H.Histogram BinD Double]
ns = zipWith (\s q -> H.map (/ s) q) ss qs
\end{code}
\end{frame}

\begin{frame}{Prior}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/prior.png}
  \end{center}
\end{frame}

\begin{frame}{After 1 Observation}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/qs1.png}
  \end{center}
\end{frame}

\begin{frame}{After 10 Observations}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/qsN.png}
  \end{center}
\end{frame}

\section{Introducing the Reverend Brown}

\begin{frame}{Robert Brown (1827)}

  \begin{itemize}
  \item You wish to emulate the famous botanist but with a difference.
  \item You have a camera which gives approximate co-ordinates of the
    pollen particle on the slide.
  \item You have a motor which can drive the slide in horizontal and
    vertical planes.
  \item How to track the camera to minimize the particle's distance
    from the centre of the microscope's field of vision?
  \end{itemize}

\end{frame}

\begin{frame}{Mathematical Model}
  We can model the pollen's motion as

  \begin{block}{Equations of Motion}
    $$
    \begin{aligned}
      \frac{\mathrm{d}^2 x_1}{\mathrm{d}t^2} &= \omega_1(t) \\
      \frac{\mathrm{d}^2 x_2}{\mathrm{d}t^2} &= \omega_2(t)
    \end{aligned}
    $$
  \end{block}

  Writing $x_3 = \mathrm{d}x_1 / \mathrm{d}t$ and
  $x_4 = \mathrm{d}x_2 / \mathrm{d}t$ this becomes

  \begin{block}{Matrix Form}
  $$
  \frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix}x_1 \\ x_2 \\ x_3 \\ x_4\end{bmatrix} =
  \begin{bmatrix}
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & 1 \\
    0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0
  \end{bmatrix}
  \begin{bmatrix}x_1 \\ x_2 \\ x_3 \\ x_4\end{bmatrix} +
  \begin{bmatrix}
    0 & 0 \\
    0 & 0 \\
    1 & 0 \\
    0 & 1
  \end{bmatrix}
  \begin{bmatrix}\omega_1 \\ \omega_2\end{bmatrix}
  $$
  \end{block}

\end{frame}

\begin{frame}{}

  \begin{block}{Discretizing at $0, \Delta t, 2\Delta t, \ldots $}
    $$
    \begin{bmatrix}x^{(k)}_1 \\ x^{(k)}_2 \\ x^{(k)}_3 \\ x^{(k)}_4\end{bmatrix} =
    \begin{bmatrix}
      1 & 0 & \Delta t & 0 \\
      0 & 1 & 0        & \Delta t \\
      0 & 0 & 1        & 0 \\
      0 & 0 & 0        & 1
    \end{bmatrix}
    \begin{bmatrix}x^{(k-1)}_1 \\ x^{(k-1)}_2 \\ x^{(k-1)}_3 \\ x^{(k-1)}_4\end{bmatrix} +
    \mathbf{q}_k
    $$
  \end{block}

  \begin{block}{In vector notation}
    $$
    \mathbf{x}_k = \mathbf{A} \mathbf{x}_{k_1} + \mathbf{q}_k
    $$
  \end{block}

\end{frame}

\section{Bayes for Pollen}

\begin{frame}{Play the Game}

%if style == newcode
\begin{code}
deltaT, sigma1, sigma2, qc1, qc2 :: Double
deltaT = 0.1
sigma1 = 1/2
sigma2 = 1/2
qc1 = 1
qc2 = 1

bigAl :: [Double]
bigAl = [1, 0 , deltaT,      0,
         0, 1,       0, deltaT,
         0, 0,       1,      0,
         0, 0,       0,      1]

bigAl2 :: [Double]
bigAl2 = [1,  deltaT,
          0,       1]

bigA :: M.Matrix Double
bigA = (4 M.>< 4) bigAl

bigA2 :: M.Matrix Double
bigA2 = (2 M.>< 2) bigAl2

bigQl :: [Double]
bigQl = [qc1 * deltaT^3 / 3,                  0, qc1 * deltaT^2 / 2,                  0,
                          0, qc2 * deltaT^3 / 3,                  0, qc2 * deltaT^2 / 2,
         qc1 * deltaT^2 / 2,                  0,       qc1 * deltaT,                  0,
                          0, qc2 * deltaT^2 / 2,                  0,       qc2 * deltaT]

bigQl2 :: [Double]
bigQl2 = [qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
          qc1 * deltaT^2 / 2,       qc1 * deltaT]

bigQ :: M.Herm Double
bigQ = M.trustSym $ (4 M.>< 4) bigQl

bigQ2 :: M.Herm Double
bigQ2 = M.trustSym $ (2 M.>< 2) bigQl2

newState :: MonadRandom m => M.Vector Double -> m (M.Vector Double)
newState xPrev = sample $ rvar (G.Normal (bigA M.#> xPrev) bigQ)

bigH2 :: M.Matrix Double
bigH2 = (1 M.>< 2) [1, 0]

bigR :: M.Herm Double
bigR = M.trustSym $ (2 M.>< 2) [sigma1^2,        0,
                                       0, sigma2^2]

bigR2 :: M.Herm Double
bigR2 = M.trustSym $ (1 M.>< 1) [sigma1^2]

m02 :: M.Vector Double
m02 = M.fromList [0, 1]
\end{code}
%endif

\begin{block}{As Before}
We can play the same game as we did before in this more complicated
setting. The derivation is similar but much longer.
\end{block}

First I sample my secret value, in this case a path followed by the
pollen particle. I simultaneously create a sample of noisy data.

\begin{code}
pollenSample :: MonadRandom m =>
             M.Vector Double ->
             m (Maybe ((M.Vector Double, M.Vector Double),
                M.Vector Double))
pollenSample xPrev = do
  xNew <- sample $ rvar (G.Normal (bigA2 M.#> xPrev) bigQ2)
  yNew <- sample $ rvar (G.Normal (bigH2 M.#> xNew) bigR2)
  return $ Just ((xNew, yNew), xNew)

pollenSamples :: [(M.Vector Double, M.Vector Double)]
pollenSamples = evalState (ML.unfoldrM pollenSample m02)
                          (pureMT 17)
\end{code}

\end{frame}

\begin{frame}{Prior}

I give you the prior $\mathbb{P}(A)$

%if style == newcode
\begin{code}
hb2 :: [(Double, Double)] ->
       HBuilder (Double, Double)
                (Histogram V.Vector (Bin2D BinD BinD) Double)
hb2 xys = forceDouble -<< mkSimple (Bin2D (binD lowerX numBins upperX)
                                          (binD lowerY numBins upperY))
  where
    xs = map fst xys
    ys = map snd xys
    (varX, xBar, _) = stats xs
    lowerX = xBar - 4.0 * sqrt varX
    upperX = xBar + 4.0 * sqrt varX
    (varY, yBar, _) = stats ys
    lowerY = yBar - 4.0 * sqrt varY
    upperY = xBar + 4.0 * sqrt varY

hist2 :: [(Double, Double)] ->
         Histogram V.Vector (Bin2D BinD BinD) Double
hist2 xys = fillBuilder (hb2 xys) xys

conv :: [M.Vector Double] -> [(Double, Double)]
conv = map (\[x, y] -> (x, y)) .  map M.toList
\end{code}
%endif

\begin{code}
bigSigma0 :: M.Herm Double
bigSigma0 = M.sym $ (2 M.>< 2) [1.0, 0.0,
                                0.0, 1.0]

prePriors :: Int -> Int -> [M.Vector Double]
prePriors seed n =
  evalState (replicateM n (sample $ G.Normal m02 bigSigma0))
            (pureMT (fromIntegral seed))


priorsPollen :: Histogram V.Vector (Bin2D BinD BinD) Double
priorsPollen = hist2 $ conv $
               prePriors 2 100000
\end{code}

\end{frame}

\begin{frame}{An Aside: Marginals}

\begin{code}
marginalX, marginalY :: H.Histogram BinD Double
marginalX = H.reduceX H.sum priorsPollen
marginalY = H.reduceY H.sum priorsPollen

m0X0 = H.sum marginalX
m1X0 = H.sum $ H.bmap (\f v -> v * f) marginalX
m2X0 = H.sum $ H.bmap (\f v -> v^2 * f) marginalX
muX0 = m1X0 / m0X0
varX0 = (m2X0 / m0X0) - muX0^2
\end{code}

\end{frame}

\begin{frame}{Model}

I give you the model

\begin{code}
newState2 :: MonadRandom m => M.Vector Double -> m (M.Vector Double)
newState2 xPrev = sample $ rvar (G.Normal (bigA2 M.#> xPrev) bigQ2)
\end{code}

\begin{code}
weightK :: M.Vector Double -> M.Vector Double -> Double
weightK a b = pdf (G.Normal (bigH2 M.#> a) bigR) b
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

Now you can use Bayes' to track the path of the pollen.

\begin{code}
posteriorizeK :: MonadRandom m =>
                 H.Histogram (Bin2D BinD BinD) Double ->
                 M.Vector Double ->
                 m (H.Histogram (Bin2D BinD BinD) Double)
posteriorizeK q d = do
  u <- mapM newState2 $ map (\(x, y) -> M.vector [x, y]) $
       map fst $ H.asList q
  let newQ = hist2 $ map (\v -> (v M.! 0, v M.! 1)) u
  return $ H.bmap (\(u, v) p -> p * weightK (M.vector [u, v]) d) newQ
\end{code}

%if style == newcode
\begin{code}
test = evalState (posteriorizeK priorsPollen (head $ map fst $ take 10 pollenSamples))
                 (pureMT 42)

marginal1X, marginal1Y :: H.Histogram BinD Double
marginal1X = H.reduceX H.sum test
marginal1Y = H.reduceY H.sum test

m0X1 = H.sum marginal1X
m1X1 = H.sum $ H.bmap (\f v -> v * f) marginal1X
m2X1 = H.sum $ H.bmap (\f v -> v^2 * f) marginal1X
muX1 = m1X1 / m0X1
varX1 = (m2X1 / m0X1) - muX1^2
\end{code}
%endif

\end{frame}

\section{The Real Kalman}

\begin{frame}{Recall the Model}

Compute power in 1968?

Recall our model:

$$
\begin{aligned}
\boldsymbol{x}_i &= \boldsymbol{A}_{i-1}\boldsymbol{x}_{i-1} + \boldsymbol{\psi}_{i-1} \\
\boldsymbol{y}_i &= \boldsymbol{H}_i\boldsymbol{x}_i + \boldsymbol{\upsilon}_i
\end{aligned}
$$

where

$$
\boldsymbol{\psi}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(x)}_i\big)
,\quad
\boldsymbol{\upsilon}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(y)}_i\big)
$$

\end{frame}

\begin{frame}{Kalman Itself}

A \textbf{lot} of algebraic manipulation give the optimal solution.

\begin{block}{Prediction Step}
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
\end{block}

\begin{block}{Correction Step}
$$
\begin{aligned}
\boldsymbol{v}_i & =
\boldsymbol{y}_i - \boldsymbol{H}_i\hat{\boldsymbol{x}}^\flat_i \\
\boldsymbol{S}_i & =
\boldsymbol{H}_i \hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top + \boldsymbol{\Sigma}^{(y)}_i \\
\boldsymbol{K}_i & = \hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top\boldsymbol{S}^{-1}_i \\
\hat{\boldsymbol{x}}^i &= \hat{\boldsymbol{x}}^\flat_i + \boldsymbol{K}_i\boldsymbol{v}_i \\
\hat{\boldsymbol{\Sigma}}_i &= \hat{\boldsymbol{\Sigma}}^\flat_i - \boldsymbol{K}_i\boldsymbol{S}_i\boldsymbol{K}^\top_i
\end{aligned}
$$
\end{block}

\end{frame}

\begin{frame}{Kalman in Haskell}

\begin{code}
kalmans muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = scanl kalman (muPrior, sigmaPrior) ys
  where
    kalman (xHatFlat, bigSigmaHatFlat) yy = (xHatFlatNew, bigSigmaHatFlatNew)
      where
        vv = yy - bigH #> xHatFlat
        bigS = bigH <> bigSigmaHatFlat <> (tr bigH) + bigSigmaY
        bigK = bigSigmaHatFlat <> (tr bigH) <> (inv bigS)
        xHat = xHatFlat + bigK #> vv
        bigSigmaHat = bigSigmaHatFlat - bigK <> bigS <> (tr bigK)
        xHatFlatNew = bigA #> xHat
        bigSigmaHatFlatNew = bigA <> bigSigmaHat <> (tr bigA) + bigSigmaX
\end{code}

\end{frame}

\begin{frame}{With Comments}

\begin{code}
-- R n -> Sq n -> L m n -> Sq m -> Sq n -> Sq n -> [R m] ->
-- [(R n, Sq n)]
kalmans' muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = scanl kalman (muPrior, sigmaPrior) ys
  where
    -- kalman :: (R n, Sq n) -> R m -> (R n, Sq n)
    kalman (xHatFlat, bigSigmaHatFlat) yy = (xHatFlatNew, bigSigmaHatFlatNew)
      where
        vv = yy - bigH #> xHatFlat -- R m
        bigS = bigH <> bigSigmaHatFlat <> (tr bigH) + bigSigmaY -- Sq m
        bigK = bigSigmaHatFlat <> (tr bigH) <> (inv bigS) -- L n m
        xHat = xHatFlat + bigK #> vv -- R n
        bigSigmaHat = bigSigmaHatFlat - bigK <> bigS <> (tr bigK) -- Sq n
        xHatFlatNew = bigA #> xHat -- R n
        bigSigmaHatFlatNew = bigA <> bigSigmaHat <> (tr bigA) + bigSigmaX -- Sq n
\end{code}

\end{frame}

\begin{frame}{Why Walk?}

%if style == newcode
\begin{code}
kalmans'' :: forall m n . (KnownNat m, KnownNat n) =>
             S.R n -> S.Sq n -> S.L m n -> S.Sq m -> S.Sq n -> S.Sq n -> [S.R m] ->
             [(S.R n, S.Sq n)]
\end{code}
%endif

\begin{code}
kalmans'' muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = scanl kalman (muPrior, sigmaPrior) ys
  where
    kalman (xHatFlat, bigSigmaHatFlat) yy = (xHatFlatNew, bigSigmaHatFlatNew)
      where
        vv = yy - bigH S.#> xHatFlat
        bigS = bigH S.<> bigSigmaHatFlat S.<> (tr bigH) + bigSigmaY
        bigK = bigSigmaHatFlat S.<> (tr bigH) S.<> (S.inv bigS)
        xHat = xHatFlat + bigK S.#> vv
        bigSigmaHat = bigSigmaHatFlat - bigK S.<> bigS S.<> (S.tr bigK)
        xHatFlatNew = bigA S.#> xHat
        bigSigmaHatFlatNew = bigA S.<> bigSigmaHat S.<> (S.tr bigA) + bigSigmaX
\end{code}

\end{frame}

\begin{frame}{Why Walk?}

\begin{code}
kalmans''' :: forall m n . (KnownNat m, KnownNat n) =>
              S.R n -> S.Sq n -> S.L m n -> S.Sq m -> S.Sq n -> S.Sq n -> [S.R m] ->
              [(S.R n, S.Sq n)]
kalmans''' muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = scanl kalman (muPrior, sigmaPrior) ys
  where
    kalman :: (S.R n, S.L n n) -> S.R m -> (S.R n, S.L n n)
    kalman (xHatFlat, bigSigmaHatFlat) yy = (xHatFlatNew, bigSigmaHatFlatNew)
      where
        vv = yy - bigH S.#> xHatFlat
        bigS = bigH S.<> bigSigmaHatFlat S.<> (tr bigH) + bigSigmaY
        bigK = bigSigmaHatFlat S.<> (tr bigH) S.<> (S.inv bigS)
        xHat = xHatFlat + bigK S.#> vv
        bigSigmaHat = bigSigmaHatFlat - bigK S.<> bigS S.<> (S.tr bigK)
        xHatFlatNew = bigA S.#> xHat
        bigSigmaHatFlatNew = bigA S.<> bigSigmaHat S.<> (S.tr bigA) + bigSigmaX
\end{code}

\end{frame}

\end{document}
