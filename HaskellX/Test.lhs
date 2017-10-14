\documentclass[presentation]{beamer}

%include polycode.fmt
%options ghci -fglasgow-exts

%format muPrior         = "\boldsymbol{\mu_0}"
%format sigmaPrior      = "\boldsymbol{\Sigma_0}"
%format bigH            = "\boldsymbol{H}"
%format bigHH           = "\boldsymbol{H}"
%format bigSigmaY       = "\boldsymbol{\Sigma}^{(y)}"
%format bigA            = "\boldsymbol{A}"
%format bigAA           = "\boldsymbol{A}"
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
%format S.<>            = "\mathbin{\texttt{<>}}"
%format S.#>            = "\mathbin{\texttt{\#>}}"
%format <>              = "\mathbin{\texttt{<>}}"
%format #>              = "\mathbin{\texttt{\#>}}"
%format mu0             = "\mu_0"
%format sigma0          = "\sigma_0"
%format x0              = "x_0"
%format sigma           = "\sigma"
%format xPrev           = "\boldsymbol{x}^\flat"
%format bigQ            = "\boldsymbol{Q}"
%format xNew            = "\boldsymbol{x}"
%format yNew            = "\boldsymbol{y}"
%format bigR            = "\boldsymbol{R}"
%format m0              = "\boldsymbol{m}_0"
%format bigSigma0       = "\boldsymbol{\Sigma}_0"
%format theta           = "\theta"
%format forall          = "\forall"

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
\author[Dominic Steinitz]{Dominic Steinitz \\ @@idontgetoutmuch \\ \url{https://idontgetoutmuch.wordpress.com}}
\hypersetup{
 pdfauthor={Dominic Steinitz},
 pdfkeywords={},
 pdflang={English}}

\usepackage{dramatist}

\renewcommand{\hscodestyle}{\small}

\newcommand{\blank}{\_}

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
    \item At time 1, I give you some information: the number with added noise.
    \item At time 2, I give you more information: the same number but with different added noise.
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


\begin{frame}{Prior}

%if style == newcode
\begin{code}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}

{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes          #-}

{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module Test ( kalmans'''
            , kalmans''
            , kalmans'
            , kalmans
            , m0X1
            , m1X1
            , m2X1
            , muX1
            , m0Y1
            , m1Y1
            , muY1
            , marginal1X
            , marginal1Y
            , varX1
            , x0
            , sigma
            , likelihood
            , ds
            , mu0
            , sigma0
            , hist
            , priors
            , main
            ) where

import Control.Monad
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State
import Control.Monad.Loops as ML

import qualified Data.Vector.Unboxed as V
import Data.Histogram.Fill
import qualified Data.Histogram as H
import Data.Histogram ( Histogram )

import qualified Control.Foldl as F

import           Numeric.LinearAlgebra.HMatrix ( (<>), tr, (#>), inv )
import qualified Numeric.LinearAlgebra.HMatrix as M
import qualified Numeric.LinearAlgebra.Static as S

import qualified Data.Random.Distribution.MultivariateNormal as G

import           GHC.TypeLits

import Diagrams.Prelude hiding ( normal, sample, (<>), inv, trace )
import Diagrams.Backend.Rasterific
import Plots hiding ( numBins, pdf )
import qualified Plots as P


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

hb :: F.Foldable f => f Double -> HBuilder Double (Histogram BinD Double)
hb xs = forceDouble -<< mkSimple (binD lower numBins upper)
  where
    (varX, xBar, _) = stats xs
    lower = xBar - 4.0 * sqrt varX
    upper = xBar + 4.0 * sqrt varX

hist :: F.Foldable f => f Double -> Histogram BinD Double
hist xs = fillBuilder (hb xs) xs
\end{code}
%endif

I give you the prior $\mathbb{P}(\theta)$ (was $\mathbb{P}(A)$) as a \textit{histogram}

\begin{code}
mu0, sigma0 :: Double
mu0 = 0.0
sigma0 = 1.00

priors :: Histogram BinD Double
priors = hist $ runSampler (normal mu0 sigma0) 2 100000
\end{code}

\note{

The histogram is generated from a normal distribution but as
you can see below it is not really normal.

}

\end{frame}

\begin{frame}{Prior}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/prior.png}
  \end{center}
\end{frame}

\begin{frame}{Secret, Model and Data}

I sample my secret value

\begin{code}
x0 :: Double
x0 = 1.00
\end{code}

\pause

I tell you the model aka the likelihood $\mathbb{P}(D \,\vert\,
\theta)$ (was $\mathbb{P}(A \,\vert\, B)$)

\begin{code}
sigma :: Double
sigma = 0.81

likelihood :: Double -> Double -> Double
likelihood bigD theta= n / d
  where
    n = exp (-(bigD - theta)^2 / (2 * sigma^2))
    d = sqrt (2 * sigma^2)
\end{code}

\pause

Finally I give you some noisy data

\begin{code}
ds :: [Double]
ds = runSampler (normal x0 sigma) 2 10
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

Now you can use Bayes' to determine my secret value

\begin{code}
posteriorize ::  Histogram BinD Double ->
                 Double ->
                 Histogram BinD Double
posteriorize q d = H.bmap bayes q
  where
    bayes theta p = p * likelihood d theta
\end{code}

\pause

\begin{code}
qs :: [Histogram BinD Double]
qs = scanl posteriorize priors ds
\end{code}

\pause

\begin{code}
ss :: [Double]
ss = map H.sum qs
\end{code}

\pause

\begin{code}
ns :: [Histogram BinD Double]
ns = zipWith (\s q -> H.map (/ s) q) ss qs
\end{code}
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

\pause

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
    \boldsymbol{Q}_k
    $$
  \end{block}

\pause

  \begin{block}{Can only observe position}
    $$
    \begin{bmatrix}y^{(k)}_1 \\ y^{(k)}_2\end{bmatrix} =
    \begin{bmatrix}
      1 & 0 & 0 & 0\\
      0 & 1 & 0 & 0
    \end{bmatrix}
    \begin{bmatrix}x^{(k-1)}_1 \\ x^{(k-1)}_2 \\ x^{(k-1)}_3 \\ x^{(k-1)}_4\end{bmatrix} +
    \boldsymbol{R}_k
    $$
  \end{block}

\pause

  \begin{block}{In vector notation}
    $$
    \begin{aligned}
    \boldsymbol{x}_k &= \boldsymbol{A} \boldsymbol{x}_{k-1} + \boldsymbol{Q}_k \\
    \boldsymbol{y}_k &= \boldsymbol{H} \boldsymbol{x}_{k}   + \boldsymbol{R}_k
    \end{aligned}
    $$
  \end{block}

\end{frame}

\section{Bayes for Pollen}

\begin{frame}{Play the Game}

%if style == newcode
\begin{code}
deltaT, sigma1, qc1 :: Double
deltaT = 0.1
sigma1 = 1/2
qc1 = 1

bigAl :: [Double]
bigAl = [1,  deltaT,
         0,       1]

bigAA :: M.Matrix Double
bigAA = (2 M.>< 2) bigAl

bigQl2 :: [Double]
bigQl2 = [qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
          qc1 * deltaT^2 / 2,       qc1 * deltaT]

bigQ :: M.Herm Double
bigQ = M.trustSym $ (2 M.>< 2) bigQl2

bigHH :: M.Matrix Double
bigHH = (1 M.>< 2) [1, 0]

bigR :: M.Herm Double
bigR = M.trustSym $ (1 M.>< 1) [sigma1^2]
\end{code}
%endif

\begin{block}{As Before}
We can play the same game as we did before in this more complicated
setting. The derivation is similar but much longer.
\end{block}

\pause

First I sample my secret value, in this case a path followed by the
pollen particle. I simultaneously create a sample of noisy data\footnote{Using the \texttt{random-fu} and \texttt{hmatrix} packages}.

\note{

I let the computer do the work; I only select the starting value which
happens to be $(0,1)$, that is starting $0$ with velocity $1.

I take the starting value

\begin{itemize}
  \item Update it according the state update part of the model
  \item Add some noise
  \item Take the state and provide a noisy observation of it
\end{itemize}

}

\pause

%{
%format G.Normal = "{\color{blue} G.Normal}"
%format bigAA    = "{\color{blue} \boldsymbol{A}}"
%format bigHH    = "{\color{blue} \boldsymbol{H}}"
%format bigQ     = "{\color{blue} \boldsymbol{Q}}"
%format bigR     = "{\color{blue} \boldsymbol{R}}"

\begin{code}
pollenSamples :: [(M.Vector Double, M.Vector Double)]
pollenSamples = evalState  (ML.unfoldrM pollenSample m0)
                           (pureMT 17)
  where
    pollenSample xPrev = do
      xNew <- sample $ rvar (G.Normal (bigAA #> xPrev) bigQ)
      yNew <- sample $ rvar (G.Normal (bigHH #> xNew) bigR)
      return $ Just ((xNew, yNew), xNew)
\end{code}
%}

\end{frame}

\begin{frame}{Prior}

I give you the prior $\mathbb{P}(A)$ again as a
\textit{histogram} but now 2D\footnote{we are just going to
track the $x$-axis so our state space is 2 dimensional}

%if style == newcode
\begin{code}
hb2 :: [(Double, Double)] ->
       HBuilder (Double, Double)
                (Histogram (Bin2D BinD BinD) Double)
hb2 xys = forceDouble -<< mkSimple (Bin2D (binD lowerX numBins upperX)
                                          (binD lowerY numBins upperY))
  where
    xs = map fst xys
    ys = map snd xys
    (varX, xBar, _) = stats xs
    lowerX = xBar - 5.0 * sqrt varX
    upperX = xBar + 5.0 * sqrt varX
    (varY, yBar, _) = stats ys
    lowerY = yBar - 5.0 * sqrt varY
    upperY = xBar + 5.0 * sqrt varY

hb2' :: [(Double, Double)] ->
       HBuilder (BinValue (Bin2D BinD BinD), Double)
                (Histogram (Bin2D BinD BinD) Double)
hb2' xys = forceDouble -<< mkWeighted (Bin2D (binD lowerX numBins upperX)
                                             (binD lowerY numBins upperY))
  where
    xs = map fst xys
    ys = map snd xys
    (varX, xBar, _) = stats xs
    lowerX = xBar - 5.0 * sqrt varX
    upperX = xBar + 5.0 * sqrt varX
    (varY, yBar, _) = stats ys
    lowerY = yBar - 5.0 * sqrt varY
    upperY = xBar + 5.0 * sqrt varY

hist2 :: [(Double, Double)] ->
         Histogram (Bin2D BinD BinD) Double
hist2 xys = fillBuilder (hb2 xys) xys

hist2' :: [((Double, Double), Double)] ->
         Histogram (Bin2D BinD BinD) Double
hist2' xys = fillBuilder (hb2' (map fst xys)) xys

conv :: [M.Vector Double] -> [(Double, Double)]
conv = map (\[x, y] -> (x, y)) .  map M.toList
\end{code}
%endif

\begin{code}
m0 :: M.Vector Double
m0 = M.fromList [0, 1]

bigSigma0 :: M.Herm Double
bigSigma0 = M.sym $ (2 M.>< 2)  [  1.0, 0.0,
                                   0.0, 1.0]

priorsPollen :: Histogram (Bin2D BinD BinD) Double
priorsPollen = hist2 $ conv $
               prePriors 2 100000
  where prePriors seed n =
          evalState  (replicateM n (sample $ G.Normal m0 bigSigma0))
                     (pureMT (fromIntegral seed))
\end{code}

\end{frame}

\begin{frame}{An Aside: Marginals}

\begin{code}
marginalX, marginalY :: Histogram BinD Double
marginalX = H.reduceX H.sum priorsPollen
marginalY = H.reduceY H.sum priorsPollen
\end{code}

\end{frame}

\begin{frame}{Marginal Prior for Position}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/marginalY.png}
  \end{center}
\end{frame}

\begin{frame}{Marginal Prior for Velocity}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/marginalX.png}
  \end{center}
\end{frame}

\begin{frame}{Model}

I give you the model

\begin{code}
newState :: MonadRandom m =>
             M.Vector Double -> m (M.Vector Double)
newState xPrev =
  sample $ rvar (G.Normal (bigAA #> xPrev) bigQ)
\end{code}

\begin{code}
weightK :: M.Vector Double -> M.Vector Double -> Double
weightK a b = pdf (G.Normal (bigHH #> a) bigR) b
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

Now you can use Bayes' to track the path of the pollen.

\begin{code}
posteriorizeK :: MonadRandom m =>
                 Histogram (Bin2D BinD BinD) Double ->
                 M.Vector Double ->
                 m (Histogram (Bin2D BinD BinD) Double)
posteriorizeK q d = do
  let xfs = H.asList q
      xs = map pair2Vec $ map fst xfs
      fs = map snd xfs
  xNew <- mapM newState xs
  let newQ = hist2' $ zip (map vec2Pair xNew) fs
  return $ H.bmap (\(u, v) p -> p * weightK (M.vector [u, v]) d) newQ
\end{code}

%if style == newcode
\begin{code}
vec2Pair :: M.Vector Double -> (Double, Double)
vec2Pair v = (v M.! 0, v M.! 1)

pair2Vec :: (Double, Double) -> M.Vector Double
pair2Vec (x, y) = M.vector [x, y]

test :: Histogram (Bin2D BinD BinD) Double
test = evalState (posteriorizeK priorsPollen (head $
                                              map (M.vector . pure . (!!0) . M.toList) $
                                              map fst $ take 10 pollenSamples))
                 (pureMT 42)
\end{code}
%endif

\end{frame}

\begin{frame}{Position Estimate After 1 Observation}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/marginal1Y.png}
  \end{center}
\end{frame}

\begin{frame}{Velocity Estimate After 1 Observation}
  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/marginal1X.png}
  \end{center}
\end{frame}

\begin{frame}{Marginals Again}

\begin{verbatim}
*Test> take 1 pollenSamples
[([9.99586e-2,0.75767],[-0.36189])]
\end{verbatim}

\begin{code}
marginal1X, marginal1Y :: Histogram BinD Double
marginal1X = H.reduceX H.sum test
marginal1Y = H.reduceY H.sum test

m0X1, m1X1, m2X1, muX1, varX1 :: Double
m0X1 = H.sum marginal1X
m1X1 = H.sum $ H.bmap (\f v -> v * f) marginal1X
\end{code}

%if style == newcode
\begin{code}
m0Y1, m1Y1, muY1 :: Double
m0Y1 = H.sum marginal1Y
m1Y1 = H.sum $ H.bmap (\f v -> v * f) marginal1Y
\end{code}
%endif

\begin{code}
muY1 = m1Y1 / m0Y1
muX1 = m1X1 / m0X1
\end{code}

\begin{verbatim}
*Test> muX1
0.98156
*Test> muY1
9.98856e-2
\end{verbatim}

\end{frame}

\section{The Real Kalman}

\begin{frame}{Recall the Model}

\begin{block}{Wait a Minute $\ldots$}
Compute power in 1968?
\end{block}

\pause

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

The model is \textbf{linear} and the errors are \textbf{Gaussian}

\end{frame}

\begin{frame}{Kalman Itself}

A \textbf{lot} of algebraic manipulation gives the optimal solution.

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

%if style == newcode
\begin{code}
kalmans :: (M.Field t, Num (M.Vector t)) =>
           M.Vector t
        -> M.Matrix t
        -> M.Matrix t
        -> M.Matrix t
        -> M.Matrix t
        -> M.Matrix t
        -> [M.Vector t]
        -> [(M.Vector t, M.Matrix t)]
\end{code}
%endif

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

%if style == newcode
\begin{code}
kalmans' :: (M.Field t, Num (M.Vector t)) =>
            M.Vector t
         -> M.Matrix t
         -> M.Matrix t
         -> M.Matrix t
         -> M.Matrix t
         -> M.Matrix t
         -> [M.Vector t]
         -> [(M.Vector t, M.Matrix t)]
\end{code}
%endif

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

\begin{code}
kalmans'' :: forall m n . (KnownNat m, KnownNat n) =>
              S.R n -> S.Sq n -> S.L m n -> S.Sq m -> S.Sq n -> S.Sq n -> [S.R m] ->
              [(S.R n, S.Sq n)]
kalmans'' muPrior sigmaPrior bigH bigSigmaY bigA bigSigmaX ys = scanl kalman (muPrior, sigmaPrior) ys
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

\section{Naperian Functors}

\begin{frame}{Even Better}

\begin{itemize}
  \item With Naperian functors we can do even better\footnote{$F \simeq \mathrm{Hom}(A, \blank)$}
  \item Come back next year for variational inference using Naperian functors
  \item See Jeremy Gibbons' "APLicative Programming with Naperian Functors"
  \item Haskell package \texttt{random-fu} for random variables
  \item Haskell package \texttt{histogram-fill} for histograms
  \item Haskell package \texttt{kalman} for kalman, extended kalman, unscented kalman, particle filtering and smoothing
  \item Haskell package \texttt{hmatrix} for statically typed matrices
\end{itemize}

%if style == newcode
\begin{code}
kSaxis :: [(Double, Double)] -> Axis B V2 Double
kSaxis xs = r2Axis &~ do
  linePlot' xs
  xMin .= Just (-5.0)
  xMax .= Just 5.0
  yMin .= Just 0.0
  yMax .= Just 4200.0

main :: IO ()
main = do
  renderRasterific "diagrams/marginalX.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ H.asList marginalX)
  renderRasterific "diagrams/marginalY.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ H.asList marginalY)
  renderRasterific "diagrams/marginal1X.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ H.asList marginal1X)
  renderRasterific "diagrams/marginal1Y.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ H.asList marginal1Y)
\end{code}
%endif

\end{frame}

\end{document}
