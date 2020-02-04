\documentclass[handout]{beamer}
\usepackage{pgfpages}

\mode<handout>{%
    \pgfpagesuselayout{4 on 1}[a4paper]
    \setbeameroption{show notes}
}

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
%format xInit           = "\boldsymbol{x}_{\mathrm{init}}"
%format sigma0          = "\sigma_0"
%format mu              = "\mu"
%% This is so ghc doesn't complain but we get the right symbol on the slide
%format nu              = "\mu"
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
%format bigD            = "D"
%format pi              = "\pi"

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
\date{Thursday 6 February 2020}
\title{Some Thoughts on Filtering}
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

\note{I believe the W-matrix is the noise for the state update for the
 Kalman filter but have been unable to get a definitive reference. I
 further think the onboard computer had very low precision, possibly 8
 bits.}

\end{frame}

\section{Bayes 101}

\begin{frame}{Setting}

\note{Before getting into the Kalman filter which I am not going to
 derive anyway, let's take a very simple example.}

  \begin{block}{Samples}
    \begin{itemize}

    \item I can sample from a normal distribution with the mean known
    only to me and with a publicly known variance.

    \item Based on your knowledge about me you have a (prior) belief
    about the mean of this distribution: it's centred at 0 and is
    probably between -3 and +3.

    \item At time 1, I give you some information: a sample from my
    distribution.

    \item At time 2, I give you more information: another sample from
    my distribution.

    \item And so on $\ldots$
    \end{itemize}
  \end{block}

\end{frame}

\begin{frame}{Bayes' Theorem}

  $$
  \mathbb{P}(\theta \,\vert\, D) \triangleq \frac{\mathbb{P}(\theta \cap D)}{\mathbb{P}(D)}
  $$

  Also

  $$
  \mathbb{P}(D \,\vert\, \theta) \triangleq \frac{\mathbb{P}(\theta \cap D)}{\mathbb{P}(\theta)}
  $$

  Thus

  $$
  \mathbb{P}(\theta \,\vert\, D) \propto {\mathbb{P}(D \,\vert\, \theta)}{\mathbb{P}(\theta)}
  $$

\end{frame}


\begin{frame}{An Application of Bayes' Theorem}

%if style == newcode
\begin{code}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}

{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RankNTypes          #-}

{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module Test ( m0X1
            , m1X1
            , muX1
            , m0Y1
            , m1Y1
            , muY1
            , marginal1X
            , marginal1Y
            , mu
            , sigma
            , likelihood
            , ds
            , mu0
            , sigma0
            , hist
            , priors
            , ns
            , main
            ) where

import Control.Monad
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State
import Control.Monad.Loops as ML

import Data.Histogram.Fill hiding ( (><) )
import qualified Data.Histogram as H
import Data.Histogram ( Histogram )

import qualified Control.Foldl as F

import           Numeric.LinearAlgebra.HMatrix ( (#>), (><) )
import qualified Numeric.LinearAlgebra.HMatrix as M

import qualified Data.Random.Distribution.MultivariateNormal as G

import Diagrams.Prelude hiding ( normal, sample, trace )
import Diagrams.Backend.Rasterific
import Plots hiding ( numBins, pdf )


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

I have told you that the data are drawn from

$$
\mathbb{P}(x \, \vert \, \mu) \propto \exp{-\frac{(x - \mu)^2}{2\sigma^2}}
$$

And, based on your knowledge about me, you have a prior

$$
\mathbb{P}(\mu) \propto \exp{-\frac{(\mu - \mu_0)^2}{2\sigma_0^2}}
$$

This gives

$$
\mathbb{P}(\mu \, \vert \, x) \propto
\exp{-\frac{(x - \mu)^2}{2\sigma^2}} \times
\exp{-\frac{(\mu - \mu_0)^2}{2\sigma_0^2}}
$$

\end{frame}

\begin{frame}{An Application of Bayes' Theorem}

\begin{align*}
\mathbb{P}(\mu \, \vert \, x) & \propto
\exp{-\frac{(x - \mu)^2}{2\sigma^2}} \times
\exp{-\frac{(\mu - \mu_0)^2}{2\sigma_0^2}} \\
 \onslide<2->{& = \exp{-[\frac{x^2 - 2\mu x + \mu^2}{2\sigma^2} +
\frac{\mu^2 - 2\mu\mu_0 + \mu_0^2}{2\sigma_0^2}]}} \\
 \onslide<3->{& = \exp{-[\mu^2(\frac{1}{2\sigma^2} + \frac{1}{2\sigma_0^2}) -
     2\mu(\frac{x}{2\sigma^2} + \frac{\mu_0}{2\sigma_0^2}) +
     (\frac{x^2}{2\sigma^2} + \frac{\mu_0^2}{2\sigma_0^2})]}} \\
\onslide<4->{&\triangleq \exp{-[\frac{1}{2\sigma_1^2}(\mu^2 -2\mu\mu_1 + \mu_1^2)]}}
\end{align*}

\end{frame}

\begin{frame}{An Application of Bayes' Theorem}

\begin{align*}
\exp{-[\mu^2(\frac{1}{2\sigma^2} + \frac{1}{2\sigma_0^2}) -
     2\mu(\frac{x}{2\sigma^2} + \frac{\mu_0}{2\sigma_0^2}) +
     (\frac{x^2}{2\sigma^2} + \frac{\mu_0^2}{2\sigma_0^2})]} &\triangleq \\
\exp{-[\frac{1}{2\sigma_1^2}(\mu^2 -2\mu\mu_1 + \mu_1^2)]} &
\end{align*}

Collecting together terms

\begin{align*}
\onslide<2->{\frac{1}{\sigma_1^2} = \frac{1}{\sigma_0^2} + \frac{1}{\sigma^2}} \\
\onslide<3->{\frac{\mu_1}{\sigma_1^2} = \frac{\mu_0}{\sigma_0^2} + \frac{x}{\sigma^2}}
\end{align*}

\note{Note this required a fair amount of algebraic manipulation albeit only that taught routinely to 16 year olds (completing the square) and there is a danger that the wood can't be seen for the trees in these sorts of extended manipulations.}

\end{frame}

\begin{frame}{Prior by Sampling}

\note{With the power of a computer to hand we can take a completely na\"{i}ve approach.}

We can also represent the prior as a \textit{histogram}

\begin{code}
mu0, sigma0 :: Double
mu0 = 0.0
sigma0 = 1.00

priors :: Histogram BinD Double
priors = hist $ runSampler (normal mu0 sigma0) 42 100000
\end{code}

\end{frame}

\begin{frame}{Prior}

\note{ The histogram is generated from a normal distribution but as
  you can see, it is not really normal.}

  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/prior.png}
  \end{center}
\end{frame}

\begin{frame}{Secret, Model and Data}

But in reality the data come from a normal distribution with mean and (known) variance

\begin{code}
mu :: Double
mu = 1.00

sigma :: Double
sigma = 0.9
\end{code}

\pause

The function in Haskell for the likelihood $\mathbb{P}(x \,\vert\,
\mu)$

\begin{code}
likelihood :: Double -> Double -> Double
likelihood x nu = n / d
  where
    n = exp (-(x - nu)^2 / (2 * sigma^2))
    d = sqrt (2 * pi * sigma^2)
\end{code}

\pause

Finally I give you some noisy data

\begin{code}
ds :: [Double]
ds = runSampler (normal mu sigma) 2 10
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

\note{Of course the posterior has to be normalized.}

Now you can use Bayes' Theorem to create a posterior

\begin{code}
posteriorize ::  Histogram BinD Double ->
                 Double ->
                 Histogram BinD Double
posteriorize h x = H.bmap bayes h
  where
    bayes :: Double -> Double -> Double
    bayes theta p = p * likelihood x theta
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

\note{After 1 observation, we haven't really been persuaded to move away from our prior.}

  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/qs1.png}
  \end{center}
\end{frame}

\begin{frame}{After 10 Observations}

\note{But after 10 observations, we are getting pretty close to the actual distribution.}

  \begin{center}
    \includegraphics[height=0.80\textheight]{./diagrams/qsN.png}
  \end{center}
\end{frame}

\section{Apollo}

\begin{frame}{Observing Apollo (in 1D)}

We want to track the lunar module

  \begin{block}{Newton's second law of motion}

  $$
  ma = F
  $$
  \end{block}

\pause

  \begin{block}{As first order equation system}
    \begin{align*}
    \frac{\mathrm{d}v}{\mathrm{d}t} &= F \\
    \frac{\mathrm{d}x}{\mathrm{d}t} &= v
    \end{align*}
  \end{block}

\end{frame}

\begin{frame}{Mathematical Model}

Let's assume there are no forces acting on the lunar module.

  \begin{block}{Discretizing}
    \begin{align*}
    \frac{v_k - v_{k-1}}{\Delta t} &= 0 \\
    \frac{x_k - x_{k-1}}{\Delta t} &= v_k
    \end{align*}
  \end{block}

\note{Our model is not perfect so let's model the imperfection by
 adding in some white noise.}

  \begin{block}{Writing $x_1$ for $x$ and $x_2$ for $v$}
    $$
    \onslide<2->{
    \begin{bmatrix}x^{(k)}_1 \\ x^{(k)}_2\end{bmatrix} =
    \begin{bmatrix}
      1 & \Delta t \\
      0 & 1
    \end{bmatrix}
    \begin{bmatrix}x^{(k-1)}_1 \\ x^{(k-1)}_2\end{bmatrix}} \onslide<3->{+
    \boldsymbol{Q}_k}
    $$
  \end{block}

\end{frame}

\begin{frame}{Can only observe position}

\note{And of course our observations are also noisy.}

  \begin{block}{Observing}
    \onslide<1->{
    $$
    \begin{bmatrix}y^{(k)}_1\end{bmatrix} =
    \begin{bmatrix}
      1 & 0
    \end{bmatrix}
    \begin{bmatrix}x^{(k)}_1 \\ x^{(k)}_2\end{bmatrix}} \onslide<2->{+
    \boldsymbol{R}_k}
    $$
  \end{block}
\end{frame}

\begin{frame}
  \begin{block}{In vector notation}
    $$
    \begin{aligned}
    \boldsymbol{x}_k &= \boldsymbol{A} \boldsymbol{x}_{k-1} + \boldsymbol{Q}_k \\
    \boldsymbol{y}_k &= \boldsymbol{H} \boldsymbol{x}_{k}   + \boldsymbol{R}_k
    \end{aligned}
    $$
  \end{block}

\end{frame}

\section{Bayes for Apollo}

\begin{frame}{Sampling}

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
bigAA = (2 >< 2) bigAl

bigQl2 :: [Double]
bigQl2 = [qc1 * deltaT^3 / 3, qc1 * deltaT^2 / 2,
          qc1 * deltaT^2 / 2,       qc1 * deltaT]

bigQ :: M.Herm Double
bigQ = M.trustSym $ (2 >< 2) bigQl2

bigHH :: M.Matrix Double
bigHH = (1 >< 2) [1, 0]

bigR :: M.Herm Double
bigR = M.trustSym $ (1 >< 1) [sigma1^2]
\end{code}
%endif

I can sample a whole path from this model from which I will reveal a
step at a time.

\pause

%{
%format G.Normal = "{\color{blue} G.Normal}"
%format bigAA    = "{\color{blue} \boldsymbol{A}}"
%format bigHH    = "{\color{blue} \boldsymbol{H}}"
%format bigQ     = "{\color{blue} \boldsymbol{Q}}"
%format bigR     = "{\color{blue} \boldsymbol{R}}"

\begin{code}
xInit :: M.Vector Double
xInit = M.fromList [0, 1]

apolloSamples :: [(M.Vector Double, M.Vector Double)]
apolloSamples = evalState  (ML.unfoldrM apolloSample xInit)
                           (pureMT 17)
  where
    apolloSample xPrev = do
      xNew <- sample $ rvar (G.Normal (bigAA #> xPrev) bigQ)
      yNew <- sample $ rvar (G.Normal (bigHH #> xNew) bigR)
      return $ Just ((xNew, yNew), xNew)
\end{code}
%}

\end{frame}

\begin{frame}{Prior}

Based on your knowledge about Apollo you have a (prior) belief about
the mean of this distribution: it's centred at $(0, 1)$ and both the
position and velocity are probably $\pm 3$.

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
bigSigma0 = M.sym $ (2 >< 2)  [  1.0, 0.0,
                                 0.0, 1.0]

priorsApollo :: Histogram (Bin2D BinD BinD) Double
priorsApollo = hist2 $ conv $
               prePriors 2 100000
  where prePriors seed n =
          evalState  (replicateM n (sample $ G.Normal m0 bigSigma0))
                     (pureMT (fromIntegral seed))
\end{code}

\end{frame}

\begin{frame}{An Aside: Marginals}

\note{Using the \textit{histogram-fill} package, we can easily
  calculate the marginals. Note that these \emph{are} un-normalized.}

\begin{code}
marginalX, marginalY :: Histogram BinD Double
marginalX = H.reduceX H.sum priorsApollo
marginalY = H.reduceY H.sum priorsApollo
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

\begin{block}{Recall the state update}

    $$
    \begin{aligned}
    \boldsymbol{x}_k &= \boldsymbol{A} \boldsymbol{x}_{k-1} + \boldsymbol{Q}_k
    \end{aligned}
    $$
\end{block}

\pause

In Haskell

\begin{code}
newState :: MonadRandom m =>
             M.Vector Double -> m (M.Vector Double)
newState xPrev = sample $ rvar (G.Normal (bigAA #> xPrev) bigQ)
\end{code}

\pause

\begin{block}{Recall the observations are given by}

    $$
    \begin{aligned}
    \boldsymbol{y}_k &= \boldsymbol{H} \boldsymbol{x}_{k}   + \boldsymbol{R}_k
    \end{aligned}
    $$
\end{block}

\pause

Thus the likelihood of observation $b$ given the state $a$ is given by

\begin{code}
likelihoodA :: M.Vector Double -> M.Vector Double -> Double
likelihoodA a b = pdf (G.Normal (bigHH #> a) bigR) b
\end{code}

\end{frame}

\begin{frame}{Apply Bayes'}

\note{Note again \emph{un-normalized.}}

Now you can use Bayes' to track the path of the apollo.

\begin{code}
posteriorizeA :: MonadRandom m =>
                 Histogram (Bin2D BinD BinD) Double ->
                 M.Vector Double ->
                 m (Histogram (Bin2D BinD BinD) Double)
posteriorizeA q d = do
  let xfs = H.asList q
      xs = map pair2Vec $ map fst xfs
      fs = map snd xfs
  xNew <- mapM newState xs
  let newQ = hist2' $ zip (map vec2Pair xNew) fs
  return $
    H.bmap (\(u, v) p -> p * likelihoodA (M.vector [u, v]) d) newQ
\end{code}

%if style == newcode
\begin{code}
vec2Pair :: M.Vector Double -> (Double, Double)
vec2Pair v = (v M.! 0, v M.! 1)

pair2Vec :: (Double, Double) -> M.Vector Double
pair2Vec (x, y) = M.vector [x, y]

test :: Histogram (Bin2D BinD BinD) Double
test = evalState (posteriorizeA priorsApollo (head $
                                              map (M.vector . pure . (!!0) . M.toList) $
                                              map snd
                                              $ take 10 apolloSamples))
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

\note{-0.36189 is the observed position and \textit{reduceX }
 marginalizes out the $x$ component.}

\begin{verbatim}
*Test> take 1 apolloSamples
[([9.99586e-2,0.75767],[-0.36189])]
\end{verbatim}

\begin{code}
marginal1X, marginal1Y :: Histogram BinD Double
marginal1X = H.reduceX H.sum test
marginal1Y = H.reduceY H.sum test

m0X1, m1X1, muX1 :: Double
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
0.94108
*Test> muY1
-0.26897
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
\boldsymbol{x}_i &= \boldsymbol{A}_{i-1}\boldsymbol{x}_{i-1} + \boldsymbol{Q}_{i-1} \\
\boldsymbol{y}_i &= \boldsymbol{H}_i\boldsymbol{x}_i + \boldsymbol{R}_i
\end{aligned}
$$

where

$$
\boldsymbol{Q}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(x)}_i\big)
,\quad
\boldsymbol{R}_i \sim {\cal{N}}\big(0,\boldsymbol{\Sigma}^{(y)}_i\big)
$$

The model is \textbf{linear} and the errors are \textbf{Gaussian}

\end{frame}

\begin{frame}{P\'olya}

\note{There are also some additional assumptions: that the past is
  independent ofthe future given the present. It's possible we could
  derive this assumption from a more general setting but we don't need
  it in any event.}

One way to solve a problem is to generalize it

$$
\begin{aligned}
x_k \sim p(x_k \, \vert \, x_{k-1}) \\
y_k \sim p(y_k \, \vert \, x_{k})
\end{aligned}
$$

\begin{itemize}

\item

State update is Markovian

\item

The past is independent of the future given the present

\item

The measurement given the current state is independent of the past
(quasi-Markovian)

\end{itemize}

\pause

Re-writing our linear Gaussian model

$$
\begin{aligned}
p(x_k \, \vert \, x_{k-1}) &= {\mathcal{N}}(x_k \, \vert \, A_{k-1}x_{k-1}, Q_k) \\
p(y_k \, \vert \, x_k) &= {\mathcal{N}}(y_k \, \vert \, H_{k}x_{k}, R_k)
\end{aligned}
$$
\end{frame}

\begin{frame}{Bayesian Filtering Equations}

\note{This is the prediction step.}

What we'd like is

$$
\begin{aligned}
p(x_k \, \vert \, y_{1:k})
\end{aligned}
$$

$$
\begin{aligned}
p(x_k, x_{k-1} \, \vert \, y_{1:k-1}) &= p(x_k \, \vert \, x_{k-1}, y_{1:k-1})p(x_{k-1} \, \vert \, y_{1:k-1}) \\
&= p(x_k \, \vert \, x_{k-1})p(x_{k-1} \, \vert \, y_{1:k-1})
\end{aligned}
$$

Marginalising
$$
p(x_k \, \vert \, y_{1:k-1}) = \int p(x_k \, \vert \, x_{k-1})p(x_{k-1} \, \vert \, y_{1:k-1}) \mathrm{d}x_{k-1}
$$

\end{frame}

\begin{frame}{A Histogram Interpretation}

\note{We are being sloppy here: the last equality is not an equality
 but only a limit if for each distribution $G = \sum
 w_i\delta(x^{(i)})$ we sample an increasingly large number of samples
 from $p(x_k \, \vert \, x^{(i)}_{k-1})$.}

We can think of our histograms as roughly being a distribution

$$
G = \sum w_i\delta(x^{(i)})
$$

where the $x_i$ are the midpoints of the cells and the $w_i$ are the
(normalised) number of observations in the cell and $\delta$ is the
Dirac delta "function".

$$
\begin{aligned}
p(x_k \, \vert \, y_{1:k-1}) &= \int p(x_k \, \vert \, x_{k-1})p(x_{k-1} \, \vert \, y_{1:k-1}) \mathrm{d}x_{k-1} \\
&= \sum p(x_k \, \vert \, x_{k-1}) w_i\delta(x_{k-1}^{(i)}) \\
&= \sum  w_i \delta(x_{k}^{(i)})
\end{aligned}
$$

where $x_k^{(i)} \sim p(x_k \, \vert \, x_{k-1}^{(i)})$

\end{frame}

\begin{frame}{Bayesian Filtering Equations}

\note{This is the prediction step and the correction step
 combined. Note that we could actually implement this and it's usually
 called particle filtering. There's a bit more to it, for example, we
 could end up with too many samples or almost all the probability
 could end up focussed on one sample; particle filtering works round
 these issues.}

$$
\begin{aligned}
p(x_k \, \vert \, y_{1:k}) &= \frac{1}{Z_k}p(y_k \, \vert\, x_k, y_{1:k-1}) p(x_k \, \vert \, y_{1:k-1}) \\
&=  \frac{1}{Z_k} p(y_k \, \vert\, x_k) p(x_k \, \vert \, y_{1:k-1}) \\
\end{aligned}
$$

In histogram terms
$$
\begin{aligned}
p(x_k^{(i)} \, \vert \, y_{1:k}) &= \frac{1}{Z_k} p(y_k \, \vert\, x_k^{(i-1)}) w_i \delta(x_{k-1}^{(i)}) \\
&= \frac{1}{Z_k} w'_i \delta(x_{k-1}^{(i)})
\end{aligned}
$$

where $w'_i = \frac{1}{Z_k} p(y_k \, \vert\, x_k^{(i-1)}) w_i$ and
$\frac{1}{Z_k} = \sum w'_i$.

\end{frame}

\begin{frame}{Kalman Itself}

\note{Remember the algebraic manipulation in the first example. Well
 there's a lot more of this to derive Kalman filtering and it really
 is very difficult to see the wood for the trees.}

A \textbf{lot} of algebraic manipulation gives the optimal solution.

\begin{block}{Prediction Step}
$$
\begin{aligned}
\hat{\boldsymbol{x}}^\flat_i &=
\boldsymbol{A}_{i-1}\hat{\boldsymbol{x}}_{i-1} \\
\hat{\boldsymbol{\Sigma}}^\flat_i &= \boldsymbol{A}_{i-1}
                                     \hat{\boldsymbol{\Sigma}}_{i-1}
                                     \boldsymbol{A}_{i-1}^\top
                                   + \boldsymbol{Q}_{i-1}
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
\boldsymbol{H}_i^\top + \boldsymbol{R}_i \\
\boldsymbol{K}_i & = \hat{\boldsymbol{\Sigma}}^\flat_i
\boldsymbol{H}_i^\top\boldsymbol{S}^{-1}_i \\
\hat{\boldsymbol{x}}_i &= \hat{\boldsymbol{x}}^\flat_i + \boldsymbol{K}_i\boldsymbol{v}_i \\
\hat{\boldsymbol{\Sigma}}_i &= \hat{\boldsymbol{\Sigma}}^\flat_i - \boldsymbol{K}_i\boldsymbol{S}_i\boldsymbol{K}^\top_i
\end{aligned}
$$
\end{block}

\end{frame}

\begin{frame}{Original Bayes' Theorem Example}

\note{We can apply the Kalman formul\ae to the original example to
 convince ourselves that Kalman is at least valid in this simple
 case.}

\begin{block}{Re-express Original}
$$
\begin{aligned}
x_k &= x_{k-1} + q_{k-1}, \quad q_{k-1} \sim {\mathcal{N}}(0, 0) \\
y_k &= x_k + r_k, \quad r_k \sim {\mathcal{N}}(0, R)
\end{aligned}
$$
\end{block}

\begin{block}{Prediction Step}
$$
\hat{{x}}^\flat_i =\hat{{x}}_{i-1} \quad \hat{{\Sigma}}^\flat_i = \hat{{\Sigma}}_{i-1}
$$
\end{block}

\begin{block}{Correction Step}
$$
\begin{aligned}
{v}_i & = {y}_i - \hat{{x}}^\flat_i \quad {S}_i = \hat{{\Sigma}}^\flat_i + {R}_i \quad
{K}_i = \frac{\hat{{\Sigma}}^\flat_i}{{S}_i} \\
\hat{{x}}_i &= \hat{{x}}^\flat_i + {K}_i{v}_i \quad
\hat{{\Sigma}}_i = \hat{{\Sigma}}^\flat_i - {K}_i{S}_i{K}
\end{aligned}
$$
\end{block}

\end{frame}

\begin{frame}{Original Bayes' Theorem Example}

\begin{block}{Prediction Step}
$$
\begin{aligned}
\hat{{x}}^\flat_i & =\hat{{x}}_{i-1} \quad \hat{{\Sigma}}^\flat_i = \hat{{\Sigma}}_{i-1} \\
\hat{{x}}^\flat_1 & =\hat{{x}}_{0} = \mu_0 \quad \hat{{\Sigma}}^\flat_1 = \hat{{\Sigma}}_{0} = \sigma_0^2
\end{aligned}
$$
\end{block}

\begin{block}{Correction Step}
$$
\begin{aligned}
{v}_i & = {y}_i - \hat{{x}}^\flat_i \quad {S}_i = \hat{{\Sigma}}^\flat_i + {R}_i \quad
{K}_i = \frac{\hat{{\Sigma}}^\flat_i}{{S}_i} \\
\hat{{x}}^i &= \hat{{x}}^\flat_i + {K}_i{v}_i \quad
\sigma_1^2 = \hat{{\Sigma}}^\flat_i - {K}_i{S}_i{K}
\end{aligned}
$$
\end{block}

\begin{block}{Correction Step}
$$
\begin{aligned}
{v} & = {y} - {\mu_0} \quad {S} = \sigma^2_0 + \sigma^2 \quad
{K} = \frac{\sigma^2_0}{\sigma^2_0 + \sigma^2} \\
\mu_1 &= \mu_0 + {K}{v} \quad
\sigma_1^2 = \sigma_0^2 - {K}{S}{K}
\end{aligned}
$$
\end{block}

\end{frame}

\begin{frame}{Original Bayes' Theorem Example}

\begin{block}{Correction Step}
$$
\begin{aligned}
{v} & = {y} - {\mu_0} \quad {S} = \sigma^2_0 + \sigma^2 \quad
{K} = \frac{\sigma^2_0}{\sigma^2_0 + \sigma^2} \\
\mu_1 &= \mu_0 + {K}{v} \quad
\sigma_1^2 = \sigma_0^2 - {K}{S}{K}
\end{aligned}
$$
\end{block}

\begin{block}{Correction Step}
$$
\begin{aligned}
\mu_1 &= \mu_0 + \frac{\sigma^2_0}{\sigma^2_0 + \sigma^2}(y - \mu_0) \quad
\sigma_1^2 = \sigma_0^2 - \frac{(\sigma_0^2)^2}{\sigma_0^2 + \sigma^2}
\end{aligned}
$$
\end{block}

\pause

\begin{block}{Correction Step}
$$
\begin{aligned}
\frac{\mu_1}{\sigma_1^2} &= \frac{\mu_0}{\sigma_0^2} + \frac{\mu}{\sigma^2} \quad
\frac{1}{\sigma_1^2} = \frac{1}{\sigma^2} + \frac{1}{\sigma_0^2}
\end{aligned}
$$
\end{block}

\end{frame}

\section{Some References}

\begin{frame}{Further Reading}

\begin{itemize}
  \item Simo S\"arkk\"a's book "Bayesian Filtering and Smoothing"
  \item Libbi: http://libbi.org
  \item Haskell package \texttt{random-fu} for random variables
  \item Haskell package \texttt{histogram-fill} for histograms
  \item Haskell package \texttt{kalman} for kalman, extended kalman, unscented kalman, particle filtering and smoothing
  \item Same in Python: https://filterpy.readthedocs.io/en/latest
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

\begin{frame}{Speculative Thoughts}

\note{We could write our system as a stochastic differential
 equation. Continuous models abstracts away from unnecessary details
 but the mathematics are harder and beyond a small number of
 dimensions become computationally infeasible. In derivatives prices,
 PDEs are preferred for speed and stability of the greeks but can't
 always be used; Monte Carlo can always be used.}

$$
d \mathbf{X}_{t}=\boldsymbol{\mu}\left(\mathbf{X}_{t}, t\right) d t+\boldsymbol{\sigma}\left(\mathbf{X}_{t}, t\right) d \mathbf{W}_{t}
$$

$$
\begin{aligned}
\frac{\partial}{\partial t} p(t, \mathbf{x})+\sum_{k=1}^{n} \frac{\partial}{\partial x_{k}}\left({\mu}_{k}(t, \mathbf{x}) p(t, \mathbf{x})\right)&= \\ \frac{1}{2} \sum_{j=1, k=1}^{n} \frac{\partial^{2}}{\partial x_{j} \partial x_{k}}\left[\left(\sigma(t, \mathbf{x}) \sigma^{T}(t, \mathbf{x})\right)_{j k} p(t, \mathbf{x})\right]
\end{aligned}
$$

\end{frame}

\end{document}
