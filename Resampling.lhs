\documentclass[12pt]{book}
%include polycode.fmt
%options ghci -fglasgow-exts

\setlength{\parskip}{\medskipamount}
\setlength{\parindent}{0pt}

\long\def\ignore#1{}

\begin{document}

We previously considered the model for population growth, the logistic equation

\begin{eqnarray}
\dot{p} & =  & rp\Big(1 - \frac{p}{k}\Big)
\end{eqnarray}

We take the parameter $k$ to be fixed and known and $r_0$ to be
sampled from some prior

\begin{eqnarray}
R_0 & \sim & {\cal{N}}(\mu_r, \sigma_r^2) \label{eq:r_prior}
\end{eqnarray}

As before we know the solution to the logistic equation

\begin{eqnarray}
p = \frac{kp_0\exp rt}{k + p_0(\exp rt - 1)}
\end{eqnarray}

where $p_0$ is the size of the population at time $t=0$.

We observe a noisy value of population at regular time intervals

\begin{eqnarray}
x_i &=& \frac{kp_0\exp r\Delta T i}{k + p_0(\exp r\Delta T i - 1)} \\
y_i &=& x_i + \upsilon_i
\end{eqnarray}

\ignore{
\begin{code}
{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
{-# OPTIONS_GHC -fno-warn-missing-methods #-}
{-# OPTIONS_GHC -fno-warn-orphans         #-}
\end{code}
}

\ignore{
\begin{code}
{-# LANGUAGE BangPatterns              #-}
\end{code}
}

We represent the prior at time $i-1$ by

\begin{eqnarray}
\hat{\Phi}_{i-1}(r) = \frac{1}{Z_{i-1}}\sum_{k=1}^N w^{(k)} {\cal{N}}(r; \mu^{(k)}_r, \sigma)
\end{eqnarray}

where the $N$ samples $\mu^{(k)}_r$ are independent samples of the
prior density $\Phi_{i-1}(r)$ and $Z_{i-1}$ is a normalizing constant.

To sample from this Kernel Density Estimate we can pick $j$ with
probability $\frac{w^{(j)}}{Z_{i-1}}$ and then sample
from ${\cal{N}}(\mu^{(j)}_r, \sigma)$.


\begin{code}
module Resampling where

import qualified Data.Vector.Unboxed as U

import Data.Bits

import Control.Monad.ST
import Control.Monad.Loops
import Control.Monad.State
import qualified Control.Monad.Writer as W

import Data.Random hiding ( gamma )
import Data.Random.Source.PureMT
import System.Random.MWC

k :: Double
k = 1.0

muPriorR, sigmaPriorR :: Double
muPriorR = 5.0
sigmaPriorR = 1e1

logit :: Double -> Double -> Double -> Double -> Double
logit p0 k r t = k * p0 * (exp (r * t)) / (k + p0 * (exp (r * t) - 1))

p0 :: Double
p0 = 0.1 * k

deltaT :: Double
deltaT = 0.001

sigma :: Double
sigma = 1e-1

initPopSir :: Int ->
              RVarT (W.Writer [(Double, U.Vector Double)])
                    (Double, U.Vector Double)

initPopSir n = do
  rs <- U.replicateM n (rvarT $ Normal muPriorR sigmaPriorR)
  return (p0, rs)
\end{code}

We create a function to perform one step of the Bayesian update.

\begin{code}
popSir :: Int -> (Double, U.Vector Double) -> Double ->
          RVarT (W.Writer [(Double, U.Vector Double)])
                (Double, U.Vector Double)
popSir n (p0Prev, rsPrev) z = do

  let ps = U.map (\r -> logit p0Prev k r deltaT) rsPrev

  let rBar = U.sum rsPrev / (fromIntegral n)
      pBar = U.sum ps / (fromIntegral n)

  lift $ W.tell [(rBar, U.take 100 rsPrev)]

  let wsRaw = U.map (\p -> pdf (Normal p sigma) z) ps
      sumWs = U.sum wsRaw
      ws    = U.map (/ sumWs) wsRaw

  let vs  :: U.Vector Double
      vs = runST (create >>= (asGenST $ \gen -> uniformVector gen n))

      cumSumWs = U.scanl (+) 0 ws

      js :: U.Vector Int
      js = myIndices (U.tail cumSumWs) vs

      rsTilde = U.map (rsPrev U.!) js

  rsNew <- U.mapM (\mu -> rvarT (Normal mu sigma)) rsTilde

  return (pBar, rsNew)
\end{code}

\begin{code}
binarySearch :: (U.Unbox a, Ord a) =>
                a -> U.Vector a -> Int
binarySearch x vec = loop 0 (U.length vec - 1)
  where
    loop !l !u
      | u <= l    = l
      | otherwise = let e = vec U.! k in if x <= e then loop l k else loop (k+1) u
      where k = (u + l) `shiftR` 1

myIndices :: U.Vector Double -> U.Vector Double -> U.Vector Int
myIndices bs xs = U.map (flip binarySearch bs) xs
\end{code}


\begin{code}
nParticles :: Int
nParticles = 10000

nObs :: Int
nObs = 1000

obsVariance :: Double
obsVariance = 1e-2

r :: Double
r = 6.0

singleSample :: Double -> RVarT (W.Writer [Double]) Double
singleSample p0 = do
  upsilon <- rvarT (Normal 0.0 obsVariance)
  let p1 = logit p0 k r deltaT
  lift $ W.tell [p1 + upsilon]
  return p1

streamSample :: RVarT (W.Writer [Double]) Double
streamSample = iterateM_ singleSample p0

samples :: [Double]
samples = take nObs $ snd $ W.runWriter (evalStateT (sample streamSample) (pureMT 2))

testPopSir :: RVarT (W.Writer [(Double, U.Vector Double)]) (Double, U.Vector Double)
testPopSir = do
  (p0, rsInit) <- initPopSir nParticles
  foldM (popSir nParticles) (p0, rsInit) samples

runPopSir :: ((Double, U.Vector Double), [(Double, U.Vector Double)])
runPopSir = W.runWriter (evalStateT (sample testPopSir) (pureMT 7))
\end{code}

\end{document}