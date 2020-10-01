\documentclass{article}
\usepackage{amsmath}
\usepackage{graphicx}
%include polycode.fmt
%options ghci

\title{Filtering Estimation Example}
\author{Dominic Steinitz}
\date{21 September 2020}

\begin{document}

\maketitle

\section{Introduction}

Apparently lynxes feed entirely on snowshoe hares. The Hudson Bay
 Company collected data on numbers of both populations as shown in Figure \ref{fig:hudsonbay_0}.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{diagrams/HudsonBay.png}
\caption{Hudson Bay Data}
\label{fig:hudsonbay_0}
\end{figure}

\section{Mathematical Model}

The Lotka-Volterra Model for a two species predator prey system is

$$
\begin{aligned}
\frac{\mathrm{d}x}{\mathrm{d}t} &= \alpha x - \beta x y \\
\frac{\mathrm{d}y}{\mathrm{d}t} &= \delta x y - \gamma y
\end{aligned}
$$

where $x$ is the number of prey, $y$ is the number of predators,
  $\alpha$ is the prey birth rate, $\beta$ is the prey death rate per
  predator, $\delta$ is the predator birth rate per prey and $\gamma$
  is the predator death rate.

\subsection{Solving Lotka-Volterra}

We can solve the Lotka-Volterra model using the Haskell bindings to
  SUNDIALS. The use of the {\em{Bogacki Shampine}} method is entirely
  arbitrary. The results are shown in Figure \ref{fig:examplelvsolution_0}.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{diagrams/ExampleLvSolution.png}
\caption{Example Model Result}
\label{fig:examplelvsolution_0}
\end{figure}

\section{Estimating the Parameters}

How might we estimate $\alpha, \beta, \delta$ and $\gamma$ from the
 observed data?

One way is to augment the model. The model does take into account that
  hare and lynx birth and death rates might vary from year to
  year. For example, the weather may determine the rates to some
  extent. Rather than try to expand the model to include the weather
  we can capture the fact that we don't know what the parameters are
  and also our feeling that the further we go into the future the less
  certain we are by modelling the parameters as Brownian Motion offset
  by some value:

$$
\mathrm{d}\mathbf{\Phi} = \mathbf{\sigma}\mathrm{d}\mathbf{W}_t
$$

where

$$
\mathbf{\Phi} = \begin{bmatrix} \alpha \\ \beta \\ \delta \\ \gamma \end{bmatrix}
$$

and $\mathbf{W}_t$ is Brownian Motion\footnote{As is well known,
 Brownian Motion is nowhere differentiable and this is a standard way
 of writing what should be an integral equation.}.

We should also accept that the deterministic part of the model may not
 explain everything. Thus the enhanced model is now

$$
\begin{aligned}
\frac{\mathrm{d}x}{\mathrm{d}t} &= \alpha x - \beta x y + \sigma_x\mathrm{d}{W}^x_t\\
\frac{\mathrm{d}y}{\mathrm{d}t} &= \delta x y - \gamma y  + \sigma_y\mathrm{d}{W}^y_t\\
\mathrm{d}\mathbf{\Phi}         &= \mathbf{\sigma}\mathrm{d}\mathbf{W}_t
\end{aligned}
$$

where $W^l_t$ and $W^h_t$ are Brownian Motion.

Of course, we can estimate the parameters in many other ways, each
 with its own, possibly augmented, model.

We still have to explain how to handle observations. The observation
 model is almost trivial:

$$
\begin{aligned}
u_t &= x_t + \sigma_u W^u_t \\
v_t &= y_t + \sigma_v W^v_t
\end{aligned}
$$

where again $W^l_t$ and $W^h_t$ are Brownian Motion, this time
 independent from the state model Brownian Motions.\footnote{There is
 no need for the noise in the state model or the observation model to
 be normal or additive but trying to make everything as general as
 possible will only obscure matters.}.

What we want to estimate is $x_t, y_t, \alpha_t, \beta_t, delta_t$ and
 $\gamma_t$ given $u_t$ and $v_t$, the numbers of hares and lynxes
 given at times $t_0, t_1, \ldots, t_N$.

Rather than write down the mathematical theory for this which involves
       a lot of technical machinery, the basic idea is to use the
       state equations to move the state forward in time and then
       apply a Bayesian update step to the prior distribution just
       before an observation to produce a posterior distribution just
       after the observation. Note that the state is actually a
       distribution. One way of approximating this distribution is to
       use a set of samples. This is generally known as particle
       filtering. We could also represent the distribution using some
       well known distribution e.g. the Normal distribution. This is
       known as Kalman filtering\footnote{Strictly speaking, we would
       have to use extended or unscented Kalman filtering as vanilla
       Kalman requires the state update is linear --- clearly not the
       case here.}.

%if False
\begin{code}
{-# LANGUAGE TypeFamilies      #-}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE TypeOperators     #-}
{-# LANGUAGE NumDecimals       #-}
{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes       #-}
{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE BangPatterns        #-}
{-# LANGUAGE DeriveFunctor        #-}
{-# LANGUAGE DeriveFoldable        #-}

{-# OPTIONS_GHC -Wall          #-}

module LotkaVolterra (main) where

import Prelude as P

import           Numeric.LinearAlgebra
import           Numeric.Sundials

import           Control.Exception
import           Katip
import           Katip.Monadic
import           System.IO (stderr)
import           GHC.Int

import           Foreign.C.Types
import           Data.Coerce

import           Control.Monad.Writer

import qualified Data.Vector as V

import           Data.List (transpose)

import qualified Language.R as R
import           Language.R.QQ

import           Frames
import           Control.Lens((^.))
import           System.IO.Unsafe (unsafePerformIO)
import qualified Data.Foldable as F

import           Numeric.Particle
import           Data.Random.Distribution.MultivariateNormal ( Normal(..) )
import qualified Data.Random as R
\end{code}
%endif

%if False
\begin{code}
tableTypes "NoisyObs" "lynx_hare_df.csv"

loadNoisyObs :: IO (Frame NoisyObs)
loadNoisyObs = inCoreAoS (readTable "lynx_hare_df.csv")

predPreyObs :: [((Int, Double), (Int, Double))]
predPreyObs = unsafePerformIO $
          do xs <- loadNoisyObs
             let ys = F.toList $ fmap (\u -> ((u ^. year), (u ^. hare))) xs
             let zs = F.toList $ fmap (\u -> ((u ^. year), (u ^. lynx))) xs
             return $ zip ys zs

predPreyObs' :: V.Vector (SystemObs Double)
predPreyObs' = unsafePerformIO $
          do xs <- loadNoisyObs
             let ys = V.fromList $ F.toList $ fmap (^. hare) xs
                 zs = V.fromList $ F.toList $ fmap (^. lynx) xs
             return $ V.zipWith SystemObs ys zs
\end{code}
%endif


%if False
\begin{code}
data Rate a = Rate { theta1 :: a, theta2 :: a, theta3 :: a, theta4 :: a }
  deriving Show

meanRate :: Floating a => Rate a
meanRate = Rate { theta1 = 0.5, theta2 = 0.025, theta3 = 0.8, theta4 = 0.025 }

dzdt :: (Show a, Floating a) => Rate a -> a -> [a] -> [a]
dzdt x _t [u, v] = [ (alpha - beta * v) * u
                   , (-gamma + delta * u) * v
                   ]
  where
    alpha = theta1 x
    beta = theta2 x
    delta = theta4 x
    gamma = theta3 x
dzdt _x _t xs = error $ show xs

lotkaVolterra :: Double -> Double -> OdeProblem
lotkaVolterra h l = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t x -> fromList (dzdt meanRate t (toList x))
  , odeJacobian = Nothing
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = vector [h, l]
  , odeSolTimes = vector us
  , odeTolerances = defaultTolerances
  }

us :: [Double]
us = map (* 0.05) $ map fromIntegral ([0 .. 400] :: [Int])

vs :: [Double]
vs = map (+ 1900) us

sol :: Double -> Double -> IO (Matrix Double)
sol h l = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra h l)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

main :: IO ()
main = do
  let hs, ks, ts :: [Double]
      ts = map fromIntegral $ map fst $ map fst predPreyObs
      hs = map snd $ map fst predPreyObs
      ks = map snd $ map snd predPreyObs
  y <- sol (hs!!0) (ks!!0)
  let rs = transpose $ toLists y
  (as1, bs1, cs1) <- testL
  let preDs1 = V.toList $ V.map V.toList $ V.map (V.map (exp . hares1)) cs1
      ds1 = concat $ zipWith (\t ws -> zip (repeat t) ws) ts preDs1
      es1 = map fst ds1
      fs1 = map snd ds1
      gammas1 = V.toList $
                V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map gamma1 lls) cs1

  R.runRegion $ do
    _ <- [r| library(ggplot2) |]
    c0 <- [r| ggplot() |]
    c1 <- [r| c0_hs + ggtitle("Hares and Lynxes") |]
    c2 <- [r| c1_hs + xlab("Year") |]
    c3 <- [r| c2_hs + ylab("Animals (000s)") |]
    c4 <- [r| c3_hs + labs(colour = "Species") |]
    c5 <- [r| c4_hs + theme(plot.title = element_text(hjust = 0.5)) |]
    _  <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 ([ ("Hares", hs, ts), ("Lynxes", ks, ts) ] :: [(String, [Double], [Double])])
    _ <-  [r| ggsave(filename="diagrams/HudsonBay.png") |]

    e1 <- [r| c0_hs + ggtitle("Hares and Lynxes") |]
    e2 <- [r| e1_hs + xlab("Year") |]
    e3 <- [r| e2_hs + ylab("Animals (000s)") |]
    e4 <- [r| e3_hs + labs(colour = "Species") |]
    e5 <- [r| e4_hs + theme(plot.title = element_text(hjust = 0.5)) |]
    e6 <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                e5 ([ ("Hares", hs, ts), ("Lynxes", ks, ts)
                    , ("HaresL", (V.toList as1), ts), ("LynxesL", (V.toList bs1), ts) ] :: [(String, [Double], [Double])])
    _ <-  [r| e6_hs + geom_point(aes(x=es1_hs, y=fs1_hs)) |]
    _ <-  [r| ggsave(filename="diagrams/HudsonBayFit.png") |]

    d1 <- [r| c0_hs + ggtitle("Hares and Lynxes") |]
    d2 <- [r| d1_hs + xlab("Year") |]
    d3 <- [r| d2_hs + ylab("Animals (000s)") |]
    d4 <- [r| d3_hs + labs(colour = "Species") |]
    d5 <- [r| d4_hs + theme(plot.title = element_text(hjust = 0.5)) |]
    _  <- foldM (\_ (n, rr, tt) -> do
                        [r| d5_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 ([("Gamma", gammas1, ts)] :: [(String, [Double], [Double])])
    _  <-  [r| ggsave(filename="diagrams/Gamma.png") |]

    _  <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 ([("Predicted Hares", rs!!0, vs), ("Predicted Lynxes", rs!!1, vs)] :: [(String, [Double], [Double])])
    _ <- [r| ggsave(filename="diagrams/ExampleLvSolution.png") |]

    return ()

defaultOpts :: method -> ODEOpts method
defaultOpts method = ODEOpts
  { maxNumSteps = 1e5
  , minStep     = 1.0e-14
  , fixedStep   = 0.0
  , maxFail     = 10
  , odeMethod   = method
  , initStep    = Nothing
  , jacobianRepr = DenseJacobian
  }

instance Element Int8

emptyOdeProblem :: OdeProblem
emptyOdeProblem = OdeProblem
      { odeRhs = error "emptyOdeProblem: no odeRhs provided"
      , odeJacobian = Nothing
      , odeInitCond = error "emptyOdeProblem: no odeInitCond provided"
      , odeEventDirections = V.empty
      , odeEventConditions = EventConditionsHaskell V.empty
      , odeTimeBasedEvents = TimeEventSpec $ return $ 1.0 / 0.0
      , odeEventHandler = nilEventHandler
      , odeMaxEvents = 100
      , odeSolTimes = error "emptyOdeProblem: no odeSolTimes provided"
      , odeTolerances = defaultTolerances
      }

nilEventHandler :: EventHandler
nilEventHandler _ _ _ = throwIO $ ErrorCall "nilEventHandler"

defaultTolerances :: Tolerances
defaultTolerances = Tolerances
  { absTolerances = Left 1.0e-5
  , relTolerance = 1.0e-10
  }
\end{code}
%endif

\subsection{Sampling Potential Paths}

We can simulate a whole set of trajectories by sampling from
 ${\mathcal{N}}(\mu, P)$ where $\mu = [\mu_h, \mu_l, m_\gamma]$ and

$$
P =
\begin{bmatrix}
\sigma^2_h &        0.0 &             0.0 \\
0.0 &        \sigma^2_l &             0.0 \\
0.0 &               0.0 &  \sigma^2_\gamma
\end{bmatrix}
$$

A set of \eval{nParticles} runs is shown in Figure~\ref{fig:lvparticles_0} with
$[\mu_h, \mu_l, m_\gamma] = \eval{m0}$.

$$
P =
\begin{bmatrix}
\eval{(unSym bigP)!0!0} &                     0.0 &                     0.0 \\
0.0                     & \eval{(unSym bigP)!1!1} &                     0.0 \\
0.0                     &                     0.0 & \eval{(unSym bigP)!2!2}
\end{bmatrix}
$$

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{diagrams/LvParticles.png}
\caption{Particle Paths}
\label{fig:lvparticles_0}
\end{figure}

%if False
\begin{code}
data SystemState1 a = SystemState1 { hares1  :: a, lynxes1  :: a, alpha1 :: a, beta1 :: a, delta1 :: a, gamma1 :: a}
  deriving (Show, Functor, Foldable)

data SystemObs a = SystemObs { obsHares  :: a, obsLynxes :: a}
  deriving Show

nParticles :: Int
nParticles = 1000

measureOpL :: Particles (SystemState1 Double) -> Particles (SystemObs Double)
measureOpL = V.map (\s -> SystemObs { obsHares = exp $ hares1 s, obsLynxes = exp $ lynxes1 s})

weight :: SystemObs Double -> SystemObs Double -> Double
weight obs predicted = R.pdf (Normal xs bigR) ys
  where
    xs = vector [obsHares obs, obsLynxes obs]
    ys = vector [obsHares predicted, obsLynxes predicted]

bigR :: Herm Double
bigR = sym $ (2><2) [ 32.0e-1, 0.0,
                      0.0,    32.0e-1
                    ]

scanMapM :: Monad m => (s -> a -> m s) -> (s -> m b) -> s -> V.Vector a -> m (V.Vector b)
scanMapM f g !s0 !xs
  | V.null xs = do
    u <- g s0
    return $ V.singleton u
  | otherwise = do
    s <- f s0 (V.head xs)
    u <- g s0
    liftM (u `V.cons`) (scanMapM f g s (V.tail xs))

lotkaVolterra3 :: [Double] -> [Double] -> OdeProblem
lotkaVolterra3 ts s = -- trace (show r ++ " " ++ show y) $
  emptyOdeProblem
  { odeRhs = odeRhsPure $ \t x -> fromList (dzdt (Rate (u!!0) (u!!1) (u!!2) (u!!3)) t (toList x))
  , odeJacobian = Nothing
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = vector y
  , odeSolTimes = vector ts
  , odeTolerances = defaultTolerances
  }
  where
    u = coerce $ take 4 $ drop 2 s
    y = take 2 $ s

sol3 :: MonadIO m => [Double] -> [Double] -> m (Matrix Double)
sol3 ts s = do
  -- handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  -- logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runNoLoggingT $ solve (defaultOpts BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra3 ts s)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

m0L :: [Double]
m0L = fmap log [0.5, 0.025, 0.8, 0.025, 30, 4]

initParticlesL :: R.MonadRandom m =>
                 m (Particles (SystemState1 Double))
initParticlesL = V.replicateM nParticles $ do
  rr <- R.sample $ R.rvar (Normal (vector m0L) (sym $ (6><6) bigQ))
  return $ fmap exp $ SystemState1 { alpha1  = rr!0
                                   , beta1   = rr!1
                                   , delta1  = rr!2
                                   , gamma1  = rr!3
                                   , hares1  = rr!4
                                   , lynxes1 = rr!5
                                   }

stateUpdateL :: Particles (SystemState1 Double) -> IO (Particles (SystemState1 Double))
stateUpdateL ps = do
  qs <-  V.mapM (sol3 (take 21 us)) (V.map (F.toList . fmap exp) ps)
  rr <- V.replicateM nParticles $
        R.sample $ R.rvar (Normal (vector (replicate 6 0.0)) (sym $ (6><6) bigQ))

  let ms :: V.Vector (Vector Double)
      ms = V.map (log . (!20)) qs
      ns = V.zipWith3 (\p q u -> SystemState1 { alpha1  = u!0 + alpha1 q
                                              , beta1   = u!1 + beta1 q
                                              , delta1  = u!2 + delta1 q
                                              , gamma1  = u!3 + gamma1 q
                                              , hares1  = u!4 + p!0
                                              , lynxes1 = u!5 + p!1
                                              })
                      ms ps rr
  return ns

testL :: IO (V.Vector Double, V.Vector Double, V.Vector (Particles (SystemState1 Double)))
testL = do
  is <- initParticlesL
  foo <- scanMapM (runPF stateUpdateL measureOpL weight) return (V.map (fmap log) is) (V.drop 1 predPreyObs')
  let as = V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map hares1 lls) foo
  let bs = V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map lynxes1 lls) foo
  return (as, bs, foo)

bigQ :: [Double]
bigQ = [ 1.0e-2, 0.0,    0.0,    0.0,    0.0,    0.0
       , 0.0,    5.0e-3, 0.0,    0.0,    0.0,    0.0
       , 0.0,    0.0,    1.0e-2, 0.0,    0.0,    0.0
       , 0.0,    0.0,    0.0,    5.0e-3, 0.0,    0.0
       , 0.0,    0.0,    0.0,    0.0,    1.0e-1, 0.0
       , 0.0,    0.0,    0.0,    0.0,    0.0,    1.0e-1
       ]
\end{code}
%endif

model {
  theta[{1, 3}] ~ normal(1, 0.5);
  theta[{2, 4}] ~ normal(0.05, 0.05);
  sigma ~ lognormal(-1, 1);
  z_init ~ lognormal(log(10), 1);
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]);
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  }
}

meanRate :: Floating a => Rate a
meanRate = Rate { theta1 = 0.5, theta2 = 0.025, theta3 = 0.8, theta4 = 0.025 }

m0 :: Vector Double
m0 = vector [30.0, 4.0, 2.5e-2]

\end{document}
