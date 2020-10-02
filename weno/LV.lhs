\documentclass{article}
\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{hyperref}
\setlength{\parskip}{\medskipamount}
\setlength{\parindent}{0pt}
%include polycode.fmt
%options ghci

%format alpha = "\alpha"
%format beta  = "\beta"
%format delta = "\delta"
%format gamma = "\gamma"
%format bigQ  = "Q"
%format m0    = "m_0"

\title{Filtering Estimation Example}
\author{Dominic Steinitz}
\date{21 September 2020}

\begin{document}

\maketitle

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
{-# LANGUAGE DeriveFunctor       #-}
{-# LANGUAGE DeriveFoldable      #-}

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

\section{Using a Filtering Library}

We can use the Haskell library
\href{https://hackage.haskell.org/package/kalman}{kalman} which contains a function:

\eval{:t runPF}

which we can read as a function which takes a state update function
 |Particles a -> m (Particles a)|, an observation function |Particles a -> Particles b|,
 a measurement function for the observations |b -> b -> Double|
 and returns a function |Particles a -> b -> m (Particles a)|
 which takes a distribution (the prior) and an observation and
 returns a new distribution (the posterior)\footnote{The symbol |m|
 can be ignored --- it restricts how |runPF| is used to prevent
 runtime errors.}.

Here's the state update notated slightly differently

\begin{gather}
 \begin{bmatrix}
    \dot{\alpha} \\
    \dot{\beta}  \\
    \dot{\delta} \\
    \dot{\gamma} \\
    \dot{x} \\
    \dot{y}
\end{bmatrix}
 =
\begin{bmatrix}
    0 \\
    0 \\
    0 \\
    0 \\
    \alpha x - \beta x y \\
    \delta x y - \gamma y
\end{bmatrix}
+
\mathbf{Q}\mathbf{W}_t
\end{gather}

where $\mathbf{W}$ is Brownian Motion and $\mathbf{Q}$ is a covariance matrix.

There is no reason to expect any correlations so informed by the
  results in
  \href{https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html}{Predator-Prey
  Population Dynamics: the Lotka-Volterra model in Stan} let us take
  the state update covariance to be

$$
\mathbf{Q} =
\begin{bmatrix}
\eval{bigQ!!0}  & \eval{bigQ!!1}  & \eval{bigQ!!2}  & \eval{bigQ!!3}  & \eval{bigQ!!4}  & \eval{bigQ!!5} \\
\eval{bigQ!!6}  & \eval{bigQ!!7}  & \eval{bigQ!!8}  & \eval{bigQ!!9}  & \eval{bigQ!!10} & \eval{bigQ!!11} \\
\eval{bigQ!!12} & \eval{bigQ!!13} & \eval{bigQ!!14} & \eval{bigQ!!15} & \eval{bigQ!!16} & \eval{bigQ!!17} \\
\eval{bigQ!!18} & \eval{bigQ!!19} & \eval{bigQ!!20} & \eval{bigQ!!21} & \eval{bigQ!!22} & \eval{bigQ!!23} \\
\eval{bigQ!!24} & \eval{bigQ!!25} & \eval{bigQ!!26} & \eval{bigQ!!27} & \eval{bigQ!!28} & \eval{bigQ!!29} \\
\eval{bigQ!!30} & \eval{bigQ!!31} & \eval{bigQ!!32} & \eval{bigQ!!33} & \eval{bigQ!!34} & \eval{bigQ!!35}
\end{bmatrix}
$$

We take the observation noise to be given by

$$
\mathbf{R} =
\begin{bmatrix}
\eval{(unSym bigR)!0!0} & \eval{(unSym bigR)!0!1} \\
\eval{(unSym bigR)!1!0} & \eval{(unSym bigR)!1!1}
\end{bmatrix}
$$

Note that in order to avoid values being negative we work in log
  space. Thus in the state update function below we first exponentiate
  the values before using a solver to move the state forwards in time
  by one time step and the log these new values before adding the
  state noise.

Note also that |alpha, beta, delta| and |gamma| are kept constant
  apart from the addition of the noise; the solver only updates the
  numbers of hares and lynxes according to the Lotka-Volterra
  equation.

\begin{code}
stateUpdate :: Particles (SystemState Double) -> IO (Particles (SystemState Double))
stateUpdate ps = do

  let timeSteps = [0.0, 1.0]
  qs <-  V.mapM (sol' timeSteps) (V.map (F.toList . fmap exp) ps)
  let rs = V.map (log . (!1)) qs

  eps <- V.replicateM nParticles $
         R.sample $ R.rvar (Normal (vector (replicate 6 0.0)) (sym $ (6><6) bigQ))

  let f x y e =
        SystemState { alpha  = e!0 + alpha y,
                      beta   = e!1 + beta  y,
                      delta  = e!2 + delta y,
                      gamma  = e!3 + gamma y,
                      hares  = e!4 + x!0,
                      lynxes = e!5 + x!1
                    }
  return $ V.zipWith3 f rs ps eps
\end{code}

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{diagrams/HudsonBayFit.png}
\caption{Fitting the State}
\label{fig:fitting_state_0}
\end{figure}

\subsection{Initial Distributions}

One thing that is still missing from the model are the initial
 distributions\footnote{This would have been obvious if we had used
 the integral formulation of what are now stochastic differential
 equations.}.

We can create initial distributions from a (multivariate) normal
 distribution with the following mean vector and covariance matrix:

\begin{code}
m0 :: [Double]
m0 = fmap log [0.5, 0.025, 0.8, 0.025, 30, 4]

initParticles :: R.MonadRandom m =>
                 m (Particles (SystemState Double))
initParticles = V.replicateM nParticles $ do
  x <- R.sample $ R.rvar (Normal (vector m0) (sym $ (6><6) bigQ))
  return $ fmap exp $
    SystemState { alpha  = x!0, beta   = x!1, delta  = x!2, gamma  = x!3, hares  = x!4, lynxes = x!5 }
\end{code}

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
  , odeSolTimes = vector uss
  , odeTolerances = defaultTolerances
  }

uss :: [Double]
uss = map (* 0.05) $ map fromIntegral ([0 .. 400] :: [Int])

vs :: [Double]
vs = map (+ 1900) uss

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
  let preDs1 = V.toList $ V.map V.toList $ V.map (V.map (exp . hares)) cs1
      ds1 = concat $ zipWith (\t ws -> zip (repeat t) ws) ts preDs1
      es1 = map fst ds1
      fs1 = map snd ds1
      gammas1 = V.toList $
                V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map gamma lls) cs1

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


%if False
\begin{code}
data SystemState a = SystemState { hares  :: a, lynxes  :: a, alpha :: a, beta :: a, delta :: a, gamma :: a}
  deriving (Show, Functor, Foldable)

data SystemObs a = SystemObs { obsHares  :: a, obsLynxes :: a}
  deriving Show

nParticles :: Int
nParticles = 1000

measureOpL :: Particles (SystemState Double) -> Particles (SystemObs Double)
measureOpL = V.map (\s -> SystemObs { obsHares = exp $ hares s, obsLynxes = exp $ lynxes s})

weight :: SystemObs Double -> SystemObs Double -> Double
weight obs predicted = R.pdf (Normal xs bigR) ys
  where
    xs = vector [obsHares obs, obsLynxes obs]
    ys = vector [obsHares predicted, obsLynxes predicted]

bigR :: Herm Double
bigR = sym $ (2><2) [ 10.0e-1, 0.0,
                      0.0,    10.0e-1
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

lotkaVolterra' :: [Double] -> [Double] -> OdeProblem
lotkaVolterra' ts s =
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

sol' :: MonadIO m => [Double] -> [Double] -> m (Matrix Double)
sol' ts s = do
  w <- runNoLoggingT $ solve (defaultOpts BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra' ts s)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

testL :: IO (V.Vector Double, V.Vector Double, V.Vector (Particles (SystemState Double)))
testL = do
  is <- initParticles
  foo <- scanMapM (runPF stateUpdate measureOpL weight) return (V.map (fmap log) is) (V.drop 1 predPreyObs')
  let as = V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map hares lls) foo
  let bs = V.map (\lls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map exp $ V.map lynxes lls) foo
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

\end{document}
