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

Let us assume that $\alpha, \beta, \delta$ are given but that we wish
 to estimate $\gamma$ from the observed data.

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

{-# OPTIONS_GHC -Wall          #-}

module LotkaVolterra (main) where

import Prelude as P

import           Numeric.LinearAlgebra
import           Numeric.Sundials

import           Control.Exception
import           Katip
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
import           Control.Monad.State ( evalState, replicateM )
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

\end{code}
%endif

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
  is <- initParticles
  ys :: V.Vector (Matrix Double) <-  V.mapM sol' is
  let ss = V.map (transpose . toLists) ys
      tt :: V.Vector [(String, [Double], [Double])]
      tt = V.map (\ss -> [("Predicted Hares",  ss!!0, vs),
                          ("Predicted Lynxes", ss!!1, vs)]) ss
      uu = concat $ V.toList tt
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
                c5 ([("Hares", hs, ts), ("Lynxes", ks, ts)] :: [(String, [Double], [Double])])
    _ <-  [r| ggsave(filename="diagrams/HudsonBay.png") |]

    _  <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 ([("Predicted Hares", rs!!0, vs), ("Predicted Lynxes", rs!!1, vs)] :: [(String, [Double], [Double])])
    _ <- [r| ggsave(filename="diagrams/ExampleLvSolution.png") |]

    _  <- foldM (\c (n, rr, tt) -> do
                [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 uu
    _ <- [r| ggsave(filename="diagrams/LvParticles.png") |]
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

A set of \eval{nParticles} runs is shown in Figure~\ref{fig:lvparticles_0}.

\begin{figure}[h]
\centering
\includegraphics[width=0.8\textwidth]{diagrams/LvParticles.png}
\caption{Particle Paths}
\label{fig:lvparticles_0}
\end{figure}

%if False
\begin{code}
data SystemState a = SystemState { hares  :: a, lynxes  :: a, gamma :: a}
  deriving Show

data SystemObs a = SystemObs { obsHares  :: a, obsLynxes :: a}
  deriving Show

nParticles :: Int
nParticles = 100

m0 :: Vector Double
m0 = vector [30.0, 4.0, 2.5e-2]

bigP :: Herm Double
bigP = sym $ (3><3) [ 1.0e-1, 0.0,    0.0,
                      0.0,    1.0e-1, 0.0,
                      0.0,    0.0,    1.0e-5
                    ]

initParticles :: R.MonadRandom m =>
                 m (Particles (SystemState Double))
initParticles = V.replicateM nParticles $ do
  rr <- R.sample $ R.rvar (Normal m0 bigP)
  let h = rr!0
      l = rr!1
      g = rr!2
  return $ SystemState { hares = h, lynxes = l, gamma = g}

lotkaVolterra' :: SystemState Double -> OdeProblem
lotkaVolterra' s = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t x -> fromList (dzdt (meanRate {theta4 = coerce $ gamma s}) t (toList x))
  , odeJacobian = Nothing
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = vector [hares s, lynxes s]
  , odeSolTimes = vector us
  , odeTolerances = defaultTolerances
  }

sol' :: SystemState Double -> IO (Matrix Double)
sol' s = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra' s)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)
\end{code}
%endif

\end{document}
