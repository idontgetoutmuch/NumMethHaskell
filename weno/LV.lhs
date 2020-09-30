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
import qualified Data.Random.Distribution.Normal as RN
import qualified Data.Random as R

import Debug.Trace
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
  is <- initParticles
  foo <- scanMapM (runPF stateUpdate measureOp weight) return is (V.tail predPreyObs')
  let as = take 21 $ V.toList $
           V.map (\ls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map hares ls) foo
      bs = take 21 $ V.toList $
           V.map (\ls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map lynxes ls) foo
      cs = take 21 $ V.toList $
           V.map (\ls -> (* (1 / (fromIntegral nParticles))) $ sum $ V.map gamma ls) foo
      preDs = V.toList $ V.map V.toList $ V.map (V.map hares) foo
      ds = concat $ zipWith (\t vs -> zip (repeat t) vs) ts preDs
      es = map fst ds
      fs = map snd ds
  (as1, bs1, cs1) <- testL
  let preDs1 = V.toList $ V.map V.toList $ V.map (V.map (exp . hares1)) cs1
      ds1 = concat $ zipWith (\t vs -> zip (repeat t) vs) ts preDs1
      es1 = map fst ds1
      fs1 = map snd ds1
  R.runRegion $ do
    _ <- [r| library(ggplot2) |]
    c0 <- [r| ggplot() |]
    c1 <- [r| c0_hs + ggtitle("Hares and Lynxes") |]
    c2 <- [r| c1_hs + xlab("Year") |]
    c3 <- [r| c2_hs + ylab("Animals (000s)") |]
    c4 <- [r| c3_hs + labs(colour = "Species") |]
    c5 <- [r| c4_hs + theme(plot.title = element_text(hjust = 0.5)) |]
    c6 <- foldM (\c (n, rr, tt) -> do
                        [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
                c5 ([ ("Hares", hs, ts), ("Lynxes", ks, ts)
                    , ("HaresL", (V.toList as1), ts), ("LynxesL", (V.toList bs1), ts) ] :: [(String, [Double], [Double])])
    _ <-  [r| c6_hs + geom_point(aes(x=es1_hs, y=fs1_hs)) |]
    _ <-  [r| ggsave(filename="diagrams/HudsonBay.png") |]

    -- c6  <- foldM (\c (n, rr, tt) -> do
    --                     [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
    --             c5 ([("Hares Obs", hs, ts), ("Lynxes Obs", ks, ts),
    --                  ("Hares Pred", as, ts), ("Lynxes Pred", bs, ts)] :: [(String, [Double], [Double])])
    -- _ <-  [r| c6_hs + geom_point(aes(x=es_hs, y=fs_hs)) |]
    -- _ <-  [r| ggsave(filename="diagrams/Posterior.png") |]

    -- d1 <- [r| c0_hs + ggtitle("Hares and Lynxes") |]
    -- d2 <- [r| d1_hs + xlab("Year") |]
    -- d3 <- [r| d2_hs + ylab("Animals (000s)") |]
    -- d4 <- [r| d3_hs + labs(colour = "Species") |]
    -- d5 <- [r| d4_hs + theme(plot.title = element_text(hjust = 0.5)) |]
    -- _  <- foldM (\c (n, rr, tt) -> do
    --                     [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
    --             c5 ([("Gamma", cs, ts)] :: [(String, [Double], [Double])])
    -- _  <-  [r| ggsave(filename="diagrams/Gamma.png") |]

    -- _  <- foldM (\c (n, rr, tt) -> do
    --                     [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
    --             c5 ([("Predicted Hares", rs!!0, vs), ("Predicted Lynxes", rs!!1, vs)] :: [(String, [Double], [Double])])
    -- _ <- [r| ggsave(filename="diagrams/ExampleLvSolution.png") |]

    -- _  <- foldM (\c (n, rr, tt) -> do
    --             [r| c_hs + geom_line(aes(x = tt_hs, y = rr_hs, colour = n_hs)) |])
    --             c5 uu
    -- _ <- [r| ggsave(filename="diagrams/LvParticles.png") |]
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
data SystemState a = SystemState { hares  :: a, lynxes  :: a, gamma :: a}
  deriving Show

data SystemState1 a = SystemState1 { hares1  :: a, lynxes1  :: a, alpha1 :: a, beta1 :: a, delta1 :: a, gamma1 :: a}
  deriving (Show, Functor, Foldable)

data SystemObs a = SystemObs { obsHares  :: a, obsLynxes :: a}
  deriving Show

nParticles :: Int
nParticles = 1000

m0 :: Vector Double
m0 = vector [30.0, 4.0, 2.5e-2]

bigP :: Herm Double
bigP = sym $ (3><3) [ 4.0e-1, 0.0,    0.0,
                      0.0,    4.0e-1, 0.0,
                      0.0,    0.0,    4.0e-4
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


lotkaVolterra'' :: [Double] -> SystemState Double -> OdeProblem
lotkaVolterra'' ts s = emptyOdeProblem
  { odeRhs = odeRhsPure $ \t x -> fromList (dzdt (meanRate {theta4 = coerce $ gamma s}) t (toList x))
  , odeJacobian = Nothing
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = vector [hares s, lynxes s]
  , odeSolTimes = vector ts
  , odeTolerances = defaultTolerances
  }

sol'' :: [Double] -> SystemState Double -> IO (Matrix Double)
sol'' ts s = do
  handleScribe <- mkHandleScribe ColorIfTerminal stderr (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts BOGACKI_SHAMPINE_4_2_3) (lotkaVolterra'' ts s)
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

stateUpdate :: Particles (SystemState Double) -> IO (Particles (SystemState Double))
stateUpdate ps = do
  qs <-  V.mapM (sol'' (take 21 us)) ps
  rr <- V.replicateM nParticles $
        R.sample $ R.rvar (RN.Normal 0.0 ((unSym bigP)!2!2))
  rrH <- V.replicateM nParticles $
         R.sample $ R.rvar (RN.Normal 0.0 ((unSym bigP)!0!0))
  rrL <- V.replicateM nParticles $
         R.sample $ R.rvar (RN.Normal 0.0 ((unSym bigP)!1!1))

  let newHLs = V.map (\m -> m!20) qs
      newGammas = V.zipWith (\p q -> gamma p + q) ps rr
      newHares  = V.zipWith (+) (V.map (!0) newHLs) rrH
      newLynxes = V.zipWith (+) (V.map (!1) newHLs) rrL
      newStates = V.zipWith3 (\h l m -> SystemState {hares = h, lynxes = l, gamma = m})
                             newHares newLynxes newGammas
  return newStates

measureOp :: Particles (SystemState Double) -> Particles (SystemObs Double)
measureOp = V.map (\s -> SystemObs { obsHares = hares s, obsLynxes = lynxes s})

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
