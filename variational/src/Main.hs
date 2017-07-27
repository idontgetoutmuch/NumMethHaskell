{-# OPTIONS_GHC -Wall #-}

module Main where

import Data.Random
import Data.Random.Distribution.Dirichlet
import Control.Monad
import Numeric.LinearAlgebra
import OldFaithful

-- hyperparams:
alpha0, beta0, v0 :: Double
alpha0 = 0.1  -- prior coefficient count (for Dirichlet)
beta0 = 1e-20 -- variance of mean (smaller: broader the means)
v0 = 20       -- degrees of freedom in inverse wishart FIXME: this
              -- needs some explanation as the Python seems to set
              -- this to 3!

m0 :: Vector Double
m0 = vector [0.0, 0.0]               -- prior mean
w0 :: Matrix Double
w0 = (2><2) [1.0e0, 0.0, 0.0, 1.0e0] -- prior cov (bigger: smaller covariance)

-- params:
bigK, bigN :: Int
bigK  = 6
bigN = length(oldFaithful)

ones :: Num a => [a]
ones = 1 : ones

preZ :: RVar [[Double]]
preZ = replicateM bigN $ dirichlet (take bigK ones)

bigZ :: RVar (Matrix Double)
bigZ = liftM (bigN >< bigK) (liftM concat preZ)

main :: IO ()
main = do
  bigR <- sample bigZ
  let nK = asRow ((rows bigR) |> ones) <> bigR
      vk = cmap (v0 +) nK
      betak = cmap (beta0 +) nK
      xbar = tr bigR <> tr (fromLists [map eruptions oldFaithful, map waiting oldFaithful])
  putStrLn $ show nK
  putStrLn $ show vk
  putStrLn $ show xbar
  putStrLn "Hello, Haskell!"
