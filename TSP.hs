{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module TSP ( test ) where

import Prelude hiding ( zipWith, sum, foldl, foldr, scanl, length, null, tail, zip, map, reverse, concat, take )
import qualified Prelude as P
import Data.Vector.Unboxed hiding ( replicateM )
import Data.Random.Source.PureMT
import Data.Random
import Control.Monad.State hiding ( modify )

import Debug.Trace

nNodes :: Int
nNodes = 6

xs, ys :: Vector Float
xs = generate nNodes ((*1.0) . fromIntegral)
ys = generate nNodes ((*0.0) . fromIntegral)

distance :: Int -> Int -> Float
distance i j = sqrt $ (xs!i - xs!j)^2 + (ys!i - ys!j)^2

totalDistance :: Vector Int -> Float
totalDistance v = sum $ map (uncurry distance) ms
  where
    ms = zip v (tail v)

reverseBetweenPair :: Int -> Int -> Vector Int -> Vector Int
reverseBetweenPair i j v = concat [beginning, reverse middle, end]
  where
    k = length v

    beginning = slice 0       i           v
    middle    = slice i       (j - i + 1) v
    end       = slice (j + 1) (k - j - 1) v

incDistance :: Vector Int -> Int -> Int -> Float
incDistance v i j = d + c - b - a
  where
    a = distance (v!(i - 1)) (v!i)
    b = distance (v!j)       (v!((j + 1) `mod` nNodes))
    c = distance (v!(i - 1)) (v!j)
    d = distance (v!i)       (v!((j + 1) `mod` nNodes))

expDv :: Int -> Int -> Vector Int -> Float -> Float -> Float -> Float
expDv i1 i2 v kB j t = exp(-j * (incDistance v i1 i2) / (kB * t))

randomUpdates :: Int -> Vector (Int, Int, Float)
randomUpdates m =
  fromList $
  evalState (replicateM m x)
  (pureMT 1)
  where
    x = do r <- sample (uniform (1 :: Int)    (nNodes - 2))
           c <- sample (uniform (r + 1)       (nNodes - 1))
           v <- sample (uniform (0 :: Float)           1.0)
           return (r, c, v)

kB, couplingConstant :: Float
kB               = 1.0
couplingConstant = 1.0

data McState = McState {
  mcRoute :: !(Vector Int)
  }
  deriving Show

initMcState :: McState
initMcState = McState {
  mcRoute = fromList [0,4,2,3,1,5,0]
  }

singleUpdate :: Float -> McState -> (Int, Int, Float) -> McState
singleUpdate t u (i, j, r) =
  trace (show i P.++ show j P.++ ": " P.++
         show (totalDistance v) P.++ ", " P.++
         show (totalDistance (reverseBetweenPair i j v)) P.++ ", " P.++
         show (incDistance v i j) P.++ ", " P.++
         show p P.++ ", " P.++
         show r) $ if incDistance v i j <= 0 || p > r
  then
    McState { mcRoute = reverseBetweenPair i j v }
  else
    u
  where
    v = mcRoute u
    p = expDv i j v kB couplingConstant t

test :: McState
test = foldl (singleUpdate 1.0) initMcState (randomUpdates 100)