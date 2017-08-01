{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE TypeOperators       #-}
{-# LANGUAGE ExplicitForAll      #-}

{-# LANGUAGE PolyKinds           #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE StandaloneDeriving  #-}

{-# LANGUAGE Rank2Types          #-}

module Main where

import Data.Random
import Data.Random.Distribution.Dirichlet
import Control.Monad
import Control.Monad.State
import Data.Random.Source.PureMT
import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.Static as S
import Numeric.LinearAlgebra.Devel
import Foreign.Storable
import OldFaithful

import GHC.TypeLits
import Data.Type.Equality
import Data.Proxy
import GHC.TypeLits.Compare

import Data.Maybe (fromJust)


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

-- data:
bigX :: Matrix Double
bigX = tr (fromLists [map eruptions oldFaithful, map waiting oldFaithful])

-- params:
bigK, bigN, bigD :: Int
bigK  = 6
bigN = fst $ size bigX
bigD = snd $ size bigX

ones :: Num a => [a]
ones = 1 : ones

preZ :: RVar [[Double]]
preZ = replicateM bigN $ dirichlet (take bigK ones)

bigZ :: RVar (Matrix Double)
bigZ = liftM (bigN >< bigK) (liftM concat preZ)

bigR :: Matrix Double
bigR = evalState (sample bigZ) (pureMT 42)

xd :: Matrix Double
xd = konst 0.0 (bigK, bigD)

foo0, foo1 :: Matrix Double
foo0 = cmap (bigR `atIndex` (0, 0) *) (bigX ? [0])
foo1 = cmap (bigR `atIndex` (1, 0) *) (bigX ? [1])

foos, bars :: Int -> Double
foos m = sum $ map (`atIndex` (0,0)) $ map (\n -> cmap (bigR `atIndex` (n, m) *) (bigX ? [n])) [0..bigN-1]
bars m = sum $ map (`atIndex` (0,1)) $ map (\n -> cmap (bigR `atIndex` (n, m) *) (bigX ? [n])) [0..bigN-1]

test2 = S.withMatrix bigR f
  where
    f :: forall n' k . (KnownNat n', KnownNat k) => S.L n' k -> ((Matrix Double, Matrix Double, Vector Double, [Vector Double]), Matrix Double)
    f r = ( S.withMatrix (fromLists [map eruptions oldFaithful, map waiting oldFaithful]) g
          , S.extract nK
          )
      where
        g :: forall d n . (KnownNat d, KnownNat n) => S.L d n -> (Matrix Double, Matrix Double, Vector Double, [Vector Double])
        g x = case sameNat (Proxy :: Proxy n') (Proxy :: Proxy n) of
                Nothing -> error "Incompatible matrices"
                Just Refl -> (S.extract xbar, S.extract ((S.tr x) - xbar1), S.extract xxbar1, map S.extract xxbars)
                               where
                                 r1 :: S.L d n
                                 r1 = repmat' (Proxy :: Proxy d) $
                                      S.row $ getFirst $ S.tr r
                                 rs :: [S.L d n]
                                 rs = map (repmat' (Proxy :: Proxy d)) $
                                       map (S.row . fromJust . S.create) $
                                       toRows $ S.extract $ S.tr r
                                 xxbar1 :: S.R d
                                 xxbar1 = fromJust $ S.create $ fromList $
                                          map (foldVector (+) 0) $
                                          toRows $ S.extract r1 .* S.extract x
                                 xxbar :: S.L d n -> S.L d n -> S.R d
                                 xxbar rr1 xx = fromJust $ S.create $ fromList $
                                               map (foldVector (+) 0) $
                                               toRows $ S.extract rr1 .* S.extract xx
                                 xxbars :: [S.R d]
                                 xxbars = zipWith xxbar rs (repeat x)
                                 xbar :: S.L k d
                                 xbar = (S.tr r) S.<> (S.tr x)
                                 xbar1 :: S.L n d
                                 xbar1 = fromJust $
                                         S.create $
                                         fromRows $
                                         replicate l (S.extract $ getFirst xbar)

        l = fst $ S.size r
        nK :: S.L 1 k
        nK = S.row (S.vector $ take l ones) S.<> r

infixl 7 .*

repmat' :: forall m n . (KnownNat m, KnownNat n) => Proxy n -> S.L 1 m -> S.L n m
repmat' p x = fromJust $ S.create z
  where
    y = S.extract x
    z :: Matrix Double
    z = foldr (===) y (replicate (n - 1) y)
    n = fromIntegral $ natVal p

class PointWise a where
  (.*) :: a -> a -> a

instance (Num a, Storable a) => PointWise (Vector a) where
  (.*) = zipVectorWith (*)

instance (Num a, Storable a, Element a) => PointWise (Matrix a) where
  (.*) = liftMatrix2 (.*)

getFirst :: forall m n . (KnownNat m, KnownNat n) => S.L m n -> S.R n
getFirst x = case (Proxy :: Proxy 1) %<=? (Proxy :: Proxy m) of
          LE  Refl -> S.unrow $ fst $ S.splitRows x
          NLE _    -> 0

main :: IO ()
main = do
  let ((xbar, _ms1, xx1, xxs), nK) = test2
      vk = cmap (v0 +) nK
      _betak = cmap (beta0 +) nK
  putStrLn $ show $ size nK
  putStrLn $ show (xbar `atIndex` (0,0))
  putStrLn $ show xx1
  putStrLn $ show xxs
  putStrLn $ show nK
  putStrLn $ show vk
  putStrLn $ show xbar
  putStrLn "Hello, Haskell!"
