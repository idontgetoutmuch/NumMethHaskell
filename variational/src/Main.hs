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

import Data.Maybe (fromJust)

import qualified Naperian as N

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
    f :: forall n' k . (KnownNat n', KnownNat k) => S.L n' k -> ((Matrix Double, Matrix Double, [Vector Double]), Matrix Double)
    f r = ( S.withMatrix (fromLists [map eruptions oldFaithful, map waiting oldFaithful]) g
          , S.extract nK
          )
      where
        g :: forall d n . (KnownNat d, KnownNat n) => S.L d n -> (Matrix Double, Matrix Double, [Vector Double])
        g x = case sameNat (Proxy :: Proxy n') (Proxy :: Proxy n) of
                Nothing -> error "Incompatible matrices"
                Just Refl -> (S.extract xbars, S.extract ((S.tr x) - (mkXbarss (xxbars!!0))), map S.extract xxbars)
                               where
                                 eek :: S.L d n
                                 eek = x
                                 rs :: [S.L d n] -- k
                                 rs = map (repmat' (Proxy :: Proxy d)) $
                                       map (S.row . fromJust . S.create) $
                                       toRows $ S.extract $ S.tr r
                                 xxbars :: [S.R d] -- k
                                 xxbars = map sumRows $ rs .* (repeat x)
                                 xxbars' :: [S.R d] -- k
                                 xxbars' = zipWith (\n u -> S.dvmap (/n) u)
                                                   (toList (S.extract $ S.unrow nK))
                                                   (map sumRows $ rs .* (repeat x))
                                 xbars :: S.L k d
                                 xbars = (S.tr r) S.<> (S.tr x)
                                 mkXbarss :: S.R d -> S.L n d
                                 mkXbarss u = fromJust $
                                              S.create $
                                              fromRows $
                                              replicate l (S.extract u)
                                 foo :: [S.L n d] -- k
                                 foo = map mkXbarss xxbars
                                 bar :: [S.L n d] -- k
                                 bar = map (S.tr x -) foo
                                 baz :: [[S.R d]] -- k x n?
                                 baz = map (map (fromJust . S.create) . toRows . S.extract) bar
                                 urk :: [S.L d d] -- k?
                                 urk = map sum $ map (map (\x -> (S.col x) S.<> (S.row x))) baz
        l = fst $ S.size r
        nK :: S.L 1 k
        nK = S.row (S.vector $ take l ones) S.<> r

test3 :: ([Vector Double], Matrix Double)
test3 = S.withMatrix bigR f
  where
    f :: forall n' k . (KnownNat n', KnownNat k) =>
                       S.L n' k -> ([Vector Double], Matrix Double)
    f r = ( S.withMatrix (fromLists [map eruptions oldFaithful, map waiting oldFaithful]) g
          , S.extract nK
          )
      where
        g :: forall d n . (KnownNat d, KnownNat n) =>
                          S.L d n -> [Vector Double]
        g x = case sameNat (Proxy :: Proxy n') (Proxy :: Proxy n) of
                Nothing -> error "Incompatible matrices"
                Just Refl -> map S.extract xbars
                               where
                                 rs :: [S.L d n] -- k
                                 rs = map (repmat' (Proxy :: Proxy d)) $
                                       map (S.row . fromJust . S.create) $
                                       toRows $ S.extract $ S.tr r
                                 xbars :: [S.R d] -- k
                                 xbars = zipWith (\n u -> S.dvmap (/n) u)
                                                 (toList (S.extract $ S.unrow nK))
                                                 (map sumRows $ rs .* (repeat x))
        l = fst $ S.size r
        nK :: S.L 1 k
        nK = S.row (S.vector $ take l ones) S.<> r

infixl 7 .*

class Pointwise a where
  (.*) :: a -> a -> a

instance (Num a, Storable a) => Pointwise (Vector a) where
  (.*) = zipVectorWith (*)

instance (Num a, Storable a, Element a) => Pointwise (Matrix a) where
  (.*) = liftMatrix2 (.*)

instance (KnownNat n, KnownNat m) => Pointwise (S.L m n) where
  x .* y = fromJust $ S.create $ (S.extract x) .* (S.extract y)

instance Pointwise a => Pointwise [a] where
  (.*) = zipWith (.*)

sumRows ::  forall d n . (KnownNat d, KnownNat n) => S.L d n -> S.R d
sumRows = fromJust . S.create . fromList .
          map (foldVector (+) 0) .
          toRows . S.extract

repmat' :: forall m n . (KnownNat m, KnownNat n) => Proxy n -> S.L 1 m -> S.L n m
repmat' p x = fromJust $ S.create z
  where
    y = S.extract x
    z :: Matrix Double
    z = foldr (===) y (replicate (n - 1) y)
    n = fromIntegral $ natVal p

bigRH :: N.Matrix 272 6 Double
bigRH = (fromJust . N.fromList) $
        map (fromJust . N.fromList) $
        map toList $
        toRows bigR

bigRH' :: N.Hyper '[N.Vector 6, N.Vector 272] Double
bigRH' = N.Prism (N.Prism (N.Scalar bigRH))

bigXH :: N.Matrix 272 2 Double
bigXH = (fromJust . N.fromList) $
        map (fromJust . N.fromList) $
        map toList $
        toRows bigX

nKK :: N.Hyper '[N.Vector 6] Double
nKK = N.reduceBy ((+),0) $ N.transposeH (N.Prism (N.Prism (N.Scalar bigRH)))

nK :: N.Hyper '[N.Vector 6] Double
nK = N.foldrH (+) 0 $
     N.transposeH (N.Prism (N.Prism (N.Scalar bigRH)))

xbar0 :: N.Hyper '[N.Vector 2] Double
xbar0 = N.foldrH (+) 0 $
        N.binary (*) (N.Prism $ N.Scalar $ N.lookup (N.transpose bigRH) (N.Fin 0))
                     (N.Prism (N.Prism (N.Scalar (N.transpose bigXH))))

bigXHs :: N.Vector 2 (N.Vector 6 (N.Vector 272 Double))
bigXHs = N.transpose $ N.replicate (N.transpose bigXH)

altBigRH :: N.Vector 6 (N.Vector 272 Double)
altBigRH = N.transpose bigRH

xbars :: N.Hyper '[N.Vector 6, N.Vector 2] Double
xbars = N.foldrH (+) 0 $
        N.binary (*) (N.Prism (N.Prism (N.Scalar altBigRH)))
                     (N.Prism $ N.Prism $ N.Prism $ N.Scalar bigXHs)

xbars' :: N.Hyper '[N.Vector 6, N.Vector 2] Double
xbars' = N.binary (/) xbars nK

diff1 :: N.Hyper '[N.Vector 272, N.Vector 6, N.Vector 2] Double
diff1 = case xbars' of
          N.Prism (N.Prism (N.Scalar xbarU)) ->
            N.binary (-) b a
            where
              a :: N.Hyper '[N.Vector 272, N.Vector 6, N.Vector 2] Double
              a = N.Prism (N.Prism (N.Prism (N.Scalar bigXHs)))
              b :: N.Hyper '[N.Vector 272, N.Vector 6, N.Vector 2] Double
              b = N.transposeH $ N.transposeH' $
                  N.Prism (N.Prism (N.Prism (N.Scalar (N.replicate xbarU))))

diff2 :: N.Hyper '[N.Vector 272, N.Vector 6, N.Vector 2] Double
diff2 = N.binary (*) a diff1
  where
    a :: N.Hyper '[N.Vector 272, N.Vector 6, N.Vector 2] Double
    a = N.transposeH $ N.Prism (N.Prism (N.Prism (N.Scalar (N.replicate bigRH))))

bigS :: N.Hyper '[N.Vector 6] (N.Vector 2 (N.Vector 2 Double))
bigS = N.binary N.matrix a b
  where
    a :: N.Hyper '[N.Vector 6] (N.Vector 2 (N.Vector 272 Double))
    a =  N.crystal $ N.crystal $ N.transposeH' diff1
    b :: N.Hyper '[N.Vector 6] (N.Vector 272 (N.Vector 2 Double))
    b = N.unary N.transpose $ N.crystal $ N.crystal $ N.transposeH' diff2

main :: IO ()
main = do
  let ((xbar, ms1, xxs), nK) = test2
      vk = cmap (v0 +) nK
      _betak = cmap (beta0 +) nK
  putStrLn $ show $ size nK
  putStrLn $ show ms1
  putStrLn $ show nK
  putStrLn $ show vk
  putStrLn $ show xbar
  putStrLn $ show xxs
  putStrLn "Hello, Haskell!"
