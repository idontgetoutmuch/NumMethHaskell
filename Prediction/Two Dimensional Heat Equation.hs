
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE OverloadedLists     #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE FlexibleInstances   #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE GADTs               #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeOperators       #-}

import Data.Maybe
import Data.Number.Symbolic
import qualified Data.Number.Symbolic as Sym
import Data.Proxy

import qualified Naperian as N
import qualified Data.Foldable as F
import           Control.Applicative ( liftA2 )
import qualified GHC.TypeLits as M
import           Data.Functor
import           Data.List.Split

import           Numeric.Sundials.ARKode.ODE
import           Numeric.LinearAlgebra

kx, ky :: Floating a => a
kx = 0.5
ky = 0.75

-- spatial mesh size
nx, ny :: Int
nx = 30
ny = 60

-- x mesh spacing
-- y mesh spacing
dx :: Floating a => a
dx = 1 / (fromIntegral nx - 1)

dy :: Floating a => a
dy = 1 / (fromIntegral ny - 1)

c1, c2 :: Floating a => a
c1 = kx/dx/dx
c2 = ky/dy/dy

cc4' :: forall b m n . (M.KnownNat m, M.KnownNat n, Num b) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] b
cc4' = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @m <&> (\(N.Fin x) ->
      N.viota @n <&> (\(N.Fin w) ->
      N.viota @m <&> (\(N.Fin v) ->
      N.viota @n <&> (\(N.Fin u) ->
      (f m n x w v u)))))
        where
          m = fromIntegral $ M.natVal (undefined :: Proxy m)
          n = fromIntegral $ M.natVal (undefined :: Proxy n)
          f m n i j k l | i == 0               = 0
                        | j == 0               = 0
                        | i == n - 1           = 0
                        | j == m - 1           = 0
                        | k == i - 1 && l == j = 1
                        | k == i     && l == j = -2
                        | k == i + 1 && l == j = 1
                        | otherwise            = 0

cc5' :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] a
cc5' = N.binary (*) (N.Scalar c2) cc4'

cc5Sym' :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a, Eq a) =>
          N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] (Sym a)
cc5Sym' = N.binary (*) (N.Scalar $ var "c2") cc4'

yy4' :: forall b m n . (M.KnownNat m, M.KnownNat n, Num b) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] b
yy4' = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @m <&> (\(N.Fin x) ->
      N.viota @n <&> (\(N.Fin w) ->
      N.viota @m <&> (\(N.Fin v) ->
      N.viota @n <&> (\(N.Fin u) ->
      (f m n x w v u)))))
        where
          m = fromIntegral $ M.natVal (undefined :: Proxy m)
          n = fromIntegral $ M.natVal (undefined :: Proxy n)
          f :: Int -> Int -> Int -> Int -> Int -> Int -> b
          f m n i j k l | i == 0                   = 0
                        | j == 0                   = 0
                        | i == n - 1               = 0
                        | j == m - 1               = 0
                        | k == i     && l == j - 1 = 1
                        | k == i     && l == j     = -2
                        | k == i     && l == j + 1 = 1
                        | otherwise                = 0

yy5' :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] a
yy5' = N.binary (*) (N.Scalar c1) yy4'

yy5Sym' :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a, Eq a) =>
           N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] (Sym a)
yy5Sym' = N.binary (*) (N.Scalar $ var "c1") yy4'

ccSym5 = cc5Sym' @Double @4 @5
yy5Sym = yy5Sym' @Double @4 @5

ccSym5

yy5Sym

fmap (N.elements . N.Prism . N.Prism . N.Scalar) $ N.elements $ N.crystal $ N.crystal $ N.binary (+) cc5Sym yy5Sym

{-# LANGUAGE DataKinds #-}
{-# LANGUAGE OverloadedLists       #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE GADTs               #-}

import Data.Maybe
import Data.Number.Symbolic
import Data.Proxy

import qualified Naperian as N
import qualified Data.Foldable as F
import           Control.Applicative ( liftA2 )
import qualified GHC.TypeLits as M

import           Numeric.Sundials.ARKode.ODE
import           Numeric.LinearAlgebra

x1, a, x2 :: Double
x1 = 0
a = 1.0
x2 = a

y1, y2 :: Double
y1 = 0.0
y2 = 1.0

bigT :: Double
bigT = 1000.0

n :: Int
n = 2

dx :: Double
dx = a / (fromIntegral n + 1)

dy :: Double
dy = a / (fromIntegral n + 1)

beta, s :: Double
beta = 1.0e-5
s = beta / (dx * dx)

kx, ky :: Double
kx = 0.5
ky = 0.75

c1, c2 :: Double
c1 = kx/dx/dx
c2 = ky/dy/dy

bigAA1 :: Matrix Double
bigAA1 = assoc (n * n, n * n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n * n - 1]
                                                      , j <- [i - n, i,  i + n]
                                                      , j `elem` [0 .. n * n -1]]
  where
    f (i, j) | i     == j = (-2.0) * c1
             | i - n == j = 1.0    * c1
             | i + n == j = 1.0    * c1
             | otherwise = error $ show (i, j)

bigAA2 :: Matrix Double
bigAA2 = diagBlock (replicate n bigA)
  where
    bigA :: Matrix Double
    bigA = assoc (n, n) 0.0 [((i, j), f (i, j)) | i <- [0 .. n - 1]
                                                , j <- [i-1..i+1]
                                                , j `elem` [0..n-1]]
      where
        f (i, j) | i     == j = (-2.0) * c2
                 | i - 1 == j = 1.0    * c2
                 | i + 1 == j = 1.0    * c2

bigAA :: Matrix Double
bigAA = bigAA1 + bigAA2

bigAA1

bigAA2

n

bigZZ1 :: Matrix Double
bigZZ1 = assoc (m * m, m * m) 0.0 [((i, j), f (i, j)) | i <- [0 .. m * m - 1]
                                                      , j <- [0 .. m * m - 1]]
  where
    m = n + 2
    f (i, j) | i     == 0     = 0.0
             | j     == 0     = 0.0
             | i     == j     = (-2.0) * c1
             | i - n == j     = 1.0    * c1
             | i + n == j     = 1.0    * c1
             | i     == n + 1 = 0.0
             | j     == n + 1 = 0.0
             | otherwise      = 0.0


bigZZ1

x :: forall m n . (M.KnownNat m, M.KnownNat n) => N.Vector n (N.Vector m (Sym Int))
x = (fromJust . N.fromList) $
    map (fromJust . N.fromList) ([[var $ (\(x,y) -> "A" ++ show x ++ "," ++ show y) (x,y) | y <- [1..m]] | x <- [1..n]] :: [[Sym Int]])
    where
      m = M.natVal (undefined :: Proxy m)
      n = M.natVal (undefined :: Proxy n)

u1 :: N.Hyper '[N.Vector 3, N.Vector 2] (Sym Int)
u1 = N.Prism $ N.Prism (N.Scalar x)

u1

y :: forall n . M.KnownNat n => N.Vector n (Sym Int)
y = (fromJust . N.fromList) $
    (map (var . ("v" ++) . show) [1..n ] :: [Sym Int])
    where
    n = M.natVal (undefined :: Proxy n)

u2 :: N.Hyper '[N.Vector 3] (Sym Int)
u2 = N.Prism (N.Scalar y)

N.innerH u1 u2
