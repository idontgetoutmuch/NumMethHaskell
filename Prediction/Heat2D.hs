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
import Data.Proxy

import qualified Naperian as N
import qualified Data.Foldable as F
import           Control.Applicative ( liftA2 )
import qualified GHC.TypeLits as M
import           Data.Functor

import           Numeric.Sundials.ARKode.ODE
import           Numeric.LinearAlgebra

import Debug.Trace

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

x :: forall m n . (M.KnownNat m, M.KnownNat n) => N.Vector n (N.Vector m (Sym Int))
x = (fromJust . N.fromList) $
    map (fromJust . N.fromList) $
    ([[var $ (\(x,y) -> "B_{" ++ show x ++ "," ++ show y ++ "}") (x,y) | y <- [1..m]] | x <- [1..n]] :: [[Sym Int]])
    where
      m = M.natVal (undefined :: Proxy m)
      n = M.natVal (undefined :: Proxy n)

u1 :: N.Hyper '[N.Vector 3, N.Vector 2] (Sym Int)
u1 = N.Prism $ N.Prism (N.Scalar x)

y :: forall n . M.KnownNat n => N.Vector n (Sym Int)
y = (fromJust . N.fromList) $
    (map (var . ("v" ++) . show) [1..n ] :: [Sym Int])
    where
    n = M.natVal (undefined :: Proxy n)

u2 :: N.Hyper '[N.Vector 3] (Sym Int)
u2 = N.Prism (N.Scalar y)

test = N.innerH u1 u2

-- f :: Int -> Int -> Int -> Int -> Int -> Int -> Sym Int
-- f m n i j k l | i == 1               = trace ("i == 1: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "0"
--               | j == 1               = trace ("j == 1: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "0"
--               | i == n               = trace ("i == n: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "0"
--               | j == m               = trace ("j == m: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "0"
--               | k == i - 1 && l == j = trace ("diag -1: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "1"
--               | k == i     && l == j = trace ("diag  0: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "2"
--               | k == i + 1 && l == j = trace ("diag +1: " ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "1"
--               | otherwise            = trace ("other: "   ++ show m ++ "," ++ show n ++ "," ++ show i ++ "," ++ show j ++ "," ++ show k ++ "," ++ show l) $ var "0"

-- f :: Int -> Int -> Int -> Int -> Int -> Int -> Sym Int
-- f m n i j k l | i == 1               = var "0"
--               | j == 1               = var "0"
--               | i == n               = var "0"
--               | j == m               = var "0"
--               | k == i - 1 && l == j = var "1"
--               | k == i     && l == j = var "2"
--               | k == i + 1 && l == j = var "1"
--               | otherwise            = var "0"

f :: Int -> Int -> Int -> Int -> Int -> Int -> Sym Int
f m n i j k l | i == 1               = 0
              | j == 1               = 0
              | i == n               = 0
              | j == m               = 0
              | k == i - 1 && l == j = 1
              | k == i     && l == j = 2
              | k == i + 1 && l == j = 1
              | otherwise            = 0

x4 :: forall m n p q . (M.KnownNat m, M.KnownNat n, M.KnownNat p, M.KnownNat q) =>
                       N.Vector q (N.Vector p (N.Vector n (N.Vector m (Sym Int))))
x4 =
  (fromJust . N.fromList) $
  map (fromJust . N.fromList) $
  map (map (fromJust . N.fromList)) $
  map (map (map (fromJust . N.fromList))) $
  map (map (map (map (var . (\(u,v,w,x) -> "A_{" ++ show x ++ "," ++ show w ++ "," ++ show v ++ "," ++ show u ++ "}"))))) $
  -- ([[[[ (x,w,v,u) | u <- [1..m]] | v <- [1..n]] | w <- [1..p]] | x <- [1..q]] :: [[[[(Int, Int, Int, Int)]]]])
  ([[[[ (u,v,w,x) | u <- [1..m]] | v <- [1..n]] | w <- [1..p]] | x <- [1..q]] :: [[[[(Int, Int, Int, Int)]]]])
  where
    m, n, p, q :: Int
    m = fromIntegral $ M.natVal (undefined :: Proxy m)
    n = fromIntegral $ M.natVal (undefined :: Proxy n)
    p = fromIntegral $ M.natVal (undefined :: Proxy p)
    q = fromIntegral $ M.natVal (undefined :: Proxy q)

x4a :: forall m n p q . (M.KnownNat m, M.KnownNat n, M.KnownNat p, M.KnownNat q) =>
                        N.Vector q (N.Vector p (N.Vector n (N.Vector m (Sym Int))))
x4a =
  (fromJust . N.fromList) $
  map (fromJust . N.fromList) $
  map (map (fromJust . N.fromList)) $
  map (map (map (fromJust . N.fromList))) $
  map (map (map (map ((\(u,v,w,x) -> f m n x w v u))))) $
  -- map (map (map (map ((\(u,v,w,x) -> f m n u v w x))))) $
  -- ([[[[ (x,w,v,u) | u <- [1..m]] | v <- [1..n]] | w <- [1..p]] | x <- [1..q]] :: [[[[(Int, Int, Int, Int)]]]])
  ([[[[ (u,v,w,x) | u <- [1..m]] | v <- [1..n]] | w <- [1..p]] | x <- [1..q]] :: [[[[(Int, Int, Int, Int)]]]])
  where
    m, n, p, q :: Int
    m = fromIntegral $ M.natVal (undefined :: Proxy m)
    n = fromIntegral $ M.natVal (undefined :: Proxy n)
    p = fromIntegral $ M.natVal (undefined :: Proxy p)
    q = fromIntegral $ M.natVal (undefined :: Proxy q)

v1 :: N.Hyper '[N.Vector 3, N.Vector 4, N.Vector 3, N.Vector 4] (Sym Int)
v1 = N.Prism $ N.Prism $ N.Prism $ N.Prism (N.Scalar x4)

v1a :: N.Hyper '[N.Vector 3, N.Vector 4, N.Vector 3, N.Vector 4] (Sym Int)
v1a = N.Prism $ N.Prism $ N.Prism $ N.Prism (N.Scalar x4a)

xa :: forall m n . (M.KnownNat m, M.KnownNat n) => N.Vector n (N.Vector m (Sym Int))
xa = (fromJust . N.fromList) $
    map (fromJust . N.fromList) $
    ([[var $ (\(x,y) -> "B_{" ++ show x ++ "," ++ show y ++ "}") (x,y) | y <- [1..m]] | x <- [1..n]] :: [[Sym Int]])
    where
      m = M.natVal (undefined :: Proxy m)
      n = M.natVal (undefined :: Proxy n)

v2 :: N.Hyper '[N.Vector 3, N.Vector 4] (Sym Int)
v2 = N.Prism $ N.Prism (N.Scalar xa)

-- x4b :: forall m n p q . (M.KnownNat m, M.KnownNat n, M.KnownNat p, M.KnownNat q) =>
--                         N.Vector q (N.Vector p (N.Vector n (N.Vector m (Sym Int))))
-- x4b = N.viota @m <$> (\(N.Fin x) ->
--       N.viota @n <$> (\(N.Fin w) ->
--       N.viota @p <$> (\(N.Fin v) ->
--       N.viota @q <$> (\(N.Fin u) ->
--       "A_{" ++ show x ++ "," ++ show w ++ "," ++ show v ++ "," ++ show u ++ "}"))))

x4b = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @2 <&> (\(N.Fin x) ->
      N.viota @3 <&> (\(N.Fin w) ->
      N.viota @2 <&> (\(N.Fin v) ->
      N.viota @3 <&> (\(N.Fin u) ->
      (var ("A_{" ++ show x ++ "," ++ show w ++ "," ++ show v ++ "," ++ show u ++ "}"))))))

y4b = N.Prism $ N.Prism $ N.Scalar $
      N.viota @2 <&> (\(N.Fin x) ->
      N.viota @3 <&> (\(N.Fin w) ->
      (var ("B_{" ++ show x ++ "," ++ show w ++ "}"))))

correct = N.foldrH (+) 0 $ N.foldrH (+) 0 $ N.binary (*) v1a v2

correctA = N.foldrH (+) 0 $ N.foldrH (+) 0 $ N.binary (*) x4b y4b

class HyperLift f fs where
  hyper :: (N.Shapely fs, N.Dimension f) => f a -> N.Hyper (f : fs) a

instance HyperLift f '[] where
  hyper = N.Prism . N.Scalar

instance (N.Shapely fs, HyperLift f fs) => HyperLift f (f : fs) where
  hyper = N.Prism . (\x -> (hyper $ pure x))

aa4 :: forall m n a . (M.KnownNat m, M.KnownNat n) =>
                      N.Hyper '[N.Vector m, N.Vector n, N.Vector m, N.Vector n] (Sym a)
aa4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @n <&> (\(N.Fin x) ->
      N.viota @m <&> (\(N.Fin w) ->
      N.viota @n <&> (\(N.Fin v) ->
      N.viota @m <&> (\(N.Fin u) ->
      (var ("A_{" ++ show x ++ "," ++ show w ++ "," ++ show v ++ "," ++ show u ++ "}"))))))

bb4 :: N.Hyper '[N.Vector 3, N.Vector 4, N.Vector 3, N.Vector 4] (Sym Int)
bb4 = aa4

aa2 :: forall m n a . (M.KnownNat m, M.KnownNat n) =>
                      N.Hyper '[N.Vector m, N.Vector n] (Sym a)
aa2 = N.Prism $ N.Prism $ N.Scalar $
      N.viota @n <&> (\(N.Fin x) ->
      N.viota @m <&> (\(N.Fin w) ->
      (var ("B_{" ++ show x ++ "," ++ show w ++ "}"))))

bb2 :: N.Hyper '[N.Vector 3, N.Vector 4] (Sym Int)
bb2 = aa2

cc4 :: forall m n a . (M.KnownNat m, M.KnownNat n) =>
                      N.Hyper '[N.Vector m, N.Vector n, N.Vector m, N.Vector n] (Sym Int)
cc4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @n <&> (\(N.Fin x) ->
      N.viota @m <&> (\(N.Fin w) ->
      N.viota @n <&> (\(N.Fin v) ->
      N.viota @m <&> (\(N.Fin u) ->
      (f m n x w v u)))))
        where
          m = fromIntegral $ M.natVal (undefined :: Proxy m)
          n = fromIntegral $ M.natVal (undefined :: Proxy n)
          f :: Int -> Int -> Int -> Int -> Int -> Int -> Sym Int
          f m n i j k l | i == 0               = 0
                        | j == 0               = 0
                        | i == n - 1           = 0
                        | j == m - 1           = 0
                        | k == i - 1 && l == j = 1
                        | k == i     && l == j = -2
                        | k == i + 1 && l == j = 1
                        | otherwise            = 0

dd4 :: N.Hyper '[N.Vector (2 M.+ 1), N.Vector (3 M.+ 1), N.Vector (2 M.+ 1), N.Vector (3 M.+ 1)] (Sym Int)
dd4 = cc4

dd4' :: N.Hyper '[N.Vector (2 M.+ 1), N.Vector (3 M.+ 1), N.Vector (2 M.+ 1), N.Vector (3 M.+ 1)] Double
dd4' = N.binary (*) (N.Scalar c1) (fmap fromIntegral dd4)


almostTheFunction = N.foldrH (+) 0 $ N.foldrH (+) 0 $ N.binary (*) dd4 bb2

yy4 :: forall m n a . (M.KnownNat m, M.KnownNat n) =>
                      N.Hyper '[N.Vector m, N.Vector n, N.Vector m, N.Vector n] (Sym Int)
yy4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @n <&> (\(N.Fin x) ->
      N.viota @m <&> (\(N.Fin w) ->
      N.viota @n <&> (\(N.Fin v) ->
      N.viota @m <&> (\(N.Fin u) ->
      (f m n x w v u)))))
        where
          m = fromIntegral $ M.natVal (undefined :: Proxy m)
          n = fromIntegral $ M.natVal (undefined :: Proxy n)
          f :: Int -> Int -> Int -> Int -> Int -> Int -> Sym Int
          f m n i j k l | i == 0                   = 0
                        | j == 0                   = 0
                        | i == n - 1               = 0
                        | j == m - 1               = 0
                        | k == i     && l == j - 1 = 1
                        | k == i     && l == j     = -2
                        | k == i     && l == j + 1 = 1
                        | otherwise                = 0

zz4 :: N.Hyper '[N.Vector (2 M.+ 1), N.Vector (3 M.+ 1), N.Vector (2 M.+ 1), N.Vector (3 M.+ 1)] (Sym Int)
zz4 = yy4

zz4' :: N.Hyper '[N.Vector (2 M.+ 1), N.Vector (3 M.+ 1), N.Vector (2 M.+ 1), N.Vector (3 M.+ 1)] Double
zz4' = N.binary (*) (N.Scalar c2) (fmap fromIntegral zz4)

almostTheFunction' = N.foldrH (+) 0 $ N.foldrH (+) 0 $ N.binary (*) zz4 bb2


ff :: N.Hyper '[N.Vector 3, N.Vector 4] Double -> N.Hyper '[N.Vector 3, N.Vector 4] Double
ff x = N.foldrH (+) 0 $ N.foldrH (+) 0 $ N.binary (*) (N.binary (+) dd4' zz4') x
