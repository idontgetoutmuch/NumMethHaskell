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

import Debug.Trace


kx, ky :: Floating a => a
kx = 0.5
ky = 0.75

-- x mesh spacing
-- y mesh spacing
dx :: Floating a => a
dx = 1 / (fromIntegral nx - 1)

dy :: Floating a => a
dy = 1 / (fromIntegral ny - 1)

c1, c2 :: Floating a => a
c1 = kx/dx/dx
c2 = ky/dy/dy

preA2 :: forall b m n . (M.KnownNat m, M.KnownNat n, Num b) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] b
preA2 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
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

a2Num :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a) =>
        N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] a
a2Num = N.binary (*) (N.Scalar c2) preA2

a2 :: forall a m n . (M.KnownNat m, M.KnownNat n, Floating a, Eq a) =>
          N.Hyper '[N.Vector n, N.Vector m, N.Vector n, N.Vector m] (Sym a)
a2 = N.binary (*) (N.Scalar $ var "c2") preA2

yy5 :: Floating a => N.Hyper '[N.Vector 5, N.Vector 4, N.Vector 5, N.Vector 4] a
yy5 = N.binary (*) (N.Scalar c1) yy4

yy5Sym :: forall a . (Floating a, Eq a) =>
          N.Hyper '[N.Vector 5, N.Vector 4, N.Vector 5, N.Vector 4] (Sym a)
yy5Sym = N.binary (*) (N.Scalar $ var "c1") yy4

yy4 :: forall b . Num b => N.Hyper '[N.Vector 5, N.Vector 4, N.Vector 5, N.Vector 4] b
yy4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
      N.viota @4 <&> (\(N.Fin x) ->
      N.viota @5 <&> (\(N.Fin w) ->
      N.viota @4 <&> (\(N.Fin v) ->
      N.viota @5 <&> (\(N.Fin u) ->
      (f m n x w v u)))))
        where
          m = 5
          n = 4
          f :: Int -> Int -> Int -> Int -> Int -> Int -> b
          f m n i j k l | i == 0                   = 0
                        | j == 0                   = 0
                        | i == n - 1               = 0
                        | j == m - 1               = 0
                        | k == i     && l == j - 1 = 1
                        | k == i     && l == j     = -2
                        | k == i     && l == j + 1 = 1
                        | otherwise                = 0

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

-- spatial mesh size
nx, ny :: Int
nx = 3
ny = 6

bigA :: Matrix Double
bigA = fromLists $
       fmap (N.elements . N.Prism . N.Prism . N.Scalar) $
       N.elements $ N.crystal $ N.crystal $ N.binary (+) cc5 yy5
  where
    cc5, yy5 :: N.Hyper '[N.Vector 3, N.Vector 6, N.Vector 3, N.Vector 6] Double
    cc5 = a2Num
    yy5 = yy5'

b :: Vector Double
b = fromList $
    N.elements (h' :: N.Hyper '[N.Vector 3, N.Vector 6] Double)

safeUpdate :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a, Show a) =>
                             N.Hyper '[N.Vector m, N.Vector n] a ->
                             N.Hyper '[N.Vector m, N.Vector n] a
safeUpdate x = -- trace (show cc5 ++ "\n" ++ show yy5 ++ "\n" ++ show x) $
               N.foldrH (+) 0 $
               N.foldrH (+) 0 $
               N.binary (*) (N.binary (+) cc5 yy5) x
  where
    cc5 = N.binary (*) (N.Scalar c2) (fmap fromIntegral cc4)
    cc4 :: Num b => N.Hyper '[N.Vector m, N.Vector n, N.Vector m, N.Vector n] b
    cc4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
          N.viota @n <&> (\(N.Fin x) ->
          N.viota @m <&> (\(N.Fin w) ->
          N.viota @n <&> (\(N.Fin v) ->
          N.viota @m <&> (\(N.Fin u) ->
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

    yy5 = N.binary (*) (N.Scalar c1) (fmap fromIntegral yy4)
    yy4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
          N.viota @n <&> (\(N.Fin x) ->
          N.viota @m <&> (\(N.Fin w) ->
          N.viota @n <&> (\(N.Fin v) ->
          N.viota @m <&> (\(N.Fin u) ->
          (f m n x w v u)))))
            where
              m = fromIntegral $ M.natVal (undefined :: Proxy m)
              n = fromIntegral $ M.natVal (undefined :: Proxy n)
              f m n i j k l | i == 0                   = 0
                            | j == 0                   = 0
                            | i == n - 1               = 0
                            | j == m - 1               = 0
                            | k == i     && l == j - 1 = 1
                            | k == i     && l == j     = -2
                            | k == i     && l == j + 1 = 1
                            | otherwise                = 0

hyperfy :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a) =>
                          [a] -> N.Hyper '[N.Vector m, N.Vector n] a
hyperfy = N.Prism . N.Prism . N.Scalar .
          (fromJust . N.fromList) .
          map (fromJust . N.fromList) .
          chunksOf m
            where
              m = fromIntegral $ M.natVal (undefined :: Proxy m)

update :: Int -> Int -> [Double] -> [Double]
update k n xs =
  case M.someNatVal $ fromIntegral n of
    Nothing -> error "static/dynamic mismatch"
    Just (M.SomeNat (_ :: Proxy n)) ->
      case M.someNatVal $ fromIntegral k of
        Nothing -> error "static/dynamic mismatch"
        Just (M.SomeNat (_ :: Proxy k)) ->
          N.elements $ safeUpdate b
          where
            b :: N.Hyper '[N.Vector k, N.Vector n] Double
            b = hyperfy xs

h :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a) =>
                    N.Hyper '[N.Vector m, N.Vector n] (Sym a)
h = N.Prism $ N.Prism $ N.Scalar $
    N.viota @n <&> (\(N.Fin x) ->
    N.viota @m <&> (\(N.Fin w) ->
    (var ("u_{" ++ show x ++ "," ++ show w ++ "}"))))

h34 :: Floating a => N.Hyper '[N.Vector 3, N.Vector 4] (Sym a)
h34 = h

bigH :: Vector Double
bigH = assoc (nx * ny) 0.0 [(i + j * nx, f (i, j)) | i <- [0 .. nx - 1]
                                                   , j <- [0 .. ny - 1]]
 where
   f (i, j) = sin (pi * (fromIntegral i) * dx)
            * sin (2 * pi * (fromIntegral j) * dy)

h' :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a) =>
                     N.Hyper '[N.Vector m, N.Vector n] a
h' = N.Prism $ N.Prism $ N.Scalar $
     N.viota @n <&> (\(N.Fin x) ->
     N.viota @m <&> (\(N.Fin w) ->
     sin (pi * (fromIntegral w) * dx)
             * sin (2 * pi * (fromIntegral x) * dy)))

h34' :: Floating a => N.Hyper '[N.Vector 3, N.Vector 4] a
h34' = h'

safeUpdate' :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a, Show a) =>
                             N.Hyper '[N.Vector m, N.Vector n] a ->
                             N.Hyper '[N.Vector m, N.Vector n] a
safeUpdate' x = N.binary (+) s $
                N.foldrH (+) 0 $
                N.foldrH (+) 0 $
                N.binary (*) (N.binary (+) cc5 yy5) x
  where

    s = N.Prism $ N.Prism $ N.Scalar $
        N.viota @n <&> (\(N.Fin x) ->
        N.viota @m <&> (\(N.Fin w) ->
        sin (pi * (fromIntegral w) * dx)
                * sin (2 * pi * (fromIntegral x) * dy)))

    cc5 = N.binary (*) (N.Scalar c2) (fmap fromIntegral cc4)
    cc4 :: Num b => N.Hyper '[N.Vector m, N.Vector n, N.Vector m, N.Vector n] b
    cc4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
          N.viota @n <&> (\(N.Fin x) ->
          N.viota @m <&> (\(N.Fin w) ->
          N.viota @n <&> (\(N.Fin v) ->
          N.viota @m <&> (\(N.Fin u) ->
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

    yy5 = N.binary (*) (N.Scalar c1) (fmap fromIntegral yy4)
    yy4 = N.Prism $ N.Prism $ N.Prism $ N.Prism $ N.Scalar $
          N.viota @n <&> (\(N.Fin x) ->
          N.viota @m <&> (\(N.Fin w) ->
          N.viota @n <&> (\(N.Fin v) ->
          N.viota @m <&> (\(N.Fin u) ->
          (f m n x w v u)))))
            where
              m = fromIntegral $ M.natVal (undefined :: Proxy m)
              n = fromIntegral $ M.natVal (undefined :: Proxy n)
              f m n i j k l | i == 0                   = 0
                            | j == 0                   = 0
                            | i == n - 1               = 0
                            | j == m - 1               = 0
                            | k == i     && l == j - 1 = 1
                            | k == i     && l == j     = -2
                            | k == i     && l == j + 1 = 1
                            | otherwise                = 0

update' :: Int -> Int -> [Double] -> [Double]
update' k n xs =
  case M.someNatVal $ fromIntegral n of
    Nothing -> error "static/dynamic mismatch"
    Just (M.SomeNat (_ :: Proxy n)) ->
      case M.someNatVal $ fromIntegral k of
        Nothing -> error "static/dynamic mismatch"
        Just (M.SomeNat (_ :: Proxy k)) ->
          N.elements $ safeUpdate' b
          where
            b :: N.Hyper '[N.Vector k, N.Vector n] Double
            b = hyperfy xs

t0, tf :: Double
t0 = 0.0
tf = 0.3

bigNt :: Int
bigNt = 20

dTout :: Double
dTout = (tf - t0) / (fromIntegral bigNt)

ts = map (dTout *) $ map fromIntegral [1..bigNt]

initH :: forall m n a . (M.KnownNat m, M.KnownNat n, Floating a) =>
                     N.Hyper '[N.Vector m, N.Vector n] a
initH = N.Prism $ N.Prism $ N.Scalar $
        N.viota @n <&> (\(N.Fin x) ->
        N.viota @m <&> (\(N.Fin w) ->
        0.0))

sol :: Matrix Double
sol = odeSolveV SDIRK_5_3_4' Nothing 1.0e-5 1.0e-10 (const (fromList . (update' nx ny) . toList)) (assoc (nx * ny) 0.0 [] :: Vector Double) (fromList ts)

sol' :: Matrix Double
sol' = odeSolveV SDIRK_5_3_4' Nothing 1.0e-5 1.0e-10 (const bigU') (assoc (nx * ny) 0.0 [] :: Vector Double) (fromList $ ts)
  where
    bigU' bigU = bigA #> bigU + b

main = do
  mapM_ (\i -> putStrLn $ show $ sqrt $ (sol'!i) <.> (sol'!i) / (fromIntegral nx) / (fromIntegral ny)) ([0 .. length ts - 1] :: [Int])

-- https://serokell.io/blog/dimensions-haskell-singletons
