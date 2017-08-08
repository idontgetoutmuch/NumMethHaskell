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

import Data.Maybe (fromJust, fromMaybe)

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

test2 :: ((Matrix Double, Matrix Double, [Vector Double]), Matrix Double)
test2 = S.withMatrix bigR f2
  where
    f2 :: forall n' k . (KnownNat n', KnownNat k) => S.L n' k -> ((Matrix Double, Matrix Double, [Vector Double]), Matrix Double)
    f2 r = ( S.withMatrix (fromLists [map eruptions oldFaithful, map waiting oldFaithful]) g2
          , S.extract nK2
          )
      where
        g2 :: forall d n . (KnownNat d, KnownNat n) => S.L d n -> (Matrix Double, Matrix Double, [Vector Double])
        g2 x = case sameNat (Proxy :: Proxy n') (Proxy :: Proxy n) of
                Nothing -> error "Incompatible matrices"
                Just Refl -> (S.extract xbars, S.extract ((S.tr x) - (mkXbarss (xxbars!!0))), map S.extract xxbars)
                               where
                                 rs :: [S.L d n] -- k
                                 rs = map (repmat' (Proxy :: Proxy d)) $
                                       map (S.row . fromJust . S.create) $
                                       toRows $ S.extract $ S.tr r
                                 xxbars :: [S.R d] -- k
                                 xxbars = map sumRows $ rs .* (repeat x)
                                 xxbars' :: [S.R d] -- k
                                 xxbars' = zipWith (\n u -> S.dvmap (/n) u)
                                                   (toList (S.extract $ S.unrow nK2))
                                                   (map sumRows $ rs .* (repeat x))
                                 xbars :: S.L k d
                                 xbars = (S.tr r) S.<> (S.tr x)
                                 mkXbarss :: S.R d -> S.L n d
                                 mkXbarss u = fromJust $
                                              S.create $
                                              fromRows $
                                              replicate l (S.extract u)
                                 _foo :: [S.L n d] -- k
                                 _foo = map mkXbarss xxbars
                                 _bar :: [S.L n d] -- k
                                 _bar = map (S.tr x -) _foo
                                 _baz :: [[S.R d]] -- k x n?
                                 _baz = map (map (fromJust . S.create) . toRows . S.extract) _bar
                                 _urk :: [S.L d d] -- k?
                                 _urk = map sum $ map (map (\x -> (S.col x) S.<> (S.row x))) _baz
        l = fst $ S.size r
        nK2 :: S.L 1 k
        nK2 = S.row (S.vector $ take l ones) S.<> r

test3 :: ([Vector Double], Matrix Double)
test3 = S.withMatrix bigR f3
  where
    f3 :: forall n' k . (KnownNat n', KnownNat k) =>
                       S.L n' k -> ([Vector Double], Matrix Double)
    f3 r = ( S.withMatrix (fromLists [map eruptions oldFaithful, map waiting oldFaithful]) g3
          , S.extract nK3
          )
      where
        g3 :: forall d n . (KnownNat d, KnownNat n) =>
                          S.L d n -> [Vector Double]
        g3 x = case sameNat (Proxy :: Proxy n') (Proxy :: Proxy n) of
                Nothing -> error "Incompatible matrices"
                Just Refl -> map S.extract xbars3
                               where
                                 rs :: [S.L d n] -- k
                                 rs = map (repmat' (Proxy :: Proxy d)) $
                                       map (S.row . fromJust . S.create) $
                                       toRows $ S.extract $ S.tr r
                                 xbars3 :: [S.R d] -- k
                                 xbars3 = zipWith (\n u -> S.dvmap (/n) u)
                                                 (toList (S.extract $ S.unrow nK3))
                                                 (map sumRows $ rs .* (repeat x))
        l = fst $ S.size r
        nK3 :: S.L 1 k
        nK3 = S.row (S.vector $ take l ones) S.<> r

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

withMatrix
    :: forall a z
     . [[a]]
    -> (forall k n . (KnownNat k, KnownNat n) => N.Hyper '[N.Vector k, N.Vector n] a -> z)
    -> z
withMatrix a f =
    case someNatVal $ fromIntegral $ length a of
       Nothing -> error "static/dynamic mismatch"
       Just (SomeNat (_ :: Proxy n)) ->
           case someNatVal $ fromIntegral $ (length . head) a of
               Nothing -> error "static/dynamic mismatch"
               Just (SomeNat (_ :: Proxy k)) ->
                 f b
                 where
                    b :: N.Hyper '[N.Vector k, N.Vector n] a
                    b = N.Prism $ N.Prism $ N.Scalar $
                        (fromJust' . N.fromList) $ map (fromJust' . N.fromList) a
                    fromJust' = fromMaybe (error "static/dynamic mismatch")

sumMatrix :: (Num a, KnownNat k, KnownNat n) => N.Hyper '[N.Vector k, N.Vector n] a -> a
sumMatrix = N.point . N.foldrH (+) 0 . N.foldrH (+) 0

-- withCuboid
--     :: forall a z
--      . [[[a]]]
--     -> (forall k n . (KnownNat d, KnownNat k, KnownNat n) =>
--         N.Hyper '[N.Vector n, N.Vector k, N.Vector d] a -> z)
--     -> z
-- withCuboid a f =
--     case someNatVal $ fromIntegral $ length a of
--        Nothing -> error "static/dynamic mismatch"
--        Just (SomeNat (_ :: Proxy k)) ->
--            case someNatVal $ fromIntegral $ (length . head) a of
--                Nothing -> error "static/dynamic mismatch"
--                Just (SomeNat (_ :: Proxy n)) ->
--                  case someNatVal $ fromIntegral $ (length . head . head) a of
--                    Nothing -> error "static/dynamic mismatch"
--                    Just (SomeNat (_ :: Proxy d)) ->
--                      f b
--                      where
--                        b :: N.Hyper '[N.Vector k, N.Vector n] a
--                        b = N.Prism $ N.Prism $ N.Scalar $
--                            (fromJust . N.fromList) $
--                            map (fromJust . N.fromList) a

bigXH :: N.Matrix 272 2 Double
bigXH = (fromJust . N.fromList) $
        map (fromJust . N.fromList) $
        map toList $
        toRows bigX

f :: (KnownNat n, KnownNat k) => N.Matrix n k Double -> N.Hyper '[N.Vector k] Double
f x = N.foldrH (+) 0 $
      N.transposeH (N.Prism (N.Prism (N.Scalar x)))

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

g :: forall d k n . (KnownNat d, KnownNat k, KnownNat n) =>
     N.Hyper '[N.Vector k, N.Vector n] Double ->
     N.Hyper '[N.Vector n, N.Vector k, N.Vector d] Double ->
     ( N.Hyper '[N.Vector k] Double
     , N.Hyper '[N.Vector k, N.Vector d] Double
     , N.Hyper '[N.Vector k, N.Vector d, N.Vector d] Double
     )
g biggR biggXs = (nnK, xxbars, biggS)
        where

          nnK = sums $ N.transposeH biggR

          xxbars = xsums ./. nnK
          xsums = sums $ (N.transposeH biggR) .*. biggXs

          difff1 = biggXs .-. b
          b :: N.Hyper '[N.Vector n, N.Vector k, N.Vector d] Double
          b = N.transposeH $ N.transposeH' $ focus N.replicate xxbars

          difff2 = c .*. difff1
          c :: N.Hyper '[N.Vector n, N.Vector k, N.Vector d] Double
          c = N.transposeH $ focus N.replicate biggR

          biggS = biggS' ./. nnK
          biggS' = N.transposeH $ N.transposeH' $
                   N.Prism $ N.Prism $
                   N.binary N.matrix d2 d1'
          d1' = N.unary N.transpose $ N.crystal $ N.crystal $ N.transposeH' difff1
          d2 = N.crystal $ N.crystal $ N.transposeH' difff2

          focus h = N.Prism . N.Prism . N.Prism . N.Scalar .
                    h . N.point . N.crystal . N.crystal

          sums = N.foldrH (+) 0

infixl 6 .+.

(.+.) :: (Num c, N.Alignable gs (N.Max fs gs),
          N.Alignable fs (N.Max fs gs), N.Compatible fs gs) =>
         N.Hyper fs c -> N.Hyper gs c -> N.Hyper (N.Max fs gs) c
(.+.) = N.binary (+)

infixl 6 .-.

(.-.) :: (Num c, N.Alignable gs (N.Max fs gs),
          N.Alignable fs (N.Max fs gs), N.Compatible fs gs) =>
         N.Hyper fs c -> N.Hyper gs c -> N.Hyper (N.Max fs gs) c
(.-.) = N.binary (-)

infixl 7 ./.

(./.) :: (Fractional c, N.Alignable gs (N.Max fs gs),
          N.Alignable fs (N.Max fs gs), N.Compatible fs gs) =>
         N.Hyper fs c -> N.Hyper gs c -> N.Hyper (N.Max fs gs) c
(./.) = N.binary (/)

infixl 7 .*.

(.*.) :: (Fractional c, N.Alignable gs (N.Max fs gs),
          N.Alignable fs (N.Max fs gs), N.Compatible fs gs) =>
         N.Hyper fs c -> N.Hyper gs c -> N.Hyper (N.Max fs gs) c
(.*.) = N.binary (*)

main :: IO ()
main = do
  let ((xbar, ms1, xxs), nnnK) = test2
      vk = cmap (v0 +) nnnK
      _betak = cmap (beta0 +) nnnK
  putStrLn $ show $ size nnnK
  -- putStrLn $ show ms1
  putStrLn $ show nnnK
  -- putStrLn $ show vk
  putStrLn $ show xbar
  -- putStrLn $ show xxs
  let (n,m,s) = g bigRH' (N.Prism $ N.Prism $ N.Prism $ N.Scalar bigXHs)
  putStrLn $ show n
  putStrLn $ show m
  putStrLn "Hello, Haskell!"
