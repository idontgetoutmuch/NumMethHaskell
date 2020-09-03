{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE OverloadedLists     #-}
{-# LANGUAGE NumDecimals         #-}
{-# LANGUAGE TypeApplications    #-}
{-# LANGUAGE TypeFamilies        #-}
{-# LANGUAGE DataKinds           #-}
{-# LANGUAGE FlexibleContexts    #-}
{-# LANGUAGE OverloadedStrings   #-}

{-# OPTIONS_GHC -Wall            #-}

module Main (main) where

import           Numeric.Sundials
import           Numeric.LinearAlgebra
import           Data.Csv
import           Data.Char
import           Data.ByteString.Lazy (putStr, writeFile)
import           Prelude hiding (putStr, writeFile)
import           Control.Exception
import           Data.Coerce
import           Katip
import           System.IO (stdout)
import           Katip.Monadic
import           GHC.Int


bigN :: Int
-- bigN = 101
bigN = 200

deltaX :: Double
deltaX = 1.0 / (fromIntegral bigN - 1)

deltaX' :: Double
deltaX' = 1.0 / (fromIntegral bigN)

bigC :: Matrix Double
bigC = assoc (bigN, bigN) 0.0 [ ((i, j), (1/(2*deltaX)) * f (i, j)) | i <- [0 .. bigN - 1]
                                                                    , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j)
            | i     == 0        =  0
            | i     == bigN - 1 =  0
            | i     == j        =  0
            | i - 1 == j        = -1
            | i + 1 == j        =  1
            | otherwise         =  0

bigD :: Matrix Double
bigD = assoc (bigN, bigN) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. bigN - 1]
                                                   , j <- [0 .. bigN - 1]
                        ]
 where
   f (i, j) | i     == j        =  1
            | otherwise         =  0

bigV :: Vector Double
bigV = assoc bigN 0.0 [(i, f i) | i <- [0 .. bigN - 1]]
  where
  f i = (sin (2 * pi * fromIntegral i * deltaX))

bigV' :: Vector Double
bigV' = assoc (bigN + 1) 0.0 [(i, f i) | i <- [0 .. bigN]]
  where
  f i = (sin (2 * pi * fromIntegral i * deltaX'))

sol' :: IO (Matrix Double)
sol' = do
  w <- runNoLoggingT $ solve (defaultOpts' HEUN_EULER_2_1_2) burgers
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

sol'' :: IO (Matrix Double)
sol'' = do
  -- w <- runNoLoggingT $ solve (defaultOpts' HEUN_EULER_2_1_2) burgersWeno
  handleScribe <- mkHandleScribe ColorIfTerminal stdout (permitItem InfoS) V2
  logEnv <- registerScribe "stdout" handleScribe defaultScribeSettings =<< initLogEnv "namespace" "devel"
  w <- runKatipT logEnv $ solve (defaultOpts' BOGACKI_SHAMPINE_4_2_3) burgersWeno
  case w of
    Left e  -> error $ show e
    Right y -> return (solutionMatrix y)

burgers :: OdeProblem
burgers = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x -> coerce (cmap negate ((bigD #> (coerce x)) * (bigC #> (coerce x))))
  , odeJacobian = Nothing
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = bigV
  , odeSolTimes = vector $ map (0.025 *) [0 .. 10]
  , odeTolerances = defaultTolerances
  }

burgersWeno :: OdeProblem
burgersWeno = emptyOdeProblem
  { odeRhs = odeRhsPure $ \_t x -> coerce (rhs' bigN (coerce x))
  , odeJacobian = Nothing
  , odeEvents = mempty
  , odeEventHandler = nilEventHandler
  , odeMaxEvents = 0
  , odeInitCond = bigV'
  , odeSolTimes = vector $ map (0.025 *) [0 .. 10]
  , odeTolerances = defaultTolerances
  }

wcL :: Double -> Double -> Double -> Double -> Double -> Double
wcL v1 v2 v3 v4 v5 = f
  where
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^(2 :: Int) + 0.25*(v1-4.0*v2+3.0*v3)^(2 :: Int)
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^(2 :: Int) + 0.25*(v2-v4)^(2 :: Int)
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^(2 :: Int) + 0.25*(3.0*v3-4.0*v4+v5)^(2 :: Int)

    c1 = 2.0e-1/((eps+s1)^(2 :: Int))
    c2 = 5.0e-1/((eps+s2)^(2 :: Int))
    c3 = 3.0e-1/((eps+s3)^(2 :: Int))

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    q1 = v1/3.0   - 7.0/6.0*v2 + 11.0/6.0*v3
    q2 =(-v2)/6.0 + 5.0/6.0*v3 + v4/3.0
    q3 = v3/3.0   + 5.0/6.0*v4 - v5/6.0

    f = (w1*q1 + w2*q2 + w3*q3)

wcR :: Double -> Double -> Double -> Double -> Double -> Double
wcR v1 v2 v3 v4 v5 = f
  where
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^(2 :: Int) + 0.25*(v1-4.0*v2+3.0*v3)^(2 :: Int)
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^(2 :: Int) + 0.25*(v2-v4)^(2 :: Int)
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^(2 :: Int) + 0.25*(3.0*v3-4.0*v4+v5)^(2 :: Int)

    c1 = 3.0e-1/(eps+s1)^(2 :: Int)
    c2 = 5.0e-1/(eps+s2)^(2 :: Int)
    c3 = 2.0e-1/(eps+s3)^(2 :: Int)

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    q1 =(-v1)/6.0    + 5.0/6.0*v2 + v3/3.0
    q2 = v2/3.0      + 5.0/6.0*v3 - v4/6.0
    q3 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0

    f = (w1*q1 + w2*q2 + w3*q3)

crwenoR :: Int -> Vector Double -> Vector Double
crwenoR n u = assoc (n + 1) 0.0 [(i, f i) | i <- [0 .. n]]
  where
    f i | i == 0    = 0.0

    f i | i == 1    = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(n-1)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == n - 1 = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!1

    f i | i == n     = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!1
        v5 = u!2

    f i | otherwise  = wcR v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

crwenoL :: Int -> Vector Double -> Vector Double
crwenoL n u = assoc (n + 1) 0.0 [(i, f i) | i <- [0 .. n]]
  where
    f i | i == 0 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(n-2)
        v2 = u!(n-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == 1 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(n-1)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

    f i | i == n - 1 = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!1

    f i | i == n     = 0.0

    f i | otherwise  = wcL v1 v2 v3 v4 v5
      where
        v1 = u!(i-2)
        v2 = u!(i-1)
        v3 = u!i
        v4 = u!(i+1)
        v5 = u!(i+2)

mLR :: Int -> Matrix Double
mLR n = assoc (n, n) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                              , j <- [0 .. n - 1]
                                              ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  1
             | i == 0 &&
               j == n - 1 = -1
             | i == j     =  1
             | i - 1 == j = -1
             | otherwise  =  0

mL :: Int -> Matrix Double
mL n = assoc (n, n + 1) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                                 , j <- [0 .. n]
                                                 ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  1
             | i == 0 &&
               j == n - 1 = -1
             | i == j     =  1
             | i - 1 == j = -1
             | otherwise  =  0

mR :: Int -> Matrix Double
mR n = assoc (n, n + 1) 0.0 [ ((i, j), f (i, j)) | i <- [0 .. n - 1]
                                                 , j <- [0 .. n]
                                                 ]
  where
    f (i, j) | i == 0 &&
               j == 0     =  0
             | i == 0 &&
               j == 1     =  1
             | i == 0 &&
               j == n     = -1
             | i == j     = -1
             | i + 1 == j =  1
             | otherwise  =  0

preRhs :: Int -> Vector Double -> Vector Double
preRhs n v = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    ll = (mL n) #> (crwenoL n v)
    rr = (mR n) #> (crwenoR n v)

    f i | v!i >= 0.0 = negate (v!i * ll!i) / deltaX'
        | otherwise  = negate (v!i * rr!i) / deltaX'

rhs' :: Int -> Vector Double -> Vector Double
rhs' n v = vjoin [preRhs n v, [v!0]]

rhs'' :: Int -> Vector Double -> Vector Double
rhs'' n v = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    ll = (mL n) #> (crwenoL n (vjoin [v, vector [0.0]]))
    rr = (mR n) #> (crwenoR n (vjoin [vector [0.0], v]))

    f i | v!i >= 0.0 = negate (v!i * ll!i) / deltaX
        | otherwise  = negate (v!i * rr!i) / deltaX

rhs :: Int -> Vector Double -> Vector Double
rhs n u = assoc n 0.0 [(i, f i) | i <- [0 .. n - 1]]
  where
    uL = crwenoL n u
    uR = crwenoR n u
    ll = (mLR n) #> (subVector 0 n uL)
    rr = (mLR n) #> (subVector 1 n uR)
    f i | u!i >= 0.0 = negate (u!i * ll!i)
        | otherwise  = negate (u!i * rr!i)

myOptions :: EncodeOptions
myOptions = defaultEncodeOptions {
      encDelimiter = fromIntegral (ord ' ')
    }

main :: IO ()
main = do
  -- y <- sol'
  -- writeFile "burgers.txt" $ encodeWith myOptions $ map toList $ toRows y
  y <- sol''
  writeFile "burgersWeno.txt" $ encodeWith myOptions $ map toList $ toRows y


defaultOpts' :: method -> ODEOpts method
defaultOpts' method = ODEOpts
  { maxNumSteps = 1e5
  , minStep     = 1.0e-14
  , fixedStep   = 0.0 -- 0.001
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
      , odeEvents = mempty
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
  { absTolerances = Left 1.0e-6
  , relTolerance = 1.0e-10
  }
