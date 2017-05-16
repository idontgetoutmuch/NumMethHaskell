{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
{-# LANGUAGE TypeFamilies     #-}

import Prelude                            as P
import Data.Array.Accelerate              as A   hiding ((^))
import Data.Array.Accelerate.LLVM.Native  as CPU
import Data.Array.Accelerate.LLVM.PTX     as GPU
import Data.Array.Accelerate.Linear              hiding (trace)
import Data.Array.Accelerate.Control.Lens
import qualified Linear                   as L

e, q10, q20, p10, p20 :: Double
e = 0.6
q10 = 1 - e
q20 = 0.0
p10 = 0.0
p20 = sqrt ((1 + e) / (1 - e))

h :: Double
h = 0.01

oneStep2 :: Double -> Exp (V2 Double, V2 Double) -> Exp (V2 Double, V2 Double)
oneStep2 hh prev = lift (qNew, pNew)
  where
    h2 = hh / 2
    hhs = lift ((pure hh) :: V2 Double)
    hh2s = (lift ((pure h2) :: V2 Double))
    pp2 = psPrev - hh2s * nablaQ' qsPrev
    qNew = qsPrev + hhs * nablaP' pp2
    pNew = pp2 - hh2s * nablaQ' qNew
    qsPrev = A.fst prev
    psPrev = A.snd prev
    nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
    nablaQ' qs = lift (V2 (qq1 / r) (qq2 / r))
      where
        qq1 = qs ^. _x
        qq2 = qs ^. _y
        r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)
    nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
    nablaP' ps = ps

oneStepH98 :: Double -> V2 (V2 Double) -> V2 (V2 Double)
oneStepH98 hh prev = V2 qNew pNew
  where
    h2 = hh / 2
    hhs = V2 hh hh
    hh2s = V2 h2 h2
    pp2 = psPrev - hh2s * nablaQ' qsPrev
    qNew = qsPrev + hhs * nablaP' pp2
    pNew = pp2 - hh2s * nablaQ' qNew
    qsPrev = prev ^. L._x
    psPrev = prev ^. L._y
    nablaQ' qs = V2 (qq1 / r) (qq2 / r)
      where
        qq1 = qs ^. L._x
        qq2 = qs ^. L._y
        r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)
    nablaP' ps = ps

-- oneStep :: Double -> Exp (V2 (V2 Double)) -> Exp (V2 (V2 Double))
-- oneStep h prev = lift $ V2 qNew pNew
--   where
--     hs  = lift (V2 h h)
--     h2s = lift (V2 h2 h2) where h2 = h / 2
--     qsPrev = prev ^. _x
--     psPrev = prev ^. _y
--     nablaQ qs = lift (V2 (qq1 / r) (qq2 / r))
--       where
--         qq1 = qs ^. _x
--         qq2 = qs ^. _y
--         r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)
--     nablaP ps = ps
--     p2 = psPrev - h2s * nablaQ qsPrev
--     qNew = qsPrev + hs * nablaP p2
--     pNew = p2 - h2s * nablaQ qNew

-- nSteps :: Int
-- nSteps = 100

dummyStart :: Exp (V2 Double, V2 Double)
dummyStart = lift (V2 q10 q20, V2 p10 p20)

dummyStart98 :: V2 (V2 Double)
dummyStart98 = V2 (V2 q10 q20) (V2 p10 p20)

-- dummyInputs :: Acc (Array DIM1 (V2 Double, V2 Double))
-- dummyInputs = A.use $ A.fromList (Z :. nSteps) $
--               P.replicate nSteps (pure 0.0 :: V2 Double, pure 0.0 :: V2 Double)

-- runSteps :: Acc (Array DIM1 (V2 Double, V2 Double))
-- runSteps = A.scanl (\s _x -> (oneStep2 h s)) dummyStart dummyInputs

nSteps :: Int
nSteps = 80000000 -- 100000000

runSteps' :: Exp (V2 Double, V2 Double) -> Exp (V2 Double, V2 Double)
runSteps' = A.iterate (lift nSteps) (oneStep2 h)

myIterate :: Int -> (a -> a) -> a -> a
myIterate n f x | n P.<= 0 = x
                | otherwise = myIterate (n-1) f $! f x
-- myIterate 0 _ a = a
-- myIterate n f a = myIterate (n-1) f (f a)
-- myIterate n f a = P.last $ P.take n $ P.iterate f a

runSteps98' :: V2 (V2 Double) -> V2 (V2 Double)
runSteps98' = myIterate nSteps (oneStepH98 h)

reallyRunSteps' :: (Array DIM1 (V2 Double, V2 Double))
reallyRunSteps' = CPU.run $
                  A.scanl (\s _x -> runSteps' s) dummyStart
                  (A.use $ A.fromList (Z :. 1) [(V2 0.0 0.0, V2 0.0 0.0)])

main :: IO ()
main = do
  putStrLn $ show $ reallyRunSteps'
  -- putStrLn $ show $ runSteps98' dummyStart98

-- Onesteph98 :: Double -> (V2 Double, V2 Double) -> (V2 Double, V2 Double)
-- oneStepH98 h prev = (qNew, pNew)
--   where
--     h2 = h / 2
--     hhs = pure h
--     hs = pure h2
--     p2 = psPrev - hs * nablaQ qsPrev
--     qNew = qsPrev + hhs * nablaP p2
--     pNew = p2 - hs * nablaQ qNew
--     qsPrev = P.fst prev
--     psPrev = P.snd prev
--     nablaQ qs = V2 (q1 / r) (q2 / r)
--       where
--         q1 = qs ^. L._x
--         q2 = qs ^. L._y
--         r   = (q1 ^ 2 + q2 ^ 2) ** (3/2)
--     nablaP ps = ps

-- bigH2BodyH98 :: (V2 Double, V2 Double) -> Double
-- bigH2BodyH98 x = ke + pe
--   where
--     pe = let V2 q1 q2 = P.fst x in negate $ recip (sqrt (q1^2 + q2^2))
--     ke = let V2 p1 p2 = P.snd x in 0.5 * (p1^2 + p2^2)
