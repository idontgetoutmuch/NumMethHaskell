{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators    #-}
{-# LANGUAGE TypeFamilies     #-}

import Prelude                            as P
import Data.Array.Accelerate              as A   hiding ((^))
import Data.Array.Accelerate.LLVM.Native  as CPU
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
oneStep2 h' prev = lift (qNew, pNew)
  where
    h2 = h' / 2
    gs = lift ((pure h') :: V2 Double)
    hs' = (lift ((pure h2) :: V2 Double))
    p2' = psPrev - hs' * nablaQ' qsPrev
    qNew = qsPrev + gs * nablaP' p2'
    pNew = p2' - hs' * nablaQ' qNew
    qsPrev = A.fst prev
    psPrev = A.snd prev
    nablaQ' :: Exp (V2 Double) -> Exp (V2 Double)
    nablaQ' qs = lift (V2 (q1' / r) (q2' / r))
      where
        q1' = qs ^. _x
        q2' = qs ^. _y
        r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
    nablaP' :: Exp (V2 Double) -> Exp (V2 Double)
    nablaP' ps = ps

-- oneStep :: Double -> Exp (V2 (V2 Double)) -> Exp (V2 (V2 Double))
-- oneStep h prev = lift $ V2 qNew pNew
--   where
--     hs  = lift (V2 h h)
--     h2s = lift (V2 h2 h2) where h2 = h / 2
--     qsPrev = prev ^. _x
--     psPrev = prev ^. _y
--     nablaQ qs = lift (V2 (q1' / r) (q2' / r))
--       where
--         q1' = qs ^. _x
--         q2' = qs ^. _y
--         r   = (q1' ^ 2 + q2' ^ 2) ** (3/2)
--     nablaP ps = ps
--     p2 = psPrev - h2s * nablaQ qsPrev
--     qNew = qsPrev + hs * nablaP p2
--     pNew = p2 - h2s * nablaQ qNew

-- nSteps :: Int
-- nSteps = 100

dummyStart :: Exp (V2 Double, V2 Double)
dummyStart = lift (V2 q10 q20, V2 p10 p20)

-- dummyInputs :: Acc (Array DIM1 (V2 Double, V2 Double))
-- dummyInputs = A.use $ A.fromList (Z :. nSteps) $
--               P.replicate nSteps (pure 0.0 :: V2 Double, pure 0.0 :: V2 Double)

-- runSteps :: Acc (Array DIM1 (V2 Double, V2 Double))
-- runSteps = A.scanl (\s _x -> (oneStep2 h s)) dummyStart dummyInputs

runSteps' :: Exp (V2 Double, V2 Double) -> Exp (V2 Double, V2 Double)
runSteps' = A.iterate (lift (100000000 :: Int)) (oneStep2 h)

-- reallyRunSteps :: (Array DIM1 (V2 Double, V2 Double))
-- reallyRunSteps = run runSteps

reallyRunSteps' :: (Array DIM1 (V2 Double, V2 Double))
reallyRunSteps' = run $
                  A.scanl (\s _x -> runSteps' s) dummyStart
                  (A.use $ A.fromList (Z :. 1) [(V2 0.0 0.0, V2 0.0 0.0)])

main :: IO ()
main = do
  putStrLn $ show $ reallyRunSteps'

-- oneStepH98 :: Double -> (V2 Double, V2 Double) -> (V2 Double, V2 Double)
-- oneStepH98 h prev = (qNew, pNew)
--   where
--     h2 = h / 2
--     gs = pure h
--     hs = pure h2
--     p2 = psPrev - hs * nablaQ qsPrev
--     qNew = qsPrev + gs * nablaP p2
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
