{-# LANGUAGE NoMonomorphismRestriction #-}

module Diagrams where

import Diagrams.Prelude
import Diagrams.TwoD.Arrow

import Data.Maybe (fromMaybe)

state = circle 1 # lw 0.05 # fc silver
fState = circle 0.85 # lw 0.05 # fc lightblue <> state

clockPoints :: Int -> [(Double, Double)]
clockPoints n = [ (7 + 6*cos x, 6 + 6*sin x)
                | m <- [0..n-1]
                , let x = 2 * pi * fromIntegral m / fromIntegral n
                ]

clockTurns n = [ x @@ turn
               | m <- [0..n-1]
               , let x = fromIntegral m / fromIntegral n
               ]

points n = map p2 $ map (\m -> (clockPoints n)!!m) [0..n-1]

ds m = map (\n -> (text (show n) <> fState) # named (show n)) [1..m]

label txt size = text txt # fontSize size

clockSize = 12

states = position (zip (points clockSize) (ds clockSize))

shaft = arc (0 @@ turn) (1/6 @@ turn)
shaft' = arc (0 @@ turn) (1/2 @@ turn) # scaleX 0.33
line = trailFromOffsets [unitX]

arrowStyle1 = (with  & arrowHead  .~ noHead & tailSize .~ 0.3
                     & arrowShaft .~ shaft  & arrowTail .~ spike'
                     & tailColor  .~ black)

arrowStyle2  = (with  & arrowHead  .~ noHead &  tailSize .~ 0.3
                      & arrowShaft .~ shaft' & arrowTail .~ spike'
                      & tailColor  .~ black)

arrowStyle3  = (with  & arrowHead  .~ noHead & tailSize .~ 0.3
                      & arrowShaft .~ line & arrowTail .~ spike'
                      & tailColor  .~ black)

clockConnectWiddershins n m = connectPerim' arrowStyle3
                              (show (1 + ((n + 1) `mod` m))) (show (1 + (n `mod` m)))
                              ((clockTurns m)!!(n `mod` m))
                              ((clockTurns m)!!((n + 1) `mod` m))

clockConnectClockwise n m = connectPerim' arrowStyle3
                            (show (1 + (n `mod` m))) (show (1 + ((n + 1) `mod` m)))
                            ((clockTurns m)!!((n + 2) `mod` m))
                            ((clockTurns m)!!((n - 1) `mod` m))

oneWay   = foldr (.) id     (map (flip clockConnectWiddershins clockSize) [0..clockSize - 1])
otherWay = foldr (.) oneWay (map (flip clockConnectClockwise clockSize)   [0..clockSize - 1])

example = otherWay states

