{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

{-# LANGUAGE NoMonomorphismRestriction     #-}
{-# LANGUAGE FlexibleContexts              #-}
{-# LANGUAGE TypeFamilies                  #-}

module Diagrams ( example
                , clockSelf
                , self
                ) where

import Diagrams.Prelude
import Diagrams.TwoD.Text

state :: Renderable (Path R2) b => Diagram b R2
state = circle 1 # lw 0.05 # fc silver

fState :: Renderable (Path R2) b => Diagram b R2
fState = circle 0.85 # lw 0.05 # fc lightblue <> state

clockPoints :: Int -> [(Double, Double)]
clockPoints n = [ (7 + 6*cos x, 6 + 6*sin x)
                | m <- [0..n-1]
                , let x = 2 * pi * fromIntegral m / fromIntegral n
                ]

clockTurns :: Integral a => a -> [Angle]
clockTurns n = [ x @@ turn
               | m <- [0..n-1]
               , let x = fromIntegral m / fromIntegral n
               ]

points :: Int -> [P2]
points n = map p2 $ map (\m -> (clockPoints n)!!m) [0..n-1]

ds :: (Renderable Text b, Renderable (Path R2) b) => Int -> [Diagram b R2]
ds m = map (\n -> (text (show (m - n + 1)) <> fState) # named (show n)) [1..m]

states :: (Renderable Text b, Renderable (Path R2) b) => Int -> Diagram b R2
states clockSize = position (zip (points clockSize) (ds clockSize))

line :: Trail R2
line = trailFromOffsets [unitX]

arrowStyle2 :: ArrowOpts
arrowStyle2  = (with  & arrowHead  .~ noHead &  tailSize .~ 0.3
                      & arrowShaft .~ shaft' & arrowTail .~ spike'
                      & shaftStyle %~ lw 0.1
                      & tailColor  .~ black)

shaft' :: (Transformable t, TrailLike t, V t ~ R2) => t
shaft' = arc (1/2 @@ turn) (0 @@ turn) # scaleX 0.33

arrowStyle3 :: ArrowOpts
arrowStyle3  = (with  & arrowHead  .~ noHead & tailSize .~ 0.3
                      & arrowShaft .~ line & arrowTail .~ spike'
                      & shaftStyle %~ lw 0.1
                      & tailColor  .~ black)

clockConnectWiddershins :: Renderable (Path R2) b =>
                           Int -> Int -> Diagram b R2 -> Diagram b R2
clockConnectWiddershins n m = connectPerim' arrowStyle3
                              (show (1 + ((n + 1) `mod` m))) (show (1 + (n `mod` m)))
                              ((clockTurns m)!!(n `mod` m))
                              ((clockTurns m)!!((n + 1) `mod` m))
clockConnectClockwise :: Renderable (Path R2) b =>
                         Int -> Int -> Diagram b R2 -> Diagram b R2
clockConnectClockwise n m = connectPerim' arrowStyle3
                            (show (1 + (n `mod` m))) (show (1 + ((n + 1) `mod` m)))
                            ((clockTurns m)!!((n + 2) `mod` m))
                            ((clockTurns m)!!((n - 1) `mod` m))

clockConnectSelf :: Renderable (Path R2) b =>
                    Int -> Int -> Diagram b R2 -> Diagram b R2
clockConnectSelf n m = connectPerim' arrowStyle2
                       (show (1 + (n `mod` m))) (show (1 + (n `mod` m)))
                       ((clockTurns m)!!((n - 1) `mod` m))
                       ((clockTurns m)!!((n + 1) `mod` m))

oneWay :: Renderable (Path R2) b =>
          Int -> Diagram b R2 -> Diagram b R2
oneWay clockSize = foldr (.) id
                   (map (flip clockConnectWiddershins clockSize) [0..clockSize - 1])

otherWay :: Renderable (Path R2) b =>
            Int -> Diagram b R2 -> Diagram b R2
otherWay clockSize = foldr (.) (oneWay clockSize)
                     (map (flip clockConnectClockwise clockSize) [0..clockSize - 1])

self :: Renderable (Path R2) b =>
        Int -> Diagram b R2 -> Diagram b R2
self clockSize = foldr (.) (otherWay clockSize)
                 (map (flip clockConnectSelf clockSize) [0..clockSize - 1])

example :: (Renderable Text b, Renderable (Path R2) b) =>
           Int -> Diagram b R2
example clockSize = (otherWay clockSize) $ states clockSize

clockSelf :: (Renderable Text b, Renderable (Path R2) b) =>
             Int -> Diagram b R2
clockSelf clockSize = (self clockSize) $ states clockSize

