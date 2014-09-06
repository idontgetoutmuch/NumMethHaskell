{-# OPTIONS_GHC -Wall                     #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing  #-}
{-# OPTIONS_GHC -fno-warn-type-defaults   #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind  #-}
{-# OPTIONS_GHC -fno-warn-missing-methods #-}
{-# OPTIONS_GHC -fno-warn-orphans         #-}

{-# LANGUAGE NoMonomorphismRestriction    #-}
{-# LANGUAGE RankNTypes                   #-}
{-# LANGUAGE QuasiQuotes                  #-}
{-# LANGUAGE LambdaCase                   #-}
{-# LANGUAGE GADTs                        #-}
{-# LANGUAGE ViewPatterns                 #-}


module Main ( main ) where

import Resampling
import ResamplingChart

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine

import qualified Data.Vector.Unboxed as U


displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

diaNoisyData :: Diagram B R2
diaNoisyData =
  diagSamps (zip (map fromIntegral [0..]) (take nObs ys))
            (zip (map fromIntegral [0..]) (take nObs ys))
  where ys = samples

diaPartEvolve2 :: Diagram B R2
diaPartEvolve2 = diagPartDist2 "Distribution of Particles Over Time"
                 (blue `withOpacity` 0.2)
                 (red `withOpacity` 0.5)
                 (concat (zipWith zip
                          (map repeat [0..])
                          xss))
                 (concat (zipWith zip
                          (map repeat [0..])
                          xss))
  where
    xss = map U.toList $ map snd $ snd runPopSir

main :: IO ()
main = do
    displayHeader "diagrams/ResamplingPopGrowth.png" diaNoisyData
    putStrLn "Hello"