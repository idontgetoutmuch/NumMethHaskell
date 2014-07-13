{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

import KalmanChart
import Filtering

import Diagrams.Prelude hiding ( render, Renderable )

import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine


dia :: Diagram B R2
dia = diagPartFilter (zip [1..1000](snd runPF)) 3

displayHeader :: FilePath -> Diagram B R2 -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  displayHeader "diagrams/SingleRvNoisy.png" dia


