{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable, (*~), Time )
import Diagrams.Backend.CmdLine

import System.IO.Unsafe

import Symplectic

import qualified Data.Array.Accelerate as A   hiding ((^))
import qualified Linear                as L

import Text.Printf


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

chart :: String ->
         String ->
         [[(Double, Double)]] ->
         Renderable ()
chart t l obss = toRenderable layout
  where

    actual x l c = plot_lines_values .~ [x]
                   $ plot_lines_style  . line_color .~ opaque c
                   -- $ plot_lines_title .~ l
                   $ plot_lines_style  . line_width .~ 1.0
                   $ def

    ls = map (\n -> "Path " ++ show n) [1..]
    cs = cycle [blue, green, red, brown, crimson]

    actuals' :: [PlotLines Double Double]
    actuals' = zipWith3 actual obss ls cs

    layout = layout_title .~ t
           $ layout_plots .~ (map toPlot actuals')
           $ layout_y_axis . laxis_title .~ l
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Positions"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def


diagrm :: String -> String -> [[(Double, Double)]] -> Diagram Cairo
diagrm t l xss = fst $ runBackendR denv (chart t l xss)

diagrmM :: String -> String -> [[(Double, Double)]] -> IO (Diagram Cairo)
diagrmM t l xss = do
  denv <- defaultEnv vectorAlignmentFns 600 500
  return $ fst $ runBackendR denv (chart t l xss)

main = do
  let pqs  = A.toList runSteps
      pq1s = map (^. L._x) pqs
      pq2s = map (^. L._y) pqs
      p1s  = map (^. L._x) pq1s
      q1s  = map (^. L._y) pq1s
      p2s  = map (^. L._x) pq2s
      q2s  = map (^. L._y) pq2s
  d <- diagrmM "Orbits" "Momenta" [zip q1s p1s, zip q2s p2s]
  displayHeader "diagrams/symplectic.png" d


