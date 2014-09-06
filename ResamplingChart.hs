{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module ResamplingChart (
    diagSamps
  , diagPartDist2
  ) where

import Control.Lens hiding ( (#) )
import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Diagrams.Backend.Cairo.CmdLine
import Diagrams.Prelude hiding ( render, Renderable )
import Data.Default.Class

import System.IO.Unsafe

denv :: DEnv
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 500 500

diagSamps :: [(Double, Double)] ->
             [(Double, Double)] ->
             QDiagram Cairo R2 Any
diagSamps xsActual xsEstimated =
  fst $ runBackend denv (render (chartSamps xsActual xsEstimated) (500, 500))

chartSamps :: [(Double, Double)] ->
             [(Double, Double)] ->
             Renderable ()
chartSamps as es = toRenderable layout
  where

    actual = plot_lines_values .~ [as]
             $ plot_lines_style  . line_color .~ opaque blue
             $ plot_lines_title .~ "Actual"
             $ def

    est   = plot_lines_values .~ [es]
            $ plot_lines_style  . line_color .~ opaque green
            $ plot_lines_title .~ "Sampled"
            $ def

    layout = layout_title .~ "Noisy Population Growth"
           $ layout_plots .~ [ toPlot actual
                             , toPlot est
                             ]
           $ def

chartPartDist2 :: String ->
                  AlphaColour Double -> AlphaColour Double ->
                  [(Double, Double)] -> [(Double, Double)] ->
                  Graphics.Rendering.Chart.Renderable ()
chartPartDist2 t c1 c2 prices1 prices2 = toRenderable layout
  where

    price1 = plot_points_style . point_color .~ c1
           $ plot_points_values .~ prices1
           $ plot_points_title .~ "Sample Particles Certain"
           $ def

    price2 = plot_points_style .~ filledCircles 1 c2
           $ plot_points_values .~ prices2
           $ plot_points_title .~ "Sample Particles Doubtful"
           $ def

    layout = layout_title .~ t
           $ layout_y_axis . laxis_title .~ "Particle Value"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ layout_plots .~ [ toPlot price1
                             , toPlot price2
                             ]
           $ layout_grid_last .~ False
           $ def

diagPartDist2 :: String ->
                 AlphaColour Double -> AlphaColour Double ->
                 [(Double, Double)] -> [(Double, Double)] ->
                 QDiagram Cairo R2 Any
diagPartDist2 t c1 c2 prices1 prices2 =
  fst $ runBackend denv (render (chartPartDist2 t c1 c2 prices1 prices2) (500, 500))
