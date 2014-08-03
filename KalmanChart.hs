{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module KalmanChart (
    diagPartFilter
  , diagEsts
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

diagPartFilter :: Double -> [(Double, Double)] -> Int -> QDiagram Cairo R2 Any
diagPartFilter x ls n =
  fst $ runBackend denv (render (chartPartFilter x ls n) (500, 500))

chartPartFilter :: Double -> [(Double, Double)] -> Int -> Renderable ()
chartPartFilter x lineVals n = toRenderable layout
  where

    fitted = plot_lines_values .~ [lineVals]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ "Trajectory"
              $ def

    actual = plot_lines_values .~ [zip (map fst lineVals) (repeat x)]
              $ plot_lines_style  . line_color .~ opaque red
              $ plot_lines_title .~ "Actual"
              $ def

    layout = layout_title .~ "Particle Filtering"
           $ layout_y_axis . laxis_generate .~ scaledAxis def (-3,3)
           $ layout_x_axis . laxis_generate .~ scaledAxis def (0,(fromIntegral n))

           $ layout_plots .~ [ toPlot fitted
                             , toPlot actual
                             ]
           $ def

diagEsts :: [(Double, Double)] ->
            [(Double, Double)] ->
            [(Double, Double)] ->
            [(Double, Double)] ->
            QDiagram Cairo R2 Any
diagEsts ls us ks ps =
  fst $ runBackend denv (render (chartEsts ls us ks ps) (500, 500))

chartEsts :: [(Double, Double)] ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             [(Double, Double)] ->
             Renderable ()
chartEsts lineVals upperVals lowerVals pointVals = toRenderable layout
  where

    fitted = plot_lines_values .~ [lineVals]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_title .~ "Estimate"
              $ def

    upper  = plot_lines_values .~ [upperVals]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ "Upper error"
              $ def

    lower  = plot_lines_values .~ [lowerVals]
              $ plot_lines_style  . line_color .~ opaque green
              $ plot_lines_title .~ "Lower error"
              $ def

    dataPts = plot_points_style .~ filledCircles 2 (opaque red)
              $ plot_points_values .~ pointVals
              $ plot_points_title .~ "Samples"
              $ def

    layout = layout_title .~ "Estimation by Kalman Filtering"
           $ layout_plots .~ [toPlot fitted,
                              toPlot upper,
                              toPlot lower,
                              toPlot dataPts]
           $ def
