{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module BrownianChart (
    diag
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

diag :: Double -> String ->
        [(Double, Double)] ->
        Diagram Cairo R2
diag yHt t ls =
  fst $ runBackend denv (render (chart yHt t ls) (500, 500))

chart :: Double -> String ->
         [(Double, Double)] ->
         Renderable ()
chart yHt t lineVals = toRenderable layout
  where

    actuals = plot_lines_values .~ [lineVals]
              $ plot_lines_style  . line_color .~ opaque blue
              $ plot_lines_style  . line_width .~ 5.0
              $ def

    layout = layout_title .~ t
           $ layout_plots .~ [toPlot actuals]
           $ layout_y_axis . laxis_generate .~ scaledAxis def (-yHt, yHt)
           $ def
