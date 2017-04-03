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

import Tribbles

import qualified Prelude as P
import Numeric.Units.Dimensional.Prelude hiding (Unit)

import Text.Printf


denv :: DEnv Double
denv = unsafePerformIO $ defaultEnv vectorAlignmentFns 600 500

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

chart :: String ->
         [[(Double, Double)]] ->
         Renderable ()
chart t obss = toRenderable layout
  where

    actual x l c = plot_lines_values .~ [x]
                   $ plot_lines_style  . line_color .~ opaque c
                   $ plot_lines_title .~ l
                   $ plot_lines_style  . line_width .~ 1.0
                   $ def

    ls = map (\n -> "Path " ++ show n) [1..]
    cs = cycle [blue, green, red, brown, crimson]

    actuals' :: [PlotLines Double Double]
    actuals' = zipWith3 actual obss ls cs

    layout = layout_title .~ t
           $ layout_plots .~ (map toPlot actuals') -- [toPlot actuals]
           $ layout_y_axis . laxis_title .~ "Value"
           $ layout_y_axis . laxis_override .~ axisGridHide
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

mus :: [Double]
mus = map (/~ day ) $ iterate (+ deltaMu) (0.0 *~ day)

nus :: [Double]
nus = map (/~ day ) $ iterate (+ deltaNu) (0.0 *~ day)

-- pss :: [Double] -> [[(Double, Double)]]
-- pss us = map (\(p, _, _, _) -> zip us (concat $ toLists p)) ys

selNth :: Int -> [a] -> [a]
selNth n xs = map snd $ filter (\z -> fst z `mod` n == 0) (zip [0..] xs)

-- pss' :: [[(Double, Double)]]
-- pss' = take 15 $ selNth 10 (pss mus)

-- mss :: [Double] -> [[(Double, Double)]]
-- mss us = map (\(_, m, _, _, _) -> zip us (concat $ toLists m)) ys

ps :: [Double]
ps = map (\(p, _, _) -> p /~ (mole / kilo gram)) ys

ms :: [Double]
ms = map (\(_, m, _) -> m /~ (mole / kilo gram)) ys

es :: [Double]
es = map (\(_, _, e) -> e /~ (mole / (metre * metre * metre))) ys

diagrm :: String -> [[(Double, Double)]] -> Diagram Cairo
diagrm t xss = fst $ runBackendR denv (chart t xss)

main :: IO ()
main = do
  mapM_ (\(x, y, z) -> printf "%4.3e, %4.3e, %4.3e\n"
                             (x /~ (mole / kilo gram))
                             (y /~ (mole / kilo gram))
                             (z /~ (mole / (metre * metre * metre))))
        ys
  displayHeader "diagrams/EPO.png"
                (diagrm "EPO" [zip (map (/~ day) ts) es])
  displayHeader "diagrams/Precursors.png"
                (diagrm "Precursors" [zip (map (/~ day) ts) ps])
  displayHeader "diagrams/Matures.png"
                (diagrm "Matures" [zip (map (/~ day) ts) (map (P.* 1e-11) ms)])
  -- displayHeader "diagrams/PrecursorPop.png"
  --               (diagrm "Precursor Population" pss')
  -- displayHeader "diagrams/InitialMature.png"
  --               (diagrm "Initial Mature" (mss nus))

  -- putStrLn "Hello"
