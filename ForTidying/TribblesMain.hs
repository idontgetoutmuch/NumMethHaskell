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
           $ layout_x_axis . laxis_title .~ "Time"
           $ layout_x_axis . laxis_override .~ axisGridHide
           $ def

mus :: [Double]
mus = map (/~ day ) $ iterate (+ deltaMu) (0.0 *~ day)

nus :: [Double]
nus = map (/~ day ) $ iterate (+ deltaNu) (0.0 *~ day)

selNth :: Int -> [a] -> [a]
selNth n xs = map snd $ filter (\z -> fst z `mod` n == 0) (zip [0..] xs)

ps :: [Double]
ps = map (\(p, _, _) -> p /~ (mole / kilo gram)) ys

ms :: [Double]
ms = map (\(_, m, _) -> m /~ (mole / kilo gram)) ys

es :: [Double]
es = map (\(_, _, e) -> e /~ muPerMl) ys

diagrm :: String -> String -> [[(Double, Double)]] -> Diagram Cairo
diagrm t l xss = fst $ runBackendR denv (chart t l xss)

diagrmM :: String -> String -> [[(Double, Double)]] -> IO (Diagram Cairo)
diagrmM t l xss = do
  denv <- defaultEnv vectorAlignmentFns 600 500
  return $ fst $ runBackendR denv (chart t l xss)

main :: IO ()
main = do
  mapM_ (\(x, y, z) -> printf "%4.3e, %4.3e, %4.3e\n"
                             (x /~ (mole / kilo gram))
                             (y /~ (mole / kilo gram))
                             (z /~ muPerMl))
        ys
  displayHeader "diagrams/EPO.png"
                (diagrm "EPO" "mUI / mL" [zip (map (/~ day) ts) es])
  displayHeader "diagrams/Precursors.png"
                (diagrm "Precursors" "Cells / Kilogram " [zip (map (/~ day) ts) ps])
  displayHeader "diagrams/Matures.png"
                (diagrm "Mature Erythrocytes" "Cells / Kilogram x 1e11" [zip (map (/~ day) ts) (map (P.* 1e-11) ms)])

