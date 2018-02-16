{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE ScopedTypeVariables #-}

import qualified Language.R as R
import Language.R (R)
import Language.R.QQ
import Data.Int

import Plots as P
import Diagrams.Prelude hiding ( normal, sample, (<>), inv, trace )
import Diagrams.Backend.Rasterific

kSaxis :: [(Double, Double)] -> P.Axis B V2 Double
kSaxis xs = P.r2Axis &~ do
  P.scatterPlot' xs

eRs :: R s [Double]
eRs = R.dynSEXP <$> [r| faithful$eruptions |]

wRs :: R s [Int32]
wRs = R.dynSEXP <$> [r| faithful$waiting |]

eAltRs :: R s [Double]
eAltRs = R.dynSEXP <$> [r| geyser$duration |]

wAltRs :: R s [Int32]
wAltRs = R.dynSEXP <$> [r| geyser$waiting |]

main = R.withEmbeddedR R.defaultConfig $ do
  es <- R.runRegion $ do _ <- [r| library(MASS) |]
                         eAltRs
  ws <- R.runRegion $ do wAltRs
  renderRasterific "diagrams/TestR.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip es (map fromIntegral ws))
  es <- R.runRegion $ do eRs
  ws <- R.runRegion $ do wRs
  renderRasterific "diagrams/TestS.png"
                   (dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ zip es (map fromIntegral ws))
  putStrLn "Finished"

