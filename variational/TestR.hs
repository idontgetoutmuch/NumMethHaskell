{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# OPTIONS_GHC -Wall #-}
{-# OPTIONS_GHC -fno-warn-unused-top-binds #-}
{-# LANGUAGE DataKinds,
 FlexibleContexts,
 TemplateHaskell,
 TypeOperators,
 TypeSynonymInstances,
 FlexibleInstances,
 QuasiQuotes,
 UndecidableInstances,
 ExplicitForAll,
 ScopedTypeVariables,
 OverloadedStrings,
 GADTs,
 TypeApplications
#-}

module Main (main) where

import Prelude as P

import Control.Lens

import Data.Foldable
import Frames hiding ( (:&) )

import Data.Vinyl (Rec(..))

import qualified Language.R as R
import Language.R (R)
import Language.R.QQ
import Data.Int

import Plots as P
import qualified Diagrams.Prelude as D
import Diagrams.Backend.Rasterific

import Data.Time.Clock(UTCTime, NominalDiffTime, diffUTCTime)
import Data.Time.Clock.POSIX(posixSecondsToUTCTime)

epochToUTC :: Integral a => a -> UTCTime
epochToUTC = posixSecondsToUTCTime . fromIntegral

tableTypes "OldFaithFul" "/Users/dom/Dropbox/Tidy/NumMethHaskell/variational/urk.csv"

declareColumn "eruptionTimeNew" ''UTCTime
declareColumn "eruptionTimeOld" ''UTCTime
declareColumn "eruptionTimeDiff" ''NominalDiffTime
declareColumn "eruptionTimeRat" ''Rational
declareColumn "eruptionTimeDoub" ''Double

loadOldFaithful :: IO (Frame OldFaithFul)
loadOldFaithful = inCoreAoS (readTable "/Users/dom/Dropbox/Tidy/NumMethHaskell/variational/urk.csv")

kSaxis :: [(Double, Double)] -> P.Axis B D.V2 Double
kSaxis xs = P.r2Axis &~ do
  P.scatterPlot' xs

getEruptionTimeNew :: OldFaithFul -> Record '[EruptionTimeNew]
getEruptionTimeNew x = pure (Col ((epochToUTC (x ^. eruptionTimeEpoch)))) :& Nil

getEruptionTimeOld :: OldFaithFul -> Record '[EruptionTimeOld]
getEruptionTimeOld x = pure (Col ((epochToUTC (x ^. eruptionTimeEpoch)))) :& Nil

getEruptionTimeDiff :: Record '[EruptionTimeOld, EruptionTimeNew] -> Record '[EruptionTimeDoub]
getEruptionTimeDiff x = pure (Col (fromRational $ toRational $ diffUTCTime (x ^. eruptionTimeNew) (x ^. eruptionTimeOld))) :& Nil

dropFrame :: Int -> Frame r -> Frame r
dropFrame n s = Frame ((frameLength s) - 1) (\i -> frameRow s (i + n))

getDuration :: OldFaithFul -> Record '["duration" :-> Text]
getDuration = (rcast @'[Duration])

eRs :: R s [Double]
eRs = R.dynSEXP <$> [r| faithful$eruptions |]

wRs :: R s [Int32]
wRs = R.dynSEXP <$> [r| faithful$waiting |]

eAltRs :: R s [Double]
eAltRs = R.dynSEXP <$> [r| geyser$duration |]

wAltRs :: R s [Int32]
wAltRs = R.dynSEXP <$> [r| geyser$waiting |]

main :: IO ()
main = do
  x <- loadOldFaithful
  let ds = fmap getDuration x
  let tn = dropFrame 1 $ fmap (\b -> (getEruptionTimeNew b)) x
  let tp = fmap (\b -> (getEruptionTimeOld b)) x
  let tpns = zipFrames tp tn
  let tds = zipFrames ds (fmap getEruptionTimeDiff tpns)
  let b = filterFrame (\u -> (u ^. eruptionTimeDoub) <= 9000) tds
  print $ frameLength b
  mapM_ print $ toList b
  R.withEmbeddedR R.defaultConfig $ do
  es <- R.runRegion $ do _ <- [r| library(MASS) |]
                         eAltRs
  ws <- R.runRegion $ do wAltRs
  renderRasterific "diagrams/TestR.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ P.zip es (P.map fromIntegral ws))
  es' <- R.runRegion $ do eRs
  ws' <- R.runRegion $ do wRs
  renderRasterific "diagrams/Test.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ P.zip es' (P.map fromIntegral ws'))
  putStrLn "Finished"

