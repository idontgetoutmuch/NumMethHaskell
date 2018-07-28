{-# LANGUAGE TemplateHaskell   #-}
{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE TypeOperators     #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RankNTypes        #-}

import           Frames hiding ((:&))
import qualified Control.Foldl as Foldl
import           Control.Lens
import           Pipes.Prelude (fold)
import           Pipes
import qualified Pipes.Prelude as P
import qualified Pipes.Core as PC
import qualified Data.Foldable as F

import           Foo

import           Data.Void

import qualified Data.Vinyl as V
import           Data.Vinyl (Rec(..))

import           Data.Text hiding (map)
import           Data.Attoparsec.Text
import           Control.Applicative



tableTypesText "ProductionData" "/Users/dom/Downloads/2017_prod_reports_from_excel.csv"

declareColumn "oil_vol_double" ''Double

getOilVol :: ProductionData -> Record '[OilVolDouble]
getOilVol x = pure (Col ((f $ (parseOnly parseDouble) (x ^. oilVol)))) :& Nil
  where
    f (Left _)  = 0.0 / 0.0
    f (Right y) = y

cleanData :: ProductionData -> Record '[ReportMonth, ReportYear, OilVolDouble]
cleanData x = pure (Col (x ^. reportMonth)) :& pure (Col (x ^. reportYear)) :& getOilVol x

readCleanData :: MonadSafe m => Producer (Record '[ReportMonth, ReportYear, OilVolDouble]) m ()
readCleanData = (readTable "/Users/dom/Downloads/2017_prod_reports_from_excel.csv") >->
                (P.map cleanData)

readOilVol :: MonadSafe m => Producer (Record '[OilVolDouble]) m ()
readOilVol = (readTable "/Users/dom/Downloads/2017_prod_reports_from_excel.csv") >->
             (P.map getOilVol)

oilVolLength :: Foldl.Fold (Record '[OilVolDouble]) Int
oilVolLength = Foldl.length

totalOilVol :: Foldl.Fold (Record '[OilVolDouble]) Double
totalOilVol = (Foldl.handles oilVolDouble) Foldl.sum

oilVolTotalAndLength :: Foldl.Fold (Record '[OilVolDouble]) (Double, Int)
oilVolTotalAndLength = (,) <$> totalOilVol <*> oilVolLength

oilVolFiltered :: Text -> Foldl.Fold (Record '[ReportMonth, ReportYear, OilVolDouble]) Double
oilVolFiltered m = (Foldl.handles (filtered (\x -> x ^. reportMonth == m)) .
                    Foldl.handles oilVolDouble) Foldl.sum

allMonths :: Foldl.Fold (Record '[ReportMonth, ReportYear, OilVolDouble]) [Double]
allMonths = sequenceA $ map oilVolFiltered $ map pack $ map show [1..12]

parseDouble :: Parser Double
parseDouble =
  do d <- double
     return d
  <|>
  do _ <- string ""
     return 0.0

main = do
  (t, l) <- runSafeT $
            Foldl.purely fold oilVolTotalAndLength readOilVol
  putStrLn $ show l ++ " records totalling " ++ show t
  f <- runSafeT $
       Foldl.purely fold allMonths readCleanData
  putStrLn $ show f

