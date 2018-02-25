% Data Sources
% Dominic Steinitz
% 18th February 2018

---
bibliography: DynSys.bib
---

Introduction
============

For the blog post still being written on variatonal methods, I
referred to the still excellent @bishop2006pattern who uses as his
example data, the data available in
[R](https://www.r-project.org) for the
geyser in Yellowstone National Park called "Old Faithful".

While explaining this to another statistician, they started to ask
about the dataset. Since I couldn't immediately answer their questions
(when was it collected? over how long? how was it collected?), I started
to look at the dataset more closely. The more I dug into where the
data came from, the less certain I became that the dataset actually
reflected how Old Faithful actually behaves.

Here's what I found. On the way, I use
[Haskell](https://www.haskell.org), [R embedded in
Haskell](https://tweag.github.io/HaskellR), [data frames in
Haskell](https://hackage.haskell.org/package/Frames) and
[nix](https://nixos.org/nix).


The Investigation
=================

First some imports and language extensions.

\begin{code}
{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-unused-top-binds #-}

{-# LANGUAGE QuasiQuotes          #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE DataKinds            #-}
{-# LANGUAGE FlexibleContexts     #-}
{-# LANGUAGE TemplateHaskell      #-}
{-# LANGUAGE TypeOperators        #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances    #-}
{-# LANGUAGE QuasiQuotes          #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE ExplicitForAll       #-}
{-# LANGUAGE ScopedTypeVariables  #-}
{-# LANGUAGE OverloadedStrings    #-}
{-# LANGUAGE GADTs                #-}
{-# LANGUAGE TypeApplications     #-}


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

import Data.Attoparsec.Text
import Control.Applicative
\end{code}

R provides many
[datasets](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/00Index.html)
presumably for ease of use. We can access the "Old Faithful" data set
in Haskell using some very simple inline R.

\begin{code}
eRs :: R s [Double]
eRs = R.dynSEXP <$> [r| faithful$eruptions |]

wRs :: R s [Int32]
wRs = R.dynSEXP <$> [r| faithful$waiting |]
\end{code}

And then plot the dataset using

\begin{code}
kSaxis :: [(Double, Double)] -> P.Axis B D.V2 Double
kSaxis xs = P.r2Axis &~ do
  P.scatterPlot' xs

plotOF :: IO ()
plotOF = do
  es' <- R.runRegion $ do eRs
  ws' <- R.runRegion $ do wRs
  renderRasterific "diagrams/Test.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ P.zip es' (P.map fromIntegral ws'))
\end{code}

![](diagrams/testS.png)

In the
[documentation](https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/faithful.html)
for this dataset, the source is given as @härdle2012smoothing. In here
the data are given as a table in Table 3 on page 201. It contains 270
observations not 272 as the R dataset does! Note that there seems to
be a correlation between inter-eruption time and duration of eruption
and eruptions seem to fall into one of two groups.

The documentation also says "See Also geyser in package MASS for the
Azzalini–Bowman version". So let's look at this and plot it.

\begin{code}
eAltRs :: R s [Double]
eAltRs = R.dynSEXP <$> [r| geyser$duration |]

wAltRs :: R s [Int32]
wAltRs = R.dynSEXP <$> [r| geyser$waiting |]

plotOF' :: IO ()
plotOF' = do
  es <- R.runRegion $ do _ <- [r| library(MASS) |]
                         eAltRs
  ws <- R.runRegion $ do wAltRs
  renderRasterific "diagrams/TestR.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ P.zip es (P.map fromIntegral ws))

\end{code}

A quick look at the first few entries suggests this is the data in
Table 1 of @azzalini but NB I didn't check each item. Let's plot this also

![](diagrams/testR.png)

It looks quite different to the *faithful* dataset. But to
answer my fellow statistician's questions: "This analysis deals with
data which were collected continuously from August 1st until August
15th, 1985".

It turns out there is a
[website](https://www.geysertimes.org/about.php) which collects data
for a large number of geysers and has [historical
datasets](https://www.geysertimes.org/archive) for them including old
faithful. The data is in a slightly awkward format and it is not clear
(to me at any rate) what some of the [fields
mean](https://github.com/glennon/GeyserTimes-Science/issues/32).

Let's examine the first 1000 observations in the data set using the
Haskell package [Frames](https://hackage.haskell.org/package/Frames).

First we ask *Frames* to automatically type the data and read it into a frame.

\begin{code}
tableTypes "OldFaithFul" "/Users/dom/Dropbox/Tidy/NumMethHaskell/variational/Old_Faithful_eruptions_1000.csv"

loadOldFaithful :: IO (Frame OldFaithFul)
loadOldFaithful = inCoreAoS (readTable "/Users/dom/Dropbox/Tidy/NumMethHaskell/variational/Old_Faithful_eruptions_1000.csv")
\end{code}

Declare some columns and helper functions to manipulate the times.

\begin{code}
declareColumn "eruptionTimeNew"  ''UTCTime
declareColumn "eruptionTimeOld"  ''UTCTime
declareColumn "eruptionTimeDiff" ''NominalDiffTime
declareColumn "eruptionTimeRat"  ''Rational
declareColumn "eruptionTimeDoub" ''Double
declareColumn "parsedDuration"   ''Double

epochToUTC :: Integral a => a -> UTCTime
epochToUTC = posixSecondsToUTCTime . fromIntegral

getEruptionTimeNew :: OldFaithFul -> Record '[EruptionTimeNew]
getEruptionTimeNew x = pure (Col ((epochToUTC (x ^. eruptionTimeEpoch)))) :& Nil

getEruptionTimeOld :: OldFaithFul -> Record '[EruptionTimeOld]
getEruptionTimeOld x = pure (Col ((epochToUTC (x ^. eruptionTimeEpoch)))) :& Nil

getEruptionTimeDiff :: Record '[EruptionTimeOld, EruptionTimeNew] -> Record '[EruptionTimeDoub]
getEruptionTimeDiff x = pure (Col ((/60.0) $ fromRational $ toRational $ diffUTCTime (x ^. eruptionTimeNew) (x ^. eruptionTimeOld))) :& Nil
\end{code}

The duration is expressed as a string mainly as e.g. "3.5m" and
"3.5min" but occasionally in other formats. Using NaNs to represent
data we can't parse is poor practice but *Frames* seems to object to
types such as *Maybe Double*.

\begin{code}
getDuration :: OldFaithFul -> Record '["duration" :-> Text]
getDuration = (rcast @'[Duration])

getParsedDuration :: Record '[Duration, EruptionTimeDoub] -> Record '[ParsedDuration]
getParsedDuration x = pure (Col ((f $ (parseOnly parseDuration) (x ^. duration)))) :& Nil
  where
    f (Left _)  = 0.0 / 0.0
    f (Right y) = y

parseDuration :: Parser Double
parseDuration =
  do d <- double
     _ <- char 'm'
     return d
  <|>
  do d <- double
     _ <- string "min"
     return d
\end{code}

To get the times between eruptions we need to zip the frame with its
tail so we define this helper function.

\begin{code}
dropFrame :: Int -> Frame r -> Frame r
dropFrame n s = Frame ((frameLength s) - 1) (\i -> frameRow s (i + n))
\end{code}

\begin{code}
main :: IO ()
main = do
  x <- loadOldFaithful
\end{code}

Get the current eruption times and the previous eruption times and put them in a frame.

\begin{code}
  let tn = dropFrame 1 $ fmap (\b -> (getEruptionTimeNew b)) x
  let tp = fmap (\b -> (getEruptionTimeOld b)) x
  let tpns = zipFrames tp tn
\end{code}

Get the durations of the eruptions and exclude any gaps when no
observations were taken; arbitrarily if the gap is greate than 1.5
hours.

\begin{code}
  let ds = fmap getDuration x
  let tds = zipFrames ds (fmap getEruptionTimeDiff tpns)
  let b = filterFrame (\u -> (u ^. eruptionTimeDoub) <= 150) tds
\end{code}

Parse some of the durations and put the durations and the
inter-eruption gap into a single frame.

\begin{code}
  let c = filterFrame (\u -> not $ isNaN $ u ^. parsedDuration) $
          fmap getParsedDuration b
  let d = zipFrames b c
  let g = fmap (\u -> (u ^. parsedDuration, u ^. eruptionTimeDoub)) d
\end{code}

Finally we can yet another data set of observations of old faithful.

\begin{code}
  renderRasterific "diagrams/TestH1000.png"
                   (D.dims2D 500.0 500.0)
                   (renderAxis $ kSaxis $ toList g)
  R.withEmbeddedR R.defaultConfig $ do
    plotOF
    plotOF'
  putStrLn "Finished"
\end{code}

![](diagrams/TestH1000.png)

To me this doesn't look like either of the other two datasets. Perhaps
the R *faithful* data set should be validated, possibly retired and
maybe replaced by something with firmer foundations?

Nix
===

I use [nix](https://nixos.org) on MACos for package management both
for Haskell and R - no more *install.packages* problems. Here's my *shell.nix* file also available here: 

``` {.nix}
{ nixpkgs ? import <nixpkgs> {}, compiler ? "ghc822", doBenchmark ? false }:

let

  inherit (nixpkgs) pkgs;

f = { mkDerivation, array, base, bytestring, cassava, containers
    , datasets, diagrams-lib, diagrams-rasterific, foldl, Frames, ghc-prim, hmatrix, hmatrix-gsl
    , inline-r, lens, mtl, pipes, plots, random-fu, R, random-source
    , stdenv, template-haskell, text, typelits-witnesses, vector, vinyl }:
mkDerivation {
  pname = "variational";
  version = "0.1.0.0";
  src = ./.;
  isLibrary = false;
  isExecutable = true;
  executableHaskellDepends = [
    array
    base
    bytestring
    cassava
    containers
    datasets
    diagrams-lib
    diagrams-rasterific
    foldl
    Frames
    ghc-prim
    hmatrix
    inline-r
    lens
    mtl
    pipes
    plots
    random-fu
    random-source
    template-haskell
    text
    typelits-witnesses
    vector
    vinyl
  ];
  executableSystemDepends = [
    R
    pkgs.rPackages.anytime
    pkgs.rPackages.ggplot2
    pkgs.rPackages.maptools
    pkgs.rPackages.reshape2
    pkgs.rPackages.rgeos
    pkgs.rPackages.rgdal
    pkgs.rPackages.rstan ];
  license = stdenv.lib.licenses.bsd3;
};

haskellPackages = if compiler == "default"
  then pkgs.haskellPackages
  else pkgs.haskell.packages.${compiler};

variant = if doBenchmark then pkgs.haskell.lib.doBenchmark else pkgs.lib.id;

drv = variant (haskellPackages.callPackage f {});

in

  if pkgs.lib.inNixShell then drv.env else drv
```

References
==========
