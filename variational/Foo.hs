module Foo where

import Frames.CSV
import Language.Haskell.TH

tableTypesText :: String -> FilePath -> DecsQ
tableTypesText n fp = tableTypesText' (rowGen fp) { rowTypeName = n }
