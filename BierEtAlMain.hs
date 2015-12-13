{-# OPTIONS_GHC -Wall                      #-}
{-# OPTIONS_GHC -fno-warn-name-shadowing   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults    #-}
{-# OPTIONS_GHC -fno-warn-unused-do-bind   #-}
{-# OPTIONS_GHC -fno-warn-missing-methods  #-}
{-# OPTIONS_GHC -fno-warn-orphans          #-}

module Main (
    main
    ) where

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Cairo.CmdLine

import BierEtAl
import BierEtAlChart

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )


-- x1 = Bier {vin =0.36, k1 =0.02, kp =6, km =15}
-- x2 = Bier {vin =0.36, k1 =0.02, kp =5, km =5}
-- x3 = Bier {vin =0.36, k1 =0.02, kp =6, km =10}
-- x4 = Bier {vin =0.2, k1 =0.02, kp =5, km =13}

-- x1 = Bier {vin =0.36, k1 =0.02, kp =4, km =15}
-- x2 = Bier {vin =0.36, k1 =0.02, kp =6, km = 7}
-- x3 = Bier {vin =0.2,  k1 =0.02, kp =5, km =13}
-- x4 = Bier {vin =0.1,  k1 =0.02, kp =6, km =13}

x1 = Bier {vin =0.36, k1 =0.01, kp =6, km =13}
x2 = Bier {vin =0.30, k1 =0.02, kp =6, km =18}
x3 = Bier {vin =0.50, k1 =0.02, kp =6, km =12}
x4 = Bier {vin =0.36, k1 =0.01, kp =7, km =13}

w1 = Bier {vin =1.10, k1 =0.02, kp =6, km =13}
w2 = Bier {vin =1.00, k1 =0.02, kp =6, km =13}
w3 = Bier {vin =0.95, k1 =0.02, kp =6, km =13}
w4 = Bier {vin =0.90, k1 =0.02, kp =6, km =13}

-- w1 = Bier {vin =0.36, k1 =0.02, kp =6, km=17}
-- w2 = Bier {vin =0.36, k1 =0.02, kp =6, km=18}
-- w3 = Bier {vin =0.36, k1 =0.02, kp =6, km=19}
-- w4 = Bier {vin =0.36, k1 =0.02, kp =6, km=20}

-- w1 = Bier {vin =0.36, k1 =0.02, kp =6, km=17}
-- w2 = Bier {vin =0.36, k1 =0.02, kp =6, km=18}
-- w3 = Bier {vin =0.36, k1 =0.02, kp =6, km=19}
-- w4 = Bier {vin =0.36, k1 =0.02, kp =6, km=20}

main :: IO ()
main = do
  let alt suffix xs = displayHeader ("diagrams/BierEtAl" ++ suffix ++ ".png")
                     (((diag ("BierEtAl" ++ suffix ++ "a")
                       (zip us' ((tol (xs!!0))!!0))
                       (zip us' ((tol (xs!!0))!!1)))
                      ===
                      (diag ("BierEtAl" ++ suffix ++ "b")
                       (zip us' ((tol (xs!!1))!!0))
                       (zip us' ((tol (xs!!1))!!1))))
                      |||
                      ((diag ("BierEtAl" ++ suffix ++ "c")
                       (zip us' ((tol (xs!!2))!!0))
                       (zip us' ((tol (xs!!2))!!1)))
                      ===
                      (diag ("BierEtAl" ++ suffix ++ "d")
                       (zip us' ((tol (xs!!3))!!0))
                       (zip us' ((tol (xs!!3))!!1)))))
  alt "11" [w1,w2,w3,w4]
  putStrLn "hello"
