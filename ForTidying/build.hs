{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE ScopedTypeVariables #-}

import Development.Shake
import Development.Shake.FilePath

main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
    want [ "diagrams" </> "symplectic" <.> "png"
         , "RunAccGPU" <.> "ll"
         , "TimeAccGPU" <.> "txt"
         , "TimeJulGPU" <.> "txt"
         , "RunJulGPU" <.> "ll"
         ]

    let compile file = do
          need [file]
          cmd "ghc --make -O2" file

    "diagrams" </> "symplectic" <.> "png" %> \_out -> do
      need ["SymplecticMain"]
      cmd "./SymplecticMain"

    "SymplecticMain" %> \out -> do
      compile (out -<.> "hs")

    "RunAccGPU" <.> "ll" %> \out -> do
      need ["RunAccGPU"]
      (Exit _code, Stdout (_stdout :: String), Stderr (stderr :: String)) <-
        cmd "./RunAccGPU +ACC -ddump-cc -dverbose -ACC"
      writeFileChanged out stderr

    "TimeAccGPU" <.> "txt" %> \out -> do
      need ["RunAccGPU"]
      (Exit _code, Stdout (_stdout :: String), Stderr (stderr :: String)) <-
        cmd "./RunAccGPU +RTS -s"
      writeFileChanged out stderr

    "RunAccGPU" %> \out -> do
      compile (out -<.> "hs")

    "TimeJulGPU" <.> "txt" %> \out -> do
      need ["JuliaCPU.jl"]
      (Exit _code, Stdout (_stdout :: String), Stderr (stderr :: String)) <-
        cmd "time /Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia JuliaCPU.jl"
      writeFileChanged out stderr

    "RunJulGPU" <.> "ll" %> \out -> do
      need ["JuliaLLVM.jl"]
      (Exit _code, Stdout (stdout :: String), Stderr (_stderr :: String)) <-
        cmd "/Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia JuliaLLVM.jl"
      writeFileChanged out stdout
