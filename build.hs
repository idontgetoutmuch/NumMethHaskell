import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
    want ["diagrams" </> "symplectic" <.> "png"]

    let compile file = do
          need [file]
          cmd "ghc --make -O2" file

    "diagrams" </> "symplectic" <.> "png" %> \out -> do
      need ["SymplecticMain"]
      cmd "./SymplecticMain"

    "SymplecticMain" %> \out -> do
      compile (out -<.> "hs")

    -- phony "clean" $ do
    --     putNormal "Cleaning files in _build"
    --     removeFilesAfter "_build" ["//*"]

    -- "_build/run" <.> exe %> \out -> do
    --     let hs = ["SymplecticMain.hs"]
    --     let os = ["_build" </> h -<.> "o" | h <- hs]
    --     need os
    --     cmd "ghc -O2" [out] os

    -- "_build/run" <.> exe %> \out -> do
    --     cs <- getDirectoryFiles "" ["//*.c"]
    --     let os = ["_build" </> c -<.> "o" | c <- cs]
    --     need os
    --     cmd "gcc -o" [out] os

    -- "_build//*.o" %> \out -> do
    --     let c = dropDirectory1 $ out -<.> "c"
    --     let m = out -<.> "m"
    --     () <- cmd "gcc -c" [c] "-o" [out] "-MMD -MF" [m]
    --     needMakefileDependencies m

