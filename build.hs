import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util

main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_build"} $ do
    want ["diagrams" </> "symplectic" <.> "png"
         ,"RunAccGPU" <.> "ll"]

    let compile file = do
          need [file]
          cmd "ghc --make -O2" file

    "diagrams" </> "symplectic" <.> "png" %> \out -> do
      need ["SymplecticMain"]
      cmd "./SymplecticMain"

    "SymplecticMain" %> \out -> do
      compile (out -<.> "hs")

    "RunAccGPU" <.> "ll" %> \out -> do
      need ["RunAccGPU"]
      (Exit code, Stdout (stdout :: String), Stderr (stderr :: String)) <-
        cmd "./RunAccGPU +ACC -ddump-cc -dverbose -ACC"
      writeFileChanged out stderr

    "RunAccGPU" %> \out -> do
      compile (out -<.> "hs")

