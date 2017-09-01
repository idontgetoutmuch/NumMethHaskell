{ nixpkgs ? import <nixpkgs> {}, compiler ? "default" }:

let

  inherit (nixpkgs) pkgs;

  f = { mkDerivation, base, datasets, hmatrix, mtl, random-fu
      , random-source, stdenv, typelits-witnesses, containers
      , ghc-prim, vector, cassava, bytestring
      }:
      mkDerivation {
        pname = "variational";
        version = "0.1.0.0";
        src = ./.;
        isLibrary = false;
        isExecutable = true;
        executableHaskellDepends = [
          base datasets hmatrix mtl random-fu random-source
          typelits-witnesses containers ghc-prim vector cassava bytestring
        ];
        license = stdenv.lib.licenses.bsd3;
      };

  haskellPackages = if compiler == "default"
                       then pkgs.haskellPackages
                       else pkgs.haskell.packages.${compiler};

  drv = haskellPackages.callPackage f {};

in

  if pkgs.lib.inNixShell then drv.env else drv
