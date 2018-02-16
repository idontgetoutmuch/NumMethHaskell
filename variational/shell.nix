{ nixpkgs ? import <nixpkgs> {}, compiler ? "ghc822", doBenchmark ? false }:

let

  inherit (nixpkgs) pkgs;

f = { mkDerivation, array, base, bytestring, cassava, containers
    , datasets, diagrams-lib, diagrams-rasterific, ghc-prim, hmatrix
    , inline-r, mtl, plots, random-fu, R, random-source
    , stdenv, typelits-witnesses, vector }:
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
    ghc-prim
    hmatrix
    inline-r
    mtl
    plots
    random-fu
    random-source
    typelits-witnesses
    vector
  ];
  executableSystemDepends = [
    R
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
