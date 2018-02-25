{ nixpkgs ? import <nixpkgs> { overlays = [(self: super: {gsl = super.gsl.overrideAttrs (o: {CFLAGS = "-DDEBUG -O3";});})]; }
, compiler ? "ghc822"
, doBenchmark ? false }:

let

  inherit (nixpkgs) pkgs;

f = { mkDerivation, array, base, bytestring, cassava, containers
    , datasets, diagrams-lib, diagrams-rasterific, foldl, Frames, ghc-prim, gsl
    , hmatrix, hmatrix-gsl
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
    gsl
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
