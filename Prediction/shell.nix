{ pkgs ? import <nixpkgs> {}, doBenchmark ? false }:

let

  hmatrix-sundials = pkgs.haskellPackages.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
    url = "https://github.com/haskell-numerics/hmatrix-sundials.git";
           rev = "9b6ec2b5fc509f74c5e61657dfc638a2c7ebced0";
  }) { sundials_arkode = pkgs.sundials; sundials_cvode = pkgs.sundials; };

  haskellDeps = ps: with ps; [
    Chart
    Chart-diagrams
    diagrams-cairo
    diagrams-lib
    diagrams-rasterific
    hmatrix
    hmatrix-sundials
    Naperian
  ];

in

  pkgs.stdenv.mkDerivation {
  name = "env";
  buildInputs = [
    (pkgs.haskellPackages.ghcWithPackages haskellDeps)
  ];
}
