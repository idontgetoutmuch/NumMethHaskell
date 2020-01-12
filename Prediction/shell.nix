{ pkgs ? import <nixpkgs> { overlays = [ ]; }, doBenchmark ? false }:

let

  hmatrix-sundials = pkgs.haskellPackages.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
    url = "https://github.com/haskell-numerics/hmatrix-sundials.git";
           rev = "9b6ec2b5fc509f74c5e61657dfc638a2c7ebced0";
  }) { sundials_arkode = pkgs.sundials; sundials_cvode = pkgs.sundials; };

  Naperian = pkgs.haskellPackages.callCabal2nix "Naperian" (builtins.fetchGit {
    url = "https://github.com/idontgetoutmuch/Naperian.git";
    rev = "54d873ffe99de865ca34e6bb3b92736e29e01619";
  }) { };

  haskellDeps = ps: with ps; [
    hmatrix
    hmatrix-sundials
    Naperian
    numbers
  ];

in

  pkgs.stdenv.mkDerivation {
  name = "env";
  buildInputs = [
    (pkgs.haskellPackages.ghcWithPackages haskellDeps)
  ];
}
