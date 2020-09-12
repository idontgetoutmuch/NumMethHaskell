let overlay1 = self: super:
{
  sundials1 = self.callPackage ./SparseSundials { };
};

myHaskellPackageOverlay = self: super: {
  myHaskellPackages = super.haskellPackages.override {
    overrides = hself: hsuper: rec {

   tasty-golden =
        let newTastyGoldenSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/tasty-golden-2.3.3/tasty-golden-2.3.3.tar.gz";
          sha256 = "0wgcs4pqr30bp801cyrg6g551i7q0vjjmd9gmk5jy44fgdhb7kkl";
          };
            tg = hself.callCabal2nix "tasty-golden" newTastyGoldenSrc {};
          in
            super.haskell.lib.dontCheck tg;

    hmatrix-sundials1 = super.haskell.lib.dontCheck (
        hself.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
    url = "file:///Users/dom/hmatrix-sundials";
          rev = "70db2c221ce2bb7f5a26a6689178ea3e1121b86b";
}) { sundials_arkode          = self.sundials1;
     sundials_cvode           = self.sundials1;
     klu                      = self.suitesparse;
     suitesparseconfig        = self.suitesparse;
     sundials_sunlinsolklu    = self.sundials1;
     sundials_sunmatrixsparse = self.sundials1;
      });
    };
  };
};

in

{ nixpkgs ? import <nixpkgs> { overlays = [ overlay1 myHaskellPackageOverlay ]; } }:

let

inherit (nixpkgs) pkgs;
inherit (pkgs) myHaskellPackages;

haskellDeps = ps: with ps; [
  base
  cassava
  hmatrix-sundials1
  (nixpkgs.haskell.lib.dontCheck inline-r)
];

ghc = myHaskellPackages.ghcWithPackages haskellDeps;

nixPackages = [
  ghc
  myHaskellPackages.cabal-install
  myHaskellPackages.stack
  myHaskellPackages.lhs2tex
  pkgs.R
  pkgs.rPackages.ggplot2
  ];

in

pkgs.mkShell {
  buildInputs = [
    nixPackages
  ];
}
