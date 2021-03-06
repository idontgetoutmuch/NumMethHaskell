# nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs-channels/archive/nixos-20.03.tar.gz minimal.nix

let

# nixpkgsRev = "bc260badaebf67442befe20fb443034d3a91f2b3"; # 20.09-beta
# nixpkgsSha256 = "1iysc4xyk88ngkfb403xfq5bs3zy29zfs83pn99kchxd45nbpb5q";

# mixpkgs = fetchTarball {
#   url = "https://github.com/nixos/nixpkgs/archive/${nixpkgsRev}.tar.gz";
#   sha256 = nixpkgsSha256;
# };

# mixpkgs = fetchTarball {
#   url = "https://github.com/NixOS/nixpkgs-channels/archive/nixos-20.03.tar.gz";
#   sha256 = "1qbs7p0mmcmpg70ibd437hl57byqx5q0pc61p1dckrkazj7kq0pc";
# };

overlay1 = self: super:
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

      Frames =
        let newFramesSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/Frames-0.7.0/Frames-0.7.0.tar.gz";
          sha256 = "1dxdwcz21rk83rp8b9l1jahdicmhywkpaa1n70fp7ihc07jghxmh";
          };
            f = hself.callCabal2nix "Frames" newFramesSrc {};
          in
            super.haskell.lib.dontCheck f;

      vinyl =
        let newVinylSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/vinyl-0.13.0/vinyl-0.13.0.tar.gz";
          sha256 = "05kj1ld70yrxfff3zlc0vr8b6r8pawikcmiyvv1cc8c40jypaji4";
          };
            v = hself.callCabal2nix "vinyl" newVinylSrc {};
          in
            super.haskell.lib.dontCheck v;

      # inline-r =
      #   let newInlineRSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/inline-r-0.10.3/inline-r-0.10.3.tar.gz";
      #     sha256 = "10v4d6ka27pxi54zcq1c3ic140v9dx6lb62l5dmxjmlb6522n3g7";
      #     };
      #       ir = hself.callCabal2nix "inline-r" newInlineRSrc {};
      #     in
      #       super.haskell.lib.dontCheck ir;

      kalman1 = super.haskell.lib.dontCheck (
        hself.callCabal2nix "kalman" (builtins.fetchGit {
    url = "file:///Users/dom/Kalman";
          rev = "2314229a14ae25143f675a9eb08d8767a9c1ac56";
}) { });

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
  Frames
  hmatrix-sundials1
  kalman1
  lens
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
