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

      BlogLiterately1 =
        let newBlogLiteratelySrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/BlogLiterately-0.8.7/BlogLiterately-0.8.7.tar.gz";
          sha256 = "10x5c3qn36fw50s0pgimxa4yhpbnx38bnqn5071wnyfh2z61clyn";
          };
            bl = hself.callCabal2nix "BlogLiterately" newBlogLiteratelySrc {};
          in
            super.haskell.lib.doJailbreak (super.haskell.lib.dontCheck bl);

      haxr =
        let newHaxrSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/haxr-3000.11.4.1/haxr-3000.11.4.1.tar.gz";
          sha256 = "1mm83k75zdnbx440zczk60ii4nkzll6dyhl3fzj1c4idb37r801r";
          };
            hr = hself.callCabal2nix "haxr" newHaxrSrc {};
          in
            super.haskell.lib.dontCheck hr;

      http-streams =
        let newHttpStreamSrc = builtins.fetchTarball { url = "https://hackage.haskell.org/package/http-streams-0.8.7.2/http-streams-0.8.7.2.tar.gz";
          sha256 = "1kz1rs89ii6mzb63h55fj3d7k7zwxi5b1ks04kak2gs9s50xykqh";
          };
            hs = hself.callCabal2nix "http-streams" newHttpStreamSrc {};
          in
            super.haskell.lib.dontCheck hs;

    hmatrix-sundials1 = super.haskell.lib.dontCheck (
        hself.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
    url = "file:///Users/dom/hmatrix-sundials";
           rev = "8333309c82b737247d3c14e66ae245febb1a6ab0";
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
  dimensional
  integration
  monad-loops
  hmatrix-sundials1
  numbers
  # Naperian
];

ghc = myHaskellPackages.ghcWithPackages haskellDeps;

my-python-packages = python-packages: with python-packages; [
  numpy
  matplotlib
  ];

python-with-my-packages = pkgs.python3.withPackages my-python-packages;

nixPackages = [
  ghc
  myHaskellPackages.cabal-install
  myHaskellPackages.stack
  myHaskellPackages.lhs2tex
  myHaskellPackages.BlogLiterately1
  ];

in

pkgs.mkShell {
  buildInputs = [
    nixPackages
    python-with-my-packages
    # pkgs.sage
  ];
  MYVARIABLE = "hi";
}
