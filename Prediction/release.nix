with import <nixpkgs> {};

let
  myHaskellPackages =
  haskellPackages.override {
  overrides = self: super: {
    hmatrix-sundials = self.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
      url = "https://github.com/haskell-numerics/hmatrix-sundials.git";
      rev = "9b6ec2b5fc509f74c5e61657dfc638a2c7ebced0";
    }) { sundials_arkode = sundials; sundials_cvode = sundials; };
  };
};

in

myHaskellPackages.callPackage ./default.nix { }
