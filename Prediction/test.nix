let overlay1 = self: super:
{
  sundials = self.callPackage ./CustomSundials { };
};

in

{ nixpkgs ? import <nixpkgs> { overlays = [ overlay1 ]; } }:

  nixpkgs.stdenv.mkDerivation {
  name = "env";
  buildInputs = [
    nixpkgs.openmpi
    nixpkgs.openssh
    nixpkgs.sundials
  ];
}
