{ mkDerivation, base, cassava, Chart, Chart-diagrams
, diagrams-cairo, diagrams-lib, diagrams-rasterific, hmatrix
, hmatrix-sundials, Naperian, stdenv
}:
mkDerivation {
  pname = "Prediction";
  version = "0.1.0.0";
  src = ./.;
  isLibrary = false;
  isExecutable = true;
  executableHaskellDepends = [
    base cassava Chart Chart-diagrams diagrams-cairo diagrams-lib
    diagrams-rasterific hmatrix hmatrix-sundials Naperian
  ];
  license = stdenv.lib.licenses.bsd3;
}
