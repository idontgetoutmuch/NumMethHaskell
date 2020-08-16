{ stdenv
, cmake
, fetchurl
, llvmPackages
, python
, liblapack
, gfortran
, suitesparse
, lapackSupport ? true }:

let liblapackShared = liblapack.override {
  shared = true;
};

in stdenv.mkDerivation rec {
  pname = "sundials";
  version = "5.0.0";

  buildInputs = stdenv.lib.optionals (lapackSupport) [ gfortran ];
  nativeBuildInputs =  [ cmake ];

  src = fetchurl {
    url = "https://computing.llnl.gov/projects/${pname}/download/${pname}-${version}.tar.gz";
    sha256 = "1lvx5pddjxgyr8kqlira36kxckz7nxwc8xilzfyx0hf607n42l9l";
  };

  patches = [
    (fetchurl {
      # https://github.com/LLNL/sundials/pull/19
      url = "https://github.com/LLNL/sundials/commit/1350421eab6c5ab479de5eccf6af2dcad1eddf30.patch";
      sha256 = "0g67lixp9m85fqpb9rzz1hl1z8ibdg0ldwq5z6flj5zl8a7cw52l";
    })
  ];

  cmakeFlags = [
    "-DEXAMPLES_INSTALL_PATH=${placeholder "out"}/share/examples"
  ] ++ stdenv.lib.optionals (lapackSupport) [

    "-DSUNDIALS_INDEX_SIZE=64"

    "-DBUILD_SHARED_LIBS=OFF"
    "-DBUILD_STATIC_LIBS=ON"

    "-DBUILD_CVODE=ON"
    "-DBUILD_CVODES=OFF"
    "-DBUILD_IDA=OFF"
    "-DBUILD_IDAS=OFF"
    "-DBUILD_ARKODE=ON"
    "-DBUILD_KINSOL=OFF"
    "-DBUILD_TESTING=ON"
    "-DEXAMPLES_ENABLE_C=OFF"
    "-DEXAMPLES_ENABLE_CXX=OFF"
    "-DEXAMPLES_INSTALL=OFF"

    "-DKLU_ENABLE=ON"
    "-DKLU_INCLUDE_DIR=${suitesparse}/include"
    "-DKLU_LIBRARY_DIR=${suitesparse}/lib"

    "-DLAPACK_ENABLE=ON"
    "-DLAPACK_LIBRARIES=${liblapackShared}/lib/liblapack${stdenv.hostPlatform.extensions.sharedLibrary};${liblapackShared}/lib/libblas${stdenv.hostPlatform.extensions.sharedLibrary}"
  ];

  doCheck = false;

  meta = with stdenv.lib; {
    description = "Suite of nonlinear differential/algebraic equation solvers";
    homepage    = https://computation.llnl.gov/projects/sundials;
    platforms   = platforms.all;
    maintainers = with maintainers; [ flokli idontgetoutmuch ];
    license     = licenses.bsd3;
  };
}
