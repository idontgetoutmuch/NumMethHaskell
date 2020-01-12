{ stdenv
, cmake
, fetchurl
, openmpi
, llvmPackages
, python
, liblapack
, gfortran
, lapackSupport ? true }:

let liblapackShared = liblapack.override {
  shared = true;
};
    openmp = llvmPackages.openmp;

in stdenv.mkDerivation rec {
  pname = "sundials";
  version = "5.0.0";

  buildInputs = [ python openmpi openmp ] ++ stdenv.lib.optionals (lapackSupport) [ gfortran ];
  nativeBuildInputs = [ cmake openmpi openmp ];

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
    (fetchurl {
      # https://github.com/LLNL/sundials/pull/20
      url = "https://github.com/LLNL/sundials/pull/20/commits/2d951bbe1ff7842fcd0dafa28c61b0aa94015f66.patch";
      sha256 = "0lcr6m4lk14yqrxah4rdscpczny5l7m1zpfsjh8bgspadfsgk512";
    })
  ];

  cmakeFlags = [
    "-DEXAMPLES_INSTALL_PATH=${placeholder "out"}/share/examples"
  ] ++ stdenv.lib.optionals (lapackSupport) [

    "-DSUNDIALS_INDEX_SIZE=64"
    "-DMPI_ENABLE=ON"
    "-DMPI_C_COMPILER=${openmpi}/bin/mpicc"
    "-DMPI_CXX_COMPILER=${openmpi}/bin/mpicxx"
    "-DMPIEXEC_EXECUTABLE=${openmpi}/bin/mpirun"

    "-DOPENMP_ENABLE=ON"

    "-DBUILD_SHARED_LIBS=OFF"
    "-DBUILD_STATIC_LIBS=ON"

    "-DBUILD_CVODE=ON"
    "-DBUILD_CVODES=OFF"
    "-DBUILD_IDA=OFF"
    "-DBUILD_IDAS=OFF"
    "-DBUILD_ARKODE=ON"
    "-DBUILD_KINSOL=OFF"
    "-DBUILD_TESTING=ON"
    "-DSUNDIALS_DEVTESTS=ON"
    "-DEXAMPLES_ENABLE_CXX=ON"

    "-DLAPACK_ENABLE=OFF"
    "-DLAPACK_LIBRARIES=${liblapackShared}/lib/liblapack${stdenv.hostPlatform.extensions.sharedLibrary};${liblapackShared}/lib/libblas${stdenv.hostPlatform.extensions.sharedLibrary}"
  ];

  doCheck = false;
  preCheck = ''
    export OMP_NUM_THREADS=8
  '';
  checkPhase = "make test";

  meta = with stdenv.lib; {
    description = "Suite of nonlinear differential/algebraic equation solvers";
    homepage    = https://computation.llnl.gov/projects/sundials;
    platforms   = platforms.all;
    maintainers = with maintainers; [ flokli idontgetoutmuch ];
    license     = licenses.bsd3;
  };
}
