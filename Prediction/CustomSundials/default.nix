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
  version = "5.0.0-dev.2";

  buildInputs = [ python openmpi openmp ] ++ stdenv.lib.optionals (lapackSupport) [ gfortran ];
  nativeBuildInputs = [ cmake openmpi openmp ];

  src = fetchurl {
    url = "https://computing.llnl.gov/projects/${pname}/download/${pname}-${version}.tar.gz";
    sha256 = "13hz8n8091kb0jdh6x18jsdv6fsv6akrsiv5ya1s991jc19vdq7y";
  };

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
