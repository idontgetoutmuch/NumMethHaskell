let
  pkgs = {
    ihaskell = builtins.fetchTarball {
      url = "https://github.com/gibiansky/IHaskell/tarball/bb2500c448c35ca79bddaac30b799d42947e8774";
      sha256 = "1n4yqxaf2xcnjfq0r1v7mzjhrizx7z5b2n6gj1kdk2yi37z672py";
    };
    # nixpkgs = builtins.fetchTarball {
    #   url = "https://github.com/NixOS/nixpkgs-channels/archive/nixos-19.09.tar.gz";
    #   sha256 = "16wdsazc7g09ibcxlqsa3kblzhbbpdpb6s29llliybw73cp37b9s";
    # };
    nixpkgs = builtins.fetchTarball {
      url = "https://github.com/NixOS/nixpkgs-channels/tarball/49dc8087a20e0d742d38be5f13333a03d171006a";
      sha256 = "1fdnqm4vyj50jb2ydcc0nldxwn6wm7qakxfhmpf72pz2y2ld55i6";
    };
  };

  rOverlay = rself: rsuper: {
    myR = rsuper.rWrapper.override {
      packages = with rsuper.rPackages; [ ggplot2 dplyr xts purrr cmaes cubature ];
    };
  };

  hmatrix-sundials = nixpkgs.haskell.lib.dontCheck (
    nixpkgs.haskellPackages.callCabal2nix "hmatrix-sundials" (builtins.fetchGit {
      url = "https://github.com/haskell-numerics/hmatrix-sundials.git";
      rev = "9b6ec2b5fc509f74c5e61657dfc638a2c7ebced0";
    }) { sundials_arkode = nixpkgs.sundials; sundials_cvode = nixpkgs.sundials; });

  haskellOverlay = self: super: {
    haskell = super.haskell // { packageOverrides =
	    hself: hsuper: {
        my-random-fu-multivariate = hself.callPackage /Users/dom/nix-config/pkgs/random-fu-multivariate { };
        my-hmatrix-sundials = hmatrix-sundials;

        Chart = super.haskell.lib.dontCheck (
          hself.callCabal2nixWithOptions "Chart" (builtins.fetchGit {
            url = "https://github.com/timbod7/haskell-chart";
            rev = "65606988864938f71ea79e3bc09872e7bab54f19";
          }) "--subpath chart" { });

        Chart-cairo = super.haskell.lib.dontCheck (
          hself.callCabal2nixWithOptions "Chart-cairo" (builtins.fetchGit {
            url = "https://github.com/timbod7/haskell-chart";
            rev = "65606988864938f71ea79e3bc09872e7bab54f19";
          }) "--subpath chart-cairo" { });

        Chart-diagrams = super.haskell.lib.dontCheck (
          hself.callCabal2nixWithOptions "Chart-diagrams" (builtins.fetchGit {
            url = "https://github.com/timbod7/haskell-chart";
            rev = "65606988864938f71ea79e3bc09872e7bab54f19";
          }) "--subpath chart-diagrams" { });

        Naperian = hself.callCabal2nix "Naperian" (builtins.fetchGit {
          url = "https://github.com/idontgetoutmuch/Naperian.git";
          rev = "54d873ffe99de865ca34e6bb3b92736e29e01619";
          }) { };
      };
    };
  };

  nixpkgs  = import pkgs.nixpkgs { overlays = [ rOverlay haskellOverlay ]; config.allowBroken = true; };

  r-libs-site = nixpkgs.runCommand "r-libs-site" {
    buildInputs = with nixpkgs; [ R rPackages.ggplot2 rPackages.dplyr rPackages.xts rPackages.purrr rPackages.cmaes rPackages.cubature ];
  } ''echo $R_LIBS_SITE > $out'';

  ihaskellEnv = (import "${pkgs.ihaskell}/release.nix" {
    compiler = "ghc864";
    nixpkgs  = nixpkgs;
  packages = self: [
    self.inline-r
    self.hmatrix
    self.my-hmatrix-sundials
    self.random-fu
    self.my-random-fu-multivariate
    self.cassava
    self.diagrams
    self.ihaskell-diagrams
    self.ihaskell-charts
    self.Chart
    self.Chart-diagrams
    self.Chart-cairo
    self.Naperian
    self.numbers
  ];
  }).passthru.ihaskellEnv;

  systemPackages = self: [ self.myR ];

  jupyterlab = nixpkgs.python3.withPackages (ps: [ ps.jupyterlab ]);

  rtsopts = "-M3g -N2";

  ihaskellJupyterCmdSh = cmd: extraArgs: nixpkgs.writeScriptBin "ihaskell-${cmd}" ''
    #! ${nixpkgs.stdenv.shell}
    export GHC_PACKAGE_PATH="$(echo ${ihaskellEnv}/lib/*/package.conf.d| tr ' ' ':'):$GHC_PACKAGE_PATH"
    export R_LIBS_SITE=${builtins.readFile r-libs-site}
    export PATH="${nixpkgs.stdenv.lib.makeBinPath ([ ihaskellEnv jupyterlab ] ++ systemPackages nixpkgs)}''${PATH:+:}$PATH"
    ${ihaskellEnv}/bin/ihaskell install \
      -l $(${ihaskellEnv}/bin/ghc --print-libdir) \
      --use-rtsopts="${rtsopts}" \
      && ${jupyterlab}/bin/jupyter ${cmd} ${extraArgs} "$@"
  '';
in
nixpkgs.buildEnv {
  name = "ihaskell-with-packages";
  buildInputs = [ nixpkgs.makeWrapper ];
  paths = [ ihaskellEnv jupyterlab ];
  postBuild = ''
    ln -s ${ihaskellJupyterCmdSh "lab" ""}/bin/ihaskell-lab $out/bin/
    ln -s ${ihaskellJupyterCmdSh "notebook" ""}/bin/ihaskell-notebook $out/bin/
    ln -s ${ihaskellJupyterCmdSh "nbconvert" ""}/bin/ihaskell-nbconvert $out/bin/
    ln -s ${ihaskellJupyterCmdSh "console" "--kernel=haskell"}/bin/ihaskell-console $out/bin/
    for prg in $out/bin"/"*;do
      if [[ -f $prg && -x $prg ]]; then
        wrapProgram $prg --set PYTHONPATH "$(echo ${jupyterlab}/lib/*/site-packages)"
      fi
    done
  '';
}
