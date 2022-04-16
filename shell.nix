let
  pkgs = import <nixpkgs> {};
  fhs = pkgs.buildFHSUserEnv {
    name = "cuda-env";
    targetPkgs = pkgs: with pkgs; [ 
      # Enable opencl through cuda
      cudatoolkit
      linuxPackages.nvidia_x11
      stdenv.cc
      binutils

      zsh
    ];
    runScript = "zsh";
    multiPkgs = pkgs: with pkgs; [ zlib ];
    profile = ''
       export CUDA_PATH=${pkgs.cudatoolkit}
       # export LD_LIBRARY_PATH=${pkgs.linuxPackages.nvidia_x11}/lib
       export EXTRA_LDFLAGS="-L/lib -L${pkgs.linuxPackages.nvidia_x11}/lib"
       export EXTRA_CCFLAGS="-I/usr/include"
       export PATH=/home/carl/scratch/thm/bins:$PATH
       export JULIA_NUM_THREADS=6
       export JULIA_LOCAL="startup-julia.jl"
     '';
  };
in pkgs.mkShell {
  nativeBuildInputs = [fhs];
  buildInputs = [pkgs.julia-bin];
  shellHook = ''
       exec cuda-env
   '';
}
