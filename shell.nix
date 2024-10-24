{pkgs ? import <nixpkgs> {}}:
pkgs.clangStdenv.mkDerivation {
  name = "shell";
  buildInputs = with pkgs; [
    cmake
    clang-tools

    cgal
    mpfr
    boost
  ];
}
