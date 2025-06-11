 {
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs?ref=nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in {
        packages = rec {
          default = xspringies;

          xspringies = pkgs.stdenv.mkDerivation rec {
            pname = "xspringies";
            version = "1.12";

            src = ./xspringies-1.12;

            postPatch = ''
              rm -v Imakefile
              cp -v Makefile.std Makefile
              substituteInPlace Makefile \
                --replace "/bin/cp" \
                          "cp" \
               --replace "DDIR = /usr/games/" \
                          "DDIR = $out/" \
                --replace "CFLAGS = " \
                          "CFLAGS = -Wno-error=format-security "
            '';

            buildPhase = ''
              make
            '';

            installPhase = ''
              mkdir -p $out/bin
              make install
            '';

            nativeBuildInputs = with pkgs; [
              xorg.imake
              xorg.gccmakedep
            ];

            buildInputs = with pkgs; [
              xorg.libX11
            ];
          };
        };
      }
    );
}

