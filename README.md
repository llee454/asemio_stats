OCaml Gnu Scientific Library Readme
===================================

This package illustrates how to use OCaml's C interface to use the Gnu Scientific Library directly.

The Gnu Scientific Library defines a comprehensive and very mature collection of functions that perform numerical computations such as numerical integration, curve fitting, and linear algebra. It has excellent documentation and is very extensive compare to other libraries available within OCaml.

OCaml's C interface is well documented in the C manual.

In this example I used dune init project to create the dune, dune-project, *.opam files and the bin, test, and other directories. I then created the gsl.c file to define wrappers around some functions provided by GSL and defined gsl.ml to define the OCaml modules that expose these functions. Compilation is handled by dune. To build run `une exec vim_merlin_test`.

You will need to install the GSL library as a shared object library (se apt install libgsl-dev) to build this example.

The generated binary will call GSL to compute an integral and a linear regressor and will compute the value and parameters of each.

Initializing the Build Environment
----------------------------------

```bash
opam switch create . ocaml-variants.4.10.0+flambda --no-install
opam update
opam install --deps-only .
dune build
dune exec bin/main.exe
```

Additional MacOS Installation Instructions
------------------------------------------

If you get the following error message: `symbol not found in flat namespace '_cblas_caxpy', reinstall conf-openblas while passing `pkgconfig` the correct install path:

```bash
PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/usr/local/Cellar/openblas/0.3.21/lib/pkgconfig opam reinstall conf-openblas
```

Running Unit Tests
------------------

```bash
dune runtest
```
