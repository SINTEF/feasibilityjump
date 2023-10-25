# Feasibility Jump

This repository contains the reference C++ implementation of the Feasibility
Jump, which is a primal heuristic for mixed-integer linear programming.

See the [open access paper](https://link.springer.com/article/10.1007/s12532-023-00234-8) for details about
the algorithm.

## Usage

To build the Feasibility Jump reference implementation with FICO Xpress
integration, Xpress must be installed and the `$XPRESSDIR` environment variable
must be set to the directory where Xpress is installed.  To build, run:

```
make
```

This should produce the executable `xpress_fj`, which solves MIP problems given
in `.mps` format. Its command line interface is:

```
Usage: xpress_fj [--timeout|t TIMEOUT] [--save-solutions|-s OUTDIR] [--verbose|-v] [--heuristic-only|-h] [--exponential-decay|-e] [--relax-continuous|-r] INFILE
```
 
For example, to use the FJ heuristic alone, and save solutions to the current
directory, run:

```
./xpress_fj -h benchmark/instances/miplib2017benchmark/academictimetablesmall.mps --save-solutions .
```

To run the integrated solver where FJ solutions are given to Xpress, run:

```
./xpress_fj --timeout 60 benchmark/instances/miplib2017benchmark/academictimetablesmall.mps --save-solutions .
```

## Benchmark

The `run_benchmark.sh` script will reproduce results used in the
[report](feasibility_jump_2022-11-07.pdf).

The script assumes that you are using Linux and that the `$XPRESSDIR`
environment variable is set. It then does the following:

 * Builds both `xpress_fj` and a baseline solver `xpress_baseline` that only
   calls the Xpress solver with standard parameters.
 * Runs each of the two solvers on 240 MPS instances from the MIPLIB 2017
   benchmark set, which needs to be placed (as `.mps` files) in the directory
   `benchmark/instances/miplib2017benchmark`. Timeouts of both 60 and 600 seconds
   are used.
 * Produces an HTML file that presents the results in `benchmark/report.html`.


