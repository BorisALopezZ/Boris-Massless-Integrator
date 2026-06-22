# Boris-Massless-Integrator

Structure-preserving Boris-type integrator for massless charged particle dynamics

## Description

This repository contains the implementation accompanying the manuscript

> **A Structure-preserving Boris-type integrator for massless charged particle dynamics**

submitted to *Computer Physics Communications*.

## Contents

```text
src/
    boris_massless.c
    rk4.c

examples/

figures/

paper/
```

## Compilation
```bash
gcc -O3 -march=native boris_massless.c -lm -o boris_massless
```
Alternatively

```bash
gcc boris_massless.c -lm -O3 -march=native
```

## Running

```bash
./boris_massless
```
or 

```bash
./a.out
```

## Numerical tests

The implementation reproduces the numerical experiments presented in the manuscript:

- Uniform magnetic fields
- Crossed electromagnetic fields
- Plane electromagnetic waves
- Long-time simulations up to 10^7 integration steps.
- Comparisons with a classical fourth-order Runge–Kutta integrator.

## Author

Boris Alberto Lopez Zubieta

## Citation

If you use this code, please cite:

B. A. Lopez Zubieta,
*A Structure-preserving Boris-type integrator for massless charged particle dynamics*,
Computer Physics Communications (under review,2026).

## Requirements

The code was tested under Ubuntu Linux using

- GCC compiler
- Standard C library
- Linux/Unix environment

## License

MIT License

## Repository status

This repository accompanies a manuscript currently under review in Computer Physics Communications.

Future updates may include:

- Predictor-corrector extensions;
- Additional benchmark problems;
- Reproducibility scripts.

## Reproducibility

The repository contains the implementations used to generate the numerical results presented in the accompanying manuscript. Additional benchmark cases and examples may be included in future updates.

## Availability

The source code is publicly available through this GitHub repository and will be maintained and updated as necessary.
