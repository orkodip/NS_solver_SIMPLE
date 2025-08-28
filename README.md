# NS_solver_SIMPLE
This repository contain the codes (written in C++) for incompressible Navier-Stokes solver in primitive variables based on the Semi-Implicit Method for Pressure Linked Equations (SIMPLE) algorithm and the Quadratic Upwind Interpolation for Convective Kinetics (QUICK) scheme with an improved deferred correction technique (second order upwind scheme with a third order correction) in structured non-uniform Cartesian mesh. The energy equation is not coupled with the momentum equation (forced convection phenomenon). Finite Volume Method is used to discretize the governing equations.

The code is designed to solve the fluid flow and thermal characteristics of confined laminar jet impingement cooling of discrete heat sources using nanofluids.

# Reference
Mookherjee, O., Pramanik, S. and Kar, U.K., 2020. Numerical investigation of a confined laminar jet impingement cooling of heat sources using nanofluids. Journal of Heat Transfer, 142(8), p.082301.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ jet.cpp -o output_flow
g++ energy.cpp -o output_en
```
To run this code

```bash
./output_flow
./output_en
```
