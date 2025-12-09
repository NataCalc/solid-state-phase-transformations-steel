# solid-state-phase-transformations-steel
Modelling the thermodynamic equilibrium and para-equilibrium between γ and α phases (using a CALPHAD-based approach)

The access to full phd thesis: https://theses.fr/2012LORR0342

Computational Models and Simulation Code:
This repository contains the code and resources developed for research on the thermodynamic and kinetic mechanisms governing the austenite (γ) to ferrite (α) phase transformation in steels.
The project focuses on modelling interface conditions, diffusion-controlled growth, and the transition between equilibrium and para-equilibrium regimes.

Project Description:
This work studies the γ → α phase transformation using a combination of thermodynamic modelling, kinetic theory, and numerical simulation.
Key themes include:

Thermodynamic interface conditions
– equilibrium and para-equilibrium between austenite and ferrite
– CALPHAD-based evaluation of chemical potentials

Kinetic modelling of interface growth
– interface mobility and diffusion effects
– transition between mixed-mode regimes

Simulation of ferrite growth
– Fe–C binary and Fe–C–X (e.g. Mn) alloys
– diffusion-controlled and interface-controlled transformations
– exploration of how carbon redistribution affects transformation paths

Interface thickening model inspired by the literature on mixed-mode transformations and extended to multicomponent steels.

This repository provides the core numerical implementation used to simulate these mechanisms.

Build & Run Instructions:
The project is compiled using a Makefile. All .cpp files in the directory are automatically built.

Requirements:

g++ or another modern C++ compiler

pthread support

POSIX environment (Linux/macOS or WSL)

Optional: headers in /opt/local/include

1. Build the executable
make

2. Run the program
./preci

3. Clean build files
make clean

Project Structure:
src/       → C++ source files used in the simulations
Makefile   → build script (compile + link all .cpp files)
preci      → generated executable (after building)
