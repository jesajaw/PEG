# About
PEG ("Parallel Efficiency of Gratings") is a tool for calculating the efficiency of diffraction gratings, particularly those used in the Soft X-ray regime.  Developed by the Materials Research Group in the Department of Physics at the University of Saskatchewan (http://beamteam.usask.ca), it was used to characterize the optical components for the REIXS XES beamline at the Canadian Light Source. It implements the _differential theory_ developed by Neviere, Vincent, and Petit [1], and is updated for stability using the _S-matrix formulation_ of Li [2].

This repository is a 2026 rework of the original PEG project by Jesaja Weintritt. The core physics and algorithm (differential theory + S-matrix formulation) are unchanged, but the codebase has been substantially modernized: the GSL dependency has been replaced with Eigen and Boost.Odeint (see Dependencies below), naming and class structure were cleaned up, and support for Transverse Magnetic (TM) polarization is being added alongside the existing Transverse Electric (TE) solver (see Limitations below). This reworked version has not been independently verified against the original; see License.


## Features
- Standard grating shape profiles: rectangular, blazed (triangular), trapezoidal, sinusoidal
- Custom grating shape profiles: any point-wise defined profile.
- Coatings: an optional coating layer on top of the grating (user-defined thickness; thick or inter-penetrating)
- Automatic lookup of complex refractive indexes for common materials from Henke data [3].
- Several built-in scanning modes:
	- over wavelength,
	- over incidence angle,
	- over wavelength maintaining constant deviation ("monochromator mode")
- Support for scalable parallel calculation:
	- fine parallelization of a single calculation, using OpenMP
	- coarse parallelization of many efficiency points, using MPI on grid computers.


## Limitations / Todo
- Transverse Electric (TE) polarization is fully supported - but Transverse Magnetic (TM) polarization is not fully implemented yet.
- Coatings:
	- The custom grating profile does not yet support coatings.
	- only single layer coatings

## Dependencies
PEG uses [Eigen](https://eigen.tuxfamily.org/) for linear algebra and [Boost.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/) for ODE integration. It no longer depends on the GNU Scientific Library (GSL). Both libraries are header-only and bundled in `third_party/`, so no separate installation is needed. They remain under their own licenses (MPL-2.0 and BSL-1.0), see `third_party/LICENSES`.

## Building
PEG is built with the included `Makefile`. You need a `C++17` compiler with OpenMP (`sudo apt install build-essential`) and, for `pegMPI`, an MPI implementation (`sudo apt install openmpi-bin libopenmpi-dev`).

```bash
make pegSerial   # single-computer version
make pegMPI      # cluster version
make all         # both
make check       # quick smoke self-test
make clean       # remove build files and binaries
```

| Option | Effect |
|---|---|
| `NATIVE=1` | Tune for the CPU you build on. |
| `DEBUG=1` | Debug build (`-O0 -g`); the default is an optimized release build (`-O3 -DNDEBUG`). |
| `EIGEN_INC=<dir>`, `BOOST_INC=<dir>` | Header locations other than `third_party/`. |

Run `make clean` after changing options and rebuild it, since object files don't do that on their own.

## Running

```bash
# single computer
./pegSerial --mode constantIncidence --min 100 --max 200 --increment 5 --eV --incidenceAngle 87 \
            --N 15 --gratingType blazed --gratingMaterial Pt --gratingPeriod 1.2 --gratingGeometry 1.2,30 \
            --outputFile results.txt --progressFile progress.txt --TE --threads 2

# cluster: same options, started through MPI
mpiexec -n <processes> ./pegMPI <same options>
```

### Command-line options

#### Grating
```text
--gratingType <rectangular|blazed|sinusoidal|trapezoidal|custom>
--gratingPeriod <um>
--gratingMaterial <name>        refractive-index file: Au, Ni, C, SiO2, ... (about 30 eV - 30 keV)
--N <n>                         truncation index: orders -N ... +N in the Fourier expansion
--gratingGeometry <v1,v2,...>   comma-separated, depends on --gratingType (um / deg):
    rectangular                 depth, valleyWidth
    blazed                      blazeAngle, antiBlazeAngle
    sinusoidal                  depth
    trapezoidal                 depth, valleyWidth, blazeAngle, antiBlazeAngle
    custom                      depth, x1,y1, x2,y2, ... the points go from (0,0) to (1,0) inclusive
				and are scaled so that (0,0)->(1,1) maps to (0,0)->(period,depth).
				An isosceles triangle of 0.1um depth is `--gratingGeometry 0.1,0,0,0.5,1,1,0`.
```

#### Operating mode
```text
--mode constantIncidence        --incidenceAngle <deg>                    # sweeps wavelength
--mode constantIncludedAngle    --includedAngle <deg> --toOrder <n>       # sweeps wavelength
--mode constantWavelength       --wavelength <um>                         # sweeps incidence angle (deg)
--min <v> --max <v> --increment <v>
```

`constantIncludedAngle` ("monochromator mode") sets the incidence angle at each wavelength so that the included angle between the incident beam and order `--toOrder` stays constant (inside orders negative, outside positive).

#### Output
```text
--outputFile <file>
--TE | --TM | --TETM = --TMTE   polarization(s) to compute (--TETM = both); at least one is required
```

#### Optional
```text
--progressFile <file>           status of a running calculation, for monitoring or front-ends
--eV                            read --min/--max/--increment/--wavelength as photon energies in eV
--coatingMaterial <name>        coating on top of the profile (requires --coatingThickness)
--coatingThickness <um>         vertical shift of the profile (requires --coatingMaterial)
--integrationTolerance <eps>    error tolerance per integration step (default 1e-5)
--rmsRoughnessNm <nm>           surface roughness correction (not verified!)
--threads <n>                   OpenMP threads (default 1)
--printDebugOutput              print intermediate results to stdout
--measureTiming                 print the time per operation category to stdout
--showLegal                     print the full copyright and license notice
```

#### runtime options
PEG parallelizes on two levels, which can be combined:

| Level | How | What it does |
|---|---|---|
| Fine (OpenMP) | `--threads <n>` | Uses several cores of one computer for a single calculation. Works in `pegSerial` and `pegMPI`. |
| Coarse (MPI) | `mpiexec -n <p> ./pegMPI` | Runs `p` processes that each calculate different points of the scan; the processes can be spread over several cluster nodes. |

```text
./pegSerial --threads 8 ...              # one computer, 8 cores
./pegSerial --threads $(nproc) ...       # all hardware threads of this computer
mpiexec -n 4 ./pegMPI --threads 2 ...    # 4 processes x 2 threads = 8 cores
```

#### Example Output
Output file (`--outputFile`), columns separated by tabs; a failed point is written as `<x>[tab]Error:<reason>`:
```text
# Input
mode=constantIncidence
incidenceAngle=87
units=eV
min=100
max=200
increment=5
gratingType=blazed
gratingPeriod=1.2
gratingGeometry=1.2,30
gratingMaterial=Pt
coatingMaterial=[none]
coatingThickness=0
N=15
integrationTolerance=1e-05
# Output
# TE
100[tab]<e-5>[tab]<e-4>[tab]<e-3>[tab]<e-2>[tab]<e-1>[tab]<e0>[tab]<e1>[tab]<e2>[tab]<e3>[tab]<e4>[tab]<e5>
105[tab]<e-5>[tab]<e-4>[tab]<e-3>[tab]<e-2>[tab]<e-1>[tab]<e0>[tab]<e1>[tab]<e2>[tab]<e3>[tab]<e4>[tab]<e5>
...
```
Progress file (`--progressFile`):
```text
# Progress
status=succeeded     (inProgress, someFailed, allFailed, succeeded)
completedSteps=41
totalSteps=41
```

The Output table lists reflected efficiencies at each (wavelength/eV/incidence angle) sequentially from the -N order to the +N order.  (Efficiencies are 0 if the orders are evanescent instead of propagating.)  Note that we use the sign convention where _inside diffraction orders_ are negative (n < 0), corresponding to the grating equation:

$$ \sin\beta = \sin\alpha + \frac{n \lambda}{d}$$


## License
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

PEG is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3. http://www.gnu.org/licenses/gpl.html

This reworked version (2026) contains substantial modifications and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently. Run any PEG executable with `--showLegal` to see the full notice.

As a modification of GPLv3-licensed software, this reworked version remains bound by the same copyleft terms, independent of the underlying numerical libraries used (previously GSL, now Eigen and Boost.Odeint).


## References
Credit first and foremost to Mark Boots -- the vast majority of the effort lies in his original implementation, [markboots/peg](https://github.com/markboots/peg)!

1. M Nevière et al, Nouvelle Revue d'Optique 5 65 (1974) http://dx.doi.org/10.1088/0335-7368/5/2/301

2. L. Li, JOSA A, Vol. 13, Issue 5, pp. 1024-1035 (1996) http://dx.doi.org/10.1364/JOSAA.13.001024

3. B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92, Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993). https://henke.lbl.gov/optical_constants/getdb2.html

4. S.K. Sinha, E.B. Sirota, S. Garoff, and H.B. Stanley. X-ray and neutron scattering from rough surfaces. Physical Review B, Vol. 38, Issue 4, 2297 (1988), Equation 4.34. http://dx.doi.org/10.1103/PhysRevB.38.2297
	Used for the optional RMS-roughness correction factor (`--rmsRoughnessNm`). TODO