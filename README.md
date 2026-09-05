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
PEG no longer depends on the GNU Scientific Library (GSL). Linear algebra is provided by [Eigen](https://eigen.tuxfamily.org/), and ODE integration by [Boost.Odeint](https://www.boost.org/doc/libs/release/libs/numeric/odeint/). Both are header-only, so no separate library build or linking step is required — just make sure the headers are reachable on your include path.

```
sudo apt install libeigen3-dev libboost-dev
```

## Building
The 'qmake' build tool from the Qt Framework can be used to generate Makefiles based on PEG_mac.pro or PEG_ubuntu.pro.  (The Qt library is not required to build PEG.)  Edit one of these files to define the include paths for Eigen (and Boost, if not already on your system include path) for your system.

This will build the single-machine command-line program 'pegSerial':

```
qmake PEG_ubuntu.pro
make
```

The pegSerial application can take advantage of fine-grained parallelization to use more than one thread (for example, as many threads as CPU cores on your machine). This speeds up the calculation of a single efficiency data point.

There is also an application to exploit coarse-grained parallelization over an arbitrary number of nodes in a cluster or grid computer, using MPI:
```
sudo apt install build-essential openmpi-bin libopenmpi-dev
```
This speeds up the calculation of many efficiency data points. To build the pegMPI program, create a makefile based on src/Makefile.example.

```
cp src/Makefile.example src/Makefile
cd src
make pegMPI
```

## Running
Single-computer version 'pegSerial':

```
./pegSerial --mode constantIncidence --min 100 --max 200 --increment 5 --eV --incidenceAngle 87 --N 15 --gratingType blazed --gratingMaterial Pt --gratingPeriod 1.2 --gratingGeometry 1.2,30 --outputFile results.txt --progressFile progress.txt --computeTE
```

Cluster-version:

```
mpiexec -n <number of nodes> ./pegMPI --mode constantIncidence --min 100 --max 200 --increment 5 --eV  --incidenceAngle 87 --N 15 --gratingType blazed --gratingMaterial Pt --gratingPeriod 1.2 --gratingGeometry 1.2,30 --outputFile results.txt --progressFile progress.txt --computeTE
```

### Command-line options

#### Grating specification:
```
--gratingType <rectangular|blazed|sinusoidal|trapezoidal|custom>
--gratingPeriod <grating period in um>
--gratingGeometry <command-delimited list of geometry parameters, in um and/or degrees>
```
- **Rectangular profile**: depth (um),valley width (um)
- **Blazed profile**: blaze angle (deg),anti-blaze angle (deg)
- **Sinusoidal profile**: depth (um)
- **Trapezoial profile**: depth (um),valley width (um),blaze angle (deg),anti-blaze angle (deg)
- **Custom profile**: depth (um), followed by list of (x,y) points that describe the surface profile y = g(x).  [x followed by y for each]. The points must go from (0,0) to (1,0) inclusive, and will be scaled so that (0,0)->(1,1) maps to (0,0)->(period,depth). For example, an isosceles triangular profile with a depth of 0.1um would need three points; the geometry argument would be: {0.1,0,0,0.5,1,1,0}.
```
--gratingMaterial <grating substrate material>
	This should be a name corresponding to a refractive index database filename, ex: Au, Ni, C, SiO2, etc.
	
--N <truncation index>
	Specifies the number of positive and negative orders to include in the Fourier expansion.
	Will also determine the number of orders that are calculated, although if you only want to calculate 3 orders, you will still need a much larger truncation index for accurate results.
```

#### Operating mode

```
--mode <constantIncidence|constantIncludedAngle|constantWavelength>
--min <min>
--max <max>
--increment <increment>

[Required, depending on the \c mode]

--incidenceAngle <incidence angle in degrees>
--includedAngle <deviation angle in degrees> --toOrder <diffraction order for the included angle>
--wavelength <wavelength in um>
```

- In **constant incidence** mode, a calculation is performed for wavelengths from --min to --max in steps of --increment, at a fixed incidence angle given by --incidenceAngle.
- In **constant included angle** mode, the incidence angle is calculated at each wavelength to ensure a constant included angle of --includedAngle between the incident light and the order specified in --toOrder. This is the operating mode for many monochromators. (Inside orders are negative, outside orders are positive.)
- In **constant wavelength** mode, a calculation is performed for incidence angles from --min to --max in steps of --increment, for a fixed wavelength given by --wavelength.

#### Output
```
--outputFile <file name>
	The calculation output will be written to this file.

--computeTE | --computeTM | --combineTETM
	If included, computes and outputs the TE-polarized efficiency (section "# TE" in the output file) || TM-polarized efficiency (section "# TM" in the output file).
	... computes both TE and TM internally and outputs their simple average as well ( $$(TE_eff + TM_eff)/2$$ per order -- section "# (TE+TM)/2" in the output file). Only points where both the TE and TM calculation succeeded are included. At least one of --computeTE, --computeTM, or --combineTETM must be specified.
```
#### Optional
```
--progressFile <file name>
	If provided, the current status of the calculation will be written in this file; it can be monitored to determine the progress of long calculations.  This provides an interface for other processes to monitor the status of this calculation (for example, a web-based or GUI front-end, etc.).
	
--eV
	If this flag is included, all wavelength inputs (--min, --max, --increment, --wavelength) will instead be interpreted as photon energies in electron volts (eV).

--coatingThickness <thickness in um>
	If provided, creates a layer of --coatingMaterial on top of the basic grating profile, by translating the profile vertically by coatingThickness um. Used to model overcoated, oxidized, or dielectric gratings.

--coatingMaterial <coating material>
	Required if a non-zero --coatingThickness is provided. This should be a name corresponding to a refractive index database filename, ex: Au, Ni, C, MgF2, etc.
	
--printDebugOutput
	If this flag is included, each calculation will print intermediate results to standard output.

--measureTiming
	If this flag is included, the solver will report the time required for each category of operations to standard output.

--integrationTolerance <tolerance>
	If provided, specifies the error tolerance (eps) required at each step of the numerical integration process. Default if not provided is 1e-5.

--showLegal
	If included, prints the full copyright and license notice (including the disclaimer about modifications) instead of the short version.
```

#### Example Output
```
# Input
mode=constantIncidence
incidenceAngle=88
units=eV
min=100
max=300
increment=5
gratingType=blazed
gratingPeriod=1.6
gratingGeometry=3.2,30.0
gratingMaterial=Au
N=5
integrationTolerance=1e-5
# Progress
status=succeeded     (inProgress, someFailed, allFailed, succeeded)
completedSteps=41
totalSteps=41
# Output
# TE
100[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
105[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
110[tab]<e-5>,<e-4>,<e-3>,<e-2>,<e-1>,<e0>,<e1>,<e2>,<e3>,<e4>,<e5>
...
```

The Output table lists reflected efficiencies at each (wavelength/eV/incidence angle) sequentially from the -N order to the +N order.  (Efficiencies are 0 if the orders are evanescent instead of propagating.)  Note that we use the sign convention where _inside diffraction orders_ are negative (n < 0), corresponding to the grating equation:

$$ \sin\beta = \sin\alpha + \frac{n \lambda}{d}$$

License
---
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

PEG is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3. http://www.gnu.org/licenses/gpl.html

This reworked version (2026) contains substantial modifications and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently. Run any PEG executable with `--showLegal` to see the full notice.

As a modification of GPLv3-licensed software, this reworked version remains bound by the same copyleft terms, independent of the underlying numerical libraries used (previously GSL, now Eigen and Boost.Odeint).


References
---

Credit first and foremost to Mark Boots — the vast majority of the effort lies in his original implementation, [markboots/peg](https://github.com/markboots/peg)!

1. M Nevière et al, Nouvelle Revue d'Optique 5 65 (1974) http://dx.doi.org/10.1088/0335-7368/5/2/301

2. L. Li, JOSA A, Vol. 13, Issue 5, pp. 1024-1035 (1996)
http://dx.doi.org/10.1364/JOSAA.13.001024

3. B.L. Henke, E.M. Gullikson, and J.C. Davis. X-ray interactions: photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, Z=1-92, Atomic Data and Nuclear Data Tables Vol. 54 (no.2), 181-342 (July 1993).

4. S.K. Sinha, E.B. Sirota, S. Garoff, and H.B. Stanley. X-ray and neutron scattering from rough surfaces. Physical Review B, Vol. 38, Issue 4, 2297 (1988), Equation 4.34. http://dx.doi.org/10.1103/PhysRevB.38.2297
	Used for the optional RMS-roughness correction factor (`--rmsRoughnessNm`).
