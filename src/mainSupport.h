/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#ifndef PEGMAINSUPPORT_H
#define PEGMAINSUPPORT_H

#include "PEG.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>

/*
This file contains the routines shared by the parallel (mainMPI.cpp) and the sequential (mainSerial.cpp) command-line app. It has three jobs:
  1. input:   'CommandLineOptions' parses, validates commandline; 'checkMaterialData()' checks the refractive-index tables.
  2. sweep:   numSweepSteps(), sweepPoint(), computeSweepPoint(), makeGrating() turn the options into the calculation of one sweep point.
  3. output:  printLegalNotice() and the writeOutputFile*() functions.
*/

/// This class holds all operating input, and is able to parse and validate it from the command-line options \c argc, \c argv.
class CommandLineOptions {
public:
	enum Mode {InvalidMode, ConstantIncidence, ConstantIncludedAngle, ConstantWavelength};
	
	/// Input variables ("not given" is marked with DBL_MAX for doubles and INT_MAX for integers):
	Mode mode;

	/// The swept range: wavelength (um) or energy (eV, if \c eV) in the modes ConstantIncidence / ConstantIncludedAngle; incidence angle (deg) in ConstantWavelength mode.
	double min, max, increment;
	/// Fixed quantity of the mode: incidence angle (deg), included angle (deg) or wavelength (um or eV, see \c eV).
	double incidenceAngle, includedAngle, wavelength;
	/// Diffraction order the monochromator is tuned to (only ConstantIncludedAngle mode).
	int toOrder;

	/// Truncation index N: the calculation uses the orders -N..+N.
	int N;
	/// Tolerance of the numerical integration of the trial solutions.
	double integrationTolerance;

	Grating::Profile profile;
	double period;
	/// Profile parameters, meaning depends on \c profile
	std::vector<double> geometry;
	std::string material, coating;
	double coatingThickness;

	/// \c progressFile is optional; \c outputFile is required
	std::string outputFile, progressFile;

	/// true: \c min, \c max, \c increment (and a fixed \c wavelength) are given in eV instead of um.
	bool eV;
	/// flagg for debug
	bool printDebugOutput;
	/// Number of threads for the fine parallelization inside one calculation.
	int threads;
	bool measureTiming;

	bool showLegal;

	double rmsRoughnessNm;

	/// Which polarizations are calculated. At least one must be set (--TE, --TM or --TETM); each one is written as its own section in the output file.
	bool TE;
	bool TM;

	
	/// Default constructor initializes all input variables to recognizable values. Doubles are set to DBL_MAX, and integers are set to INT_MAX.
	CommandLineOptions() {
		init();
	}
	
	/// This constructor parses immediately from the command-line input arguments. Check for validity after with isValid().
	CommandLineOptions(int argc, char** argv) {
		init();
		parseFromCommandLine(argc, argv);
	}
	
	/// Sets options based on command-line input arguments. Returns isValid().
	bool parseFromCommandLine(int argc, char** argv);
	
	/// Returns true if all required input options have been provided correctly. Sets firstErrorMessage() if any problems found.
	bool isValid();
	
	/// If !isValid(), returns a description of the first error found during validation.
	std::string firstErrorMessage() const { return firstErrorMessage_; }
	
protected:
	/// Initializes all input variables to recognizable values. Doubles are set to DBL_MAX, and integers are set to INT_MAX.
	void init();
	/// A description of the first error found during validation.
	std::string firstErrorMessage_;
};

/// Number of points in the sweep min..max in steps of increment, robust against floating-point truncation.
int numSweepSteps(const CommandLineOptions& io);

/// Checks that refractive-index data exists for the grating and coating material across the whole sweep. On failure, error describes the problem.
bool checkMaterialData(const CommandLineOptions& io, std::string& error);

/// Creates the grating described by the options. Returns nullptr for an unknown profile (cannot happen after a successful CommandLineOptions::isValid()).
std::unique_ptr<Grating> makeGrating(const CommandLineOptions& io);

/// Hands io.threads to the libraries behind the solvers (Eigen). The OpenMP loops get the thread count directly from the solver call. Call once, before the first calculation.
void applyThreadSettings(const CommandLineOptions& io);

/// Wavelength (um) and incidence angle (deg) of sweep point \c step (0 <= step < numSweepSteps(io)), depending on the mode and on eV / um.
/// Returns false if there is no usable geometry for this point (non-finite angle, or |angle| >= 90 deg). This happens in ConstantIncludedAngle mode when the order \c toOrder does not exist at this wavelength. The outputs are set in either case.
bool sweepPoint(const CommandLineOptions& io, int step, double& wavelength, double& incidenceDeg);

/// The requested polarization results of one sweep point; the one not requested stays at Result::InactiveCalculation.
struct PointResults {
	Result te = Result(Result::InactiveCalculation);
	Result tm = Result(Result::InactiveCalculation);
};

/// Calculates sweep point \c step for the polarizations requested in io.TE / io.TM. Never throws: a point without usable geometry or with an exception in the solver becomes Result::OtherFailure (with wavelength and incidence set, so the 'Error:' row has its x-value).
/// \c printDebugOutput is a parameter (and not read from \c io), because the MPI version only wants it on rank 0.
PointResults computeSweepPoint(const Grating& grating, const CommandLineOptions& io, const MathOptions& mathOptions, int step, bool printDebugOutput);

/// Prints the copyright notice; \c full = the complete legal text (--showLegal), otherwise the short one-line version.
void printLegalNotice(std::ostream& os, bool full);

/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& outputFileStream, const CommandLineOptions& io);

/// This helper function appends a single efficiency result (one line) to the output file stream
void writeOutputFileResult(std::ostream& outputFileStream, const Result& result, const CommandLineOptions& io);

/// Rewrites the whole result part ('# Output', then the sections '# TE' and '# TM' if requested, one line per point) starting at stream position \c position, which is the position directly behind the header. Flushes at the end.
void writeOutputFileResults(std::ostream& outputFileStream, std::streampos position, const std::vector<Result>& resultsTE, const std::vector<Result>& resultsTM, const CommandLineOptions& io);

/// This helper function appends the progress description to the given output stream
void writeOutputFileProgress(std::ostream& outputFileStream, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures);

#endif
