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
#include <fstream>
#include <string>
#include <getopt.h>
#include <cfloat>
#include <climits>
#include <vector>
#include <cmath>

/// This file contains common support routines for both the parallel and sequential versions of the PEG command-line app, related to input and output handling.

/// This structure specifies the operating input, and is able to parse it from command-line options \c argc, \c argv.
class CommandLineOptions {
public:
	enum Mode {InvalidMode, ConstantIncidence, ConstantIncludedAngle, ConstantWavelength};
	
	// Input variables:
	////////////////////////////////
	Mode mode;

	double min, max, increment, incidenceAngle, includedAngle, wavelength;
	int toOrder;

	int N;
	double integrationTolerance;

	Grating::Profile profile;
	double period;
	std::vector<double> geometry;
	std::string material, coating;
	double coatingThickness;

	std::string outputFile, progressFile;

	bool eV;
	bool printDebugOutput;
	int threads;
	bool measureTiming;

	bool showLegal;

	double rmsRoughnessNm;

	// Polarization calculation options
	bool computeTE;
	bool computeTM;
	bool combineTETM;

	////////////////////////////////
	
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

/// Combines an already-computed TE result and TM result into a single result, using the simple average (TE_eff + TM_eff) / 2 per diffraction order.
/// Does not run any new calculation itself; \c teResult and \c tmResult must come from Grating::getEffTE() / Grating::getEffTM() for the same incidence/wavelength/N.
/// \pre Callers must ensure both teResult.status and tmResult.status are Result::Success before calling this.
Result getEffCombined(const Result& teResult, const Result& tmResult);

/// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& outputFileStream, const CommandLineOptions& io);

/// This helper function appends a single efficiency result to the output file stream
void writeOutputFileResult(std::ostream& outputFileStream, const Result& result, const CommandLineOptions& io);

/// This helper function appends the progress description to the given output stream
void writeOutputFileProgress(std::ostream& outputFileStream, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures);

#endif
