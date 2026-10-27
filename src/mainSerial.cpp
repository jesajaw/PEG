/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

#include "PEG.h"
#include "mainSupport.h"

#include <iostream>
#include <fstream>
#include <string>
//#include <cfloat>
//#include <climits>
#include <vector>
//#include <algorithm>
#include <cmath>
#include <memory>
#include <Eigen/Core>

/*
This main program provides a command-line interface to run a series of sequential grating efficiency calculations. The results are written to an output file and (optionally) a second progress file. This file is only responsible for input processing and output; all numerical details are structured within Grating and both Solvers.
*/

int main(int argc, char** argv) {

	// parse command line options into input variables / 'io'
	CommandLineOptions io;
	if(!io.parseFromCommandLine(argc, argv)) {
		std::cerr << "Invalid command-line options: " << io.firstErrorMessage() << std::endl;
		return -1;
	}

	// is the chosen material available?
	std::string dataError;
	if(!checkMaterialData(io, dataError)) {
		std::cerr << dataError << std::endl;
		return -1;
	}

	applyThreadSettings(io);

	// copyright notice: the full legal text with --showLegal, otherwise the short one
	printLegalNotice(std::cout, io.showLegal);

	// How many steps do we have? number of points of the sweep 'min'..'max' in steps of 'increment'
	int totalSteps = numSweepSteps(io);

	// create the grating object.
	std::unique_ptr<Grating> grating = makeGrating(io);
	if(!grating) {
		// this should never happen; input validation assures one of the valid grating types.
		std::cerr << "Unknown Grating-Profil." << std::endl;
		return -1;
	}

	// math options: numerical settings shared by all points: truncation index 'N' and ODE tolerance
	MathOptions mathOptions(io.N, io.integrationTolerance);

	// Open the output file: rewrite every step
	std::ofstream outputFile(io.outputFile.c_str(), std::ios::out | std::ios::trunc);
	if(!outputFile.is_open()) {
		std::cerr << "Could not open output file " << io.outputFile << std::endl;
		return -1;
	}
	// Write the file header: '# Input' block
	writeOutputFileHeader(outputFile, io);
	// Remember this position in the output file; it is where we will (repeatedly) write the result lines.
	std::streampos outputFilePosition = outputFile.tellp();

	// open progress file (if provided) and write the first progress entry
	///\todo also visible to extern exc?
	if(!io.progressFile.empty()) {
		std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
		if(!progressFile.is_open()) {
			std::cerr << "Could not open progress file " << io.progressFile << std::endl;
			return -1;
		}
		writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
	}

	// output data stored here:
	bool anyFailures = false;
	bool anySuccesses = false;
	std::vector<Result> resultsTE;
	std::vector<Result> resultsTM;

	// sequential loop over calculation steps
	for(int i=0; i<totalSteps; ++i) {

		// one complete solve per requested polarization: efficiencies of all orders -N..+N for this sweep point.
		PointResults point = computeSweepPoint(*grating, io, mathOptions, i, io.printDebugOutput);
		Result& resultTE = point.te;
		Result& resultTM = point.tm;

		// failed points are kept too (written as an 'Error:...' row), so every section has one row per point.
		if(io.TE) {
			resultsTE.push_back(resultTE);
			if(resultTE.status == Result::Success) anySuccesses = true; else anyFailures = true;
		}
		if(io.TM) {
			resultsTM.push_back(resultTM);
			if(resultTM.status == Result::Success) anySuccesses = true; else anyFailures = true;
		}

		// Print results to output file
		writeOutputFileResults(outputFile, outputFilePosition, resultsTE, resultsTM, io);

		// Update progress in progressFile, if provided.
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			writeOutputFileProgress(progressFile, i+1, totalSteps, anySuccesses, anyFailures);
		}

	} // end of calculation loop.
	outputFile.close();
	return 0;
}