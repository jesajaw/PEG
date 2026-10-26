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

	// set threads
	Eigen::setNbThreads(io.threads);

	// if license-flagg included, print full legal
	if(io.showLegal) {
		std::cout <<
		"Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de)\n"
		"and 2012 Mark Boots (mark.boots@usask.ca).\n\n"

		"This program was originally implemented as a part of the Parallel Efficiency of Gratings project (\"PEG\")\n"
		"and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it\n"
		"under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.\n"
		"See http://www.gnu.org/licenses/gpl.html for details.\n\n"

		"This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been\n"
		"independently verified against the original. It is distributed in the hope that it will be useful,\n"
		"but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.\n\n";
	}
	else {
		std::cout <<
		"PEG Copyright (C)\n"
		"2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de)\n"
		"2012 Mark Boots (mark.boots@usask.ca)\n"
		"This free software comes with ABSOLUTELY NO WARRANTY; you are welcome to redistribute it under certain conditions; run with --showLegal for details.\n";
	}

	// How many steps do we have? number of points of the sweep 'min'..'max' in steps of 'increment'
	int totalSteps = numSweepSteps(io);

	// create the grating object.
	///\todo shared factory function in mainSupport
	std::unique_ptr<Grating> grating;
	switch(io.profile) {
		case Grating::RectangularProfile:
			grating = std::make_unique<RectangularGrating>(io.period, io.geometry[0], io.geometry[1], io.material, io.coating, io.coatingThickness);
			break;
		case Grating::BlazedProfile:
			grating = std::make_unique<BlazedGrating>(io.period, io.geometry[0], io.geometry[1], io.material, io.coating, io.coatingThickness);
			break;
		case Grating::SinusoidalProfile:
			grating = std::make_unique<SinusoidalGrating>(io.period, io.geometry[0], io.material, io.coating, io.coatingThickness);
			break;
		case Grating::TrapezoidalProfile:
			grating = std::make_unique<TrapezoidalGrating>(io.period, io.geometry[0], io.geometry[1], io.geometry[2], io.geometry[3], io.material, io.coating, io.coatingThickness);
			break;
		case Grating::CustomProfile:
			///\todo , io.coating, io.coatingThickness
			grating = std::make_unique<CustomProfileGrating>(io.period, io.geometry, io.material);
			break;
		default:
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
	///\todo this section could go into shared sup
	for(int i=0; i<totalSteps; ++i) {

		// value of the swept quantity: wavelength/energy or incidence angle in degrees
		double currentValue = io.min + io.increment*i;

		// determine wavelength (um): fixed (--wavelength) in constant-wavelength mode, swept in other modes
		double wavelength = (io.mode == CommandLineOptions::ConstantWavelength) ? io.wavelength : currentValue;
		// interpret input wavelength as eV instead, and convert to actual wavelength.
		// Formula: wavelength = hc / eV. => hc = 1.23984172 eV * um.
		if(io.eV)
			wavelength = M_HC / wavelength;

		// determine incidence angle: depends on mode and possibly wavelength.
		double incidenceAngle;
		switch(io.mode) {
			case CommandLineOptions::ConstantIncidence: // ConstantIncidence : fixed --incidenceAngle
				incidenceAngle = io.incidenceAngle;
				break;
			case CommandLineOptions::ConstantIncludedAngle: { //   ConstantIncludedAngle : monochromator mode
				double ciaRad = io.includedAngle * M_PI / 180;
				// formula for constant included angle: satisfies alpha + beta = cia, and grating equation io.toOrder*wavelength/d = sin(beta) - sin(alpha).
				///\todo add isinfinit => shared sup
				incidenceAngle = (asin(-io.toOrder*wavelength/2/io.period/cos(ciaRad/2)) + ciaRad/2) * 180 / M_PI;
				break;
				}
			case CommandLineOptions::ConstantWavelength: // ConstantWavelength : the swept value itself is the incidence angle
				incidenceAngle = currentValue;
				break;
			default:
				// never happens: input validation assures valid mode.
				incidenceAngle = 0;
				break;
		}

		// results start as 'InactiveCalculation' (= not requested); a requested polarization overwrites it.
		Result resultTE = Result(Result::InactiveCalculation);
		Result resultTM = Result(Result::InactiveCalculation);

		// one complete solve: efficiencies of all orders -N..+N for this wavelength and angle.
		if(io.TE)
			resultTE = grating->getEffTE(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, io.printDebugOutput, io.threads, io.measureTiming);
		if(io.TM)
			resultTM = grating->getEffTM(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, io.printDebugOutput, io.threads, io.measureTiming);

		// failed points are kept too
		if(io.TE) {
			resultsTE.push_back(resultTE);
			if(resultTE.status == Result::Success) anySuccesses = true; else anyFailures = true;
		}
		if(io.TM) {
			resultsTM.push_back(resultTM);
			if(resultTM.status == Result::Success) anySuccesses = true; else anyFailures = true;
		}

		// Print results to output file
		outputFile.seekp(outputFilePosition);
		outputFile << "# Output" << std::endl;
		if(io.TE){
			outputFile << "# TE" << std::endl;
			for(const auto& result : resultsTE)
				writeOutputFileResult(outputFile, result, io);
		}
		if(io.TM){
			outputFile << "# TM" << std::endl;
			for(const auto& result : resultsTM)
				writeOutputFileResult(outputFile, result, io);
		}

		// Update progress in progressFile, if provided.
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			writeOutputFileProgress(progressFile, i+1, totalSteps, anySuccesses, anyFailures);
		}

	} // end of calculation loop.
	outputFile.close();
	return 0;
}