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
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

/// This main program provides a command-line interface to run a series of sequential grating efficiency calculations.
/// The results are written to an output file, and (optionally) a second file is written to provide information on the
/// status of the calculation. This file is only responsible for input processing and output; all numerical details
/// are structured within Grating and PESolver.

int main(int argc, char** argv) {

	CommandLineOptions io;
	if(!io.parseFromCommandLine(argc, argv)) {
		std::cerr << "Invalid command-line options: " << io.firstErrorMessage() << std::endl;
		return -1;
	}

	if(io.showLegal) {
		std::cout 
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
		std::cout 
		"PEG Copyright (C)\n"
		"2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de)\n"
		"2012 Mark Boots (mark.boots@usask.ca)\n"
		"This free software comes with ABSOLUTELY NO WARRANTY; you are welcome to redistribute it under certain conditions; run with --showLegal for details.\n";
	}

	// Open the output file:
	std::ofstream outputFile(io.outputFile.c_str(), std::ios::out | std::ios::trunc);
	if(!outputFile.is_open()) {
		std::cerr << "Could not open output file " << io.outputFile << std::endl;
		return -1;
	}

	// Check that we can open the progress file, if provided:
	if(!io.progressFile.empty()) {
		std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
		if(!progressFile.is_open()) {
			std::cerr << "Could not open progress file " << io.progressFile << std::endl;
			return -1;
		}
	}

	// Write the file header:
	writeOutputFileHeader(outputFile, io);
	// Remember this position in the output file; it is where we will (repeatedly) write the result lines.
	// outputFile itself never gets a progress entry -- use --progressFile for that.
	std::streampos outputFilePosition = outputFile.tellp();

	// How many steps do we have?
	int totalSteps = int((io.max - io.min)/io.increment) + 1;

	// Write the initial progress to progressFile, if provided:
	if(!io.progressFile.empty()) {
		std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
		writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
	}

	// create the grating object.
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
		grating = std::make_unique<CustomProfileGrating>(io.period, io.geometry, io.material, io.coating, io.coatingThickness);
		break;
	default:
		// this should never happen; input validation assures one of the valid grating types.
		std::cerr << "Unknown Grating-Profil." << std::endl;
		return -1;
	}

	// set math options: truncation index from input.
	MathOptions mathOptions(io.N, io.integrationTolerance);

	// output data stored here:
	bool anyFailures = false;
	bool anySuccesses = false;
	std::vector<TEResult> resultsTE;
	std::vector<TMResult> resultsTM;
	std::vector<TEResult> resultsCombined; // reuses TEResult's layout for the averaged spectrum

	// sequential loop over calculation steps
	for(int i=0; i<totalSteps; ++i) {

		double currentValue = io.min + io.increment*i;

		// determine wavelength (um): depends on mode and eV/um setting.
		double wavelength = (io.mode == CommandLineOptions::ConstantWavelength) ? io.wavelength : currentValue;
		// interpret input wavelength as eV instead, and convert to actual wavelength.
		// Formula: wavelength = hc / eV. => hc = 1.23984172 eV * um.
		if(io.eV)
			wavelength = M_HC / wavelength;

		// determine incidence angle: depends on mode and possibly wavelength.
		double incidenceAngle;
		switch(io.mode) {
		case CommandLineOptions::ConstantIncidence:
			incidenceAngle = io.incidenceAngle;
			break;
		case CommandLineOptions::ConstantIncludedAngle: {
			double ciaRad = io.includedAngle * M_PI / 180;
			// formula for constant included angle: satisfies alpha + beta = cia, and grating equation io.toOrder*wavelength/d = sin(beta) - sin(alpha).
			incidenceAngle = (asin(-io.toOrder*wavelength/2/io.period/cos(ciaRad/2)) + ciaRad/2) * 180 / M_PI;
			break;
			}
		case CommandLineOptions::ConstantWavelength:
			incidenceAngle = currentValue;
			break;
		default:
			// never happens: input validation assures valid mode.
			incidenceAngle = 0;
			break;
		}

		// run the calculation(s), depending on which polarizations are requested.
		TEResult resultTE = TEResult(TEResult::InactiveCalculation);
		TMResult resultTM = TMResult(TMResult::InactiveCalculation);

		if(io.computeTE || io.combineTETM)
			resultTE = grating->getEffTE(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, io.printDebugOutput, io.threads, io.measureTiming);
		if(io.computeTM || io.combineTETM)
			resultTM = grating->getEffTM(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, io.printDebugOutput, io.threads, io.measureTiming);

		if(io.computeTE) {
			resultsTE.push_back(resultTE);
			if(resultTE.status == TEResult::Success) anySuccesses = true; else anyFailures = true;
		}
		if(io.computeTM) {
			resultsTM.push_back(resultTM);
			if(resultTM.status == TMResult::Success) anySuccesses = true; else anyFailures = true;
		}
		if(io.combineTETM && resultTE.status == TEResult::Success && resultTM.status == TMResult::Success) {
			TEResult combined = resultTE;
			for(std::size_t k = 0; k < combined.eff.size(); ++k)
				combined.eff[k] = (resultTE.eff[k] + resultTM.eff[k]) / 2.0;
			resultsCombined.push_back(combined);
		}

		// Print results to output file (no progress written here -- see --progressFile).
		outputFile.seekp(outputFilePosition);
		outputFile << "# Output" << std::endl;
		if(io.computeTE){
			outputFile << "# TE" << std::endl;
			for(const auto& result : resultsTE)
				writeOutputFileResult(outputFile, result, io);
		}
		if(io.computeTM){
			outputFile << "# TM" << std::endl;
			for(const auto& result : resultsTM)
				writeOutputFileResult(outputFile, result, io);
		}
		if(io.combineTETM){
			outputFile << "# (TE+TM)/2" << std::endl;
			for(const auto& result : resultsCombined)
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