#include "PEG.h"
#include "MainSupport.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cfloat>
#include <climits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>

#include "mpi.h"

/*
This mainstructur provides a command-line interface to run the calculations. The results are written to an output file and optionally a second file for information on the status of the calculation.
This file is only responsible for input processing and output; all physical and numerical details are structured somewhere else.
*/

int main(int argc, char** argv) {
	
	// Initialize MPI
	MPI_Init(&argc, &argv);
	
	// Wait for all processes to be up and running, and then start timing
	MPI_Barrier(MPI_COMM_WORLD);
	double startTime = MPI_Wtime();
	
	int rank, commandSize;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &commandSize);
	
	// parse command line options into input variables.
	CommandLineOptions io;

	// Prevent getopt() from printing error messages except on Process 0.
	if(rank != 0) opterr = 0;
	if(!io.parseFromCommandLine(argc, argv)) {
		if(rank == 0) std::cerr << "Invalid command-line options: " << io.firstErrorMessage() << std::endl;
		MPI_Finalize();
		return -1;
	}
	
	// Legal
	if(rank == 0) {
		if(io.showLegal) {
			std::cout 
			"Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) \n"
			"and 2012 Mark Boots (mark.boots@usask.ca).\n\n"

			"This program was originally implemented as a part of the Parallel Efficiency of Gratings project (\"PEG\")\n"
			"and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it\n"
			"under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.\n"
			"See <http://www.gnu.org/licenses/> for details.\n\n"

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
	}

	// On Process 0: Open the output file:
	std::ofstream outputFile;
	std::streampos outputFilePosition;
	if(rank == 0) {
		outputFile.open(io.outputFile.c_str(), std::ios::out | std::ios::trunc);

		if(!outputFile.is_open()) {
			std::cerr << "Could not open output file " << io.outputFile << std::endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// Check that we can open the progress file, if provided:
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			if(!progressFile.is_open()) {
				std::cerr << "Could not open progress file " << io.progressFile << std::endl;
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
		}

		// Write the file header:
		writeOutputFileHeader(outputFile, io);

		// Remember this position in the output file; it is where we will write the progress and output lines
		outputFilePosition = outputFile.tellp();
	}
	
	// How many steps do we have: int(min to max / increment)?
	int totalSteps = int((io.max - io.min)/io.increment) + 1;
	
	// On Process 0: Write the initial progress:
	if(rank == 0) {
		writeOutputFileProgress(outputFile, 0, totalSteps, false, false);
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
		}
	}
	
	// create the grating object.
	///\todo use array for coating & thickness -> multilayer
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
		if(rank == 0){
			// this should never be the case through the MPI_Abort, just wrote it to be safe
			std::cerr << "Unknown Grating-Profil." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
			return -1;
		}
	}

	// set math options: truncation index from input.
	MathOptions mathOptions(io.N, io.integrationTolerance);

	// On Process 0: output data will be stored here:
	bool anyFailures = false;
	bool anySuccesses = false;

	// On process 0: create a buffer for receiving results from other processes.
	// Both TE and TM results share the same array layout (eff[2N+1] + 4 scalar fields),
	// so we pack one TE result followed by one TM result into each rank's send buffer.
	const std::size_t resultSize = static_cast<std::size_t>(2*io.N + 1 + 4);
	const std::size_t pairSize = 2 * resultSize; // [TE fields][TM fields]

	std::vector<double> mpiReceiveBuffer;
	if(rank == 0)
		mpiReceiveBuffer.resize(pairSize * commandSize);

	std::vector<double> mpiSendBuffer(pairSize);

	// On Process 0: separate result histories per polarization / combination.
	std::vector<TEResult> resultsTE;
	std::vector<TMResult> resultsTM;
	std::vector<TEResult> resultsCombined; // reuses TEResult's layout for the averaged spectrum

	for(int i=0; i<totalSteps; i+=commandSize) {
		...
		TEResult resultTE = TEResult(TEResult::InactiveCalculation);
		TMResult resultTM = TMResult(TMResult::InactiveCalculation);

		if(i+rank < totalSteps){
			if(io.computeTE || io.combineTETM)
				resultTE = grating->getEffTE(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, (io.printDebugOutput && rank == 0), io.threads);
			if(io.computeTM || io.combineTETM)
				resultTM = grating->getEffTM(incidenceAngle, wavelength, io.rmsRoughnessNm, mathOptions, (io.printDebugOutput && rank == 0), io.threads);
		}

		// Pack both results into this rank's send buffer.
		resultTE.toDoubleArray(mpiSendBuffer.data());
		resultTM.toDoubleArray(mpiSendBuffer.data() + resultSize);

		int err = MPI_Gather(mpiSendBuffer.data(), static_cast<int>(pairSize), MPI_DOUBLE,
                      mpiReceiveBuffer.data(), static_cast<int>(pairSize), MPI_DOUBLE,
                      0, MPI_COMM_WORLD);

		if(err != MPI_SUCCESS) {
			char errStr[MPI_MAX_ERROR_STRING];
			int errLen = 0;
			MPI_Error_string(err, errStr, &errLen);
			std::cerr << "MPI_Gather failed on rank " << rank << ": " << errStr << std::endl;
			MPI_Abort(MPI_COMM_WORLD, err);
		}

		// On Process 0: collect results and output.
		if(rank == 0) {
			for(int j=0; j<commandSize; ++j) {
				TEResult resultTEj;
				TMResult resultTMj;
				resultTEj.fromDoubleArray(mpiReceiveBuffer.data() + j*pairSize);
				resultTMj.fromDoubleArray(mpiReceiveBuffer.data() + j*pairSize + resultSize);

				// indicates non-calculation for inactive process on last round
				if(resultTEj.status != TEResult::InactiveCalculation) {
					if(io.computeTE) {
						resultsTE.push_back(resultTEj);
						if(resultTEj.status == TEResult::Success) anySuccesses = true; else anyFailures = true;
					}
					if(io.computeTM) {
						resultsTM.push_back(resultTMj);
						if(resultTMj.status == TMResult::Success) anySuccesses = true; else anyFailures = true;
					}
					if(io.combineTETM && resultTEj.status == TEResult::Success && resultTMj.status == TMResult::Success) {
						TEResult combined = resultTEj;
						for(std::size_t k = 0; k < combined.eff.size(); ++k)
							combined.eff[k] = (resultTEj.eff[k] + resultTMj.eff[k]) / 2.0;
						resultsCombined.push_back(combined);
					}
				}
			}

			// Print progress and results to output file.
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
				writeOutputFileProgress(progressFile, std::min(i+commandSize, totalSteps), totalSteps, anySuccesses, anyFailures);
			}
		}

	} // end of calculation loop.

	// Timing: We know we're synchronized here because the last MPI_Gather has ensured that we have everyone's results.
	double runTime = MPI_Wtime() - startTime;
	if(rank == 0)
		std::cout << "Run time (s): " << runTime << std::endl;

	outputFile.close();

	// Finalize MPI
	MPI_Finalize();
	return 0;
}