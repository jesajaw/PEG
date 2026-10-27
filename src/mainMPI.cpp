/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/

#include "PEG.h"
#include "mainSupport.h"
#include "mpi.h"

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <memory>

/*
This main program provides a command-line interface to run a series of grating efficiency calculations in parallel with MPI. The results are written to an output file and (optionally) a second progress file. This file is only responsible for input processing, distributing the sweep points and output; everything that is identical to the sequential version (options, sweep point, grating, output format) lives in mainSupport, and all numerical details are structured within Grating and both Solvers.

Distribution: in every round each of the 'commandSize' processes calculates one sweep point (step i+rank); process 0 gathers all results of the round, appends them in step order and rewrites the output and progress file. All processes take part in every MPI_Gather, so none of them may leave the loop or throw out of it (computeSweepPoint never throws).
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
	
	// is the chosen material available?
	std::string dataError;
	if(!checkMaterialData(io, dataError)) {
		if(rank == 0) std::cerr << dataError << std::endl;
		MPI_Finalize();
		return -1;
	}
	applyThreadSettings(io);
	
	// copyright notice: the full legal text with --showLegal, otherwise the short one; once, from process 0
	if(rank == 0)
		printLegalNotice(std::cout, io.showLegal);

	// How many steps do we have? number of points of the sweep 'min'..'max' in steps of 'increment'
	int totalSteps = numSweepSteps(io);

	// On Process 0: Open the output file, check the progress file, write the header
	std::ofstream outputFile;
	std::streampos outputFilePosition;
	if(rank == 0) {
		outputFile.open(io.outputFile.c_str(), std::ios::out | std::ios::trunc);

		if(!outputFile.is_open()) {
			std::cerr << "Could not open output file " << io.outputFile << std::endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// Write the initial progress:
		///\todo also visible to extern exc?
		if(!io.progressFile.empty()) {
			std::ofstream progressFile(io.progressFile.c_str(), std::ios::out | std::ios::trunc);
			if(!progressFile.is_open()) {
				std::cerr << "Could not open progress file " << io.progressFile << std::endl;
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			writeOutputFileProgress(progressFile, 0, totalSteps, false, false);
		}

		// Write the file header: '# Input' block
		writeOutputFileHeader(outputFile, io);

		// Remember this position in the output file; it is where we will (repeatedly) write the result lines
		outputFilePosition = outputFile.tellp();
	}
	
	
	// create the grating object.
	std::unique_ptr<Grating> grating = makeGrating(io);
	if(!grating) {
		// this should never happen; input validation assures one of the valid grating types (and all processes get here together).
		if(rank == 0) std::cerr << "Unknown Grating-Profil." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
		return -1;
	}

	// math options: numerical settings shared by all points: truncation index 'N' and ODE tolerance
	MathOptions mathOptions(io.N, io.integrationTolerance);

	// On Process 0: output data will be stored here:
	bool anyFailures = false;
	bool anySuccesses = false;

	// ... separate result histories per polarization / combination.
	std::vector<Result> resultsTE;
	std::vector<Result> resultsTM;

	// ... create a buffer for receiving results from other processes.
	// Both TE and TM results share the same array layout (eff[2N+1] + 4 scalar fields), so we pack one TE result followed by one TM result into each rank's send buffer.
	const std::size_t resultSize = static_cast<std::size_t>(2*io.N + 1 + 4);
	const std::size_t pairSize = 2 * resultSize; // [TE fields][TM fields]

	std::vector<double> mpiReceiveBuffer;
	if(rank == 0)
		mpiReceiveBuffer.resize(pairSize * commandSize);
	std::vector<double> mpiSendBuffer(pairSize);


	// calculation loop: one round = 'commandSize' sweep points, one per process
	for(int i=0; i<totalSteps; i+=commandSize) {

		// results start as 'InactiveCalculation' (= nothing to do / not requested); the last round may have fewer points than processes.
		PointResults point;
		if(i+rank < totalSteps)
			point = computeSweepPoint(*grating, io, mathOptions, i+rank, io.printDebugOutput && rank == 0);

		// Pack both results into this rank's send buffer.
		point.te.toDoubleArray(mpiSendBuffer.data());
		point.tm.toDoubleArray(mpiSendBuffer.data() + resultSize);

		int err = MPI_Gather(mpiSendBuffer.data(), static_cast<int>(pairSize), MPI_DOUBLE, mpiReceiveBuffer.data(), static_cast<int>(pairSize), MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
				Result resultTEj;
				Result resultTMj;
				resultTEj.fromDoubleArray(mpiReceiveBuffer.data() + j*pairSize);
				resultTMj.fromDoubleArray(mpiReceiveBuffer.data() + j*pairSize + resultSize);

				// indicates non-calculation for inactive process on last round -- whichever result was actually requested tells us whether this rank had work this round
				bool wasActive = io.TE ? (resultTEj.status != Result::InactiveCalculation)
				                       : (resultTMj.status != Result::InactiveCalculation);
				if(wasActive) {
					// failed points are kept too (written as an 'Error:...' row), so every section has one row per point.
					if(io.TE) {
						resultsTE.push_back(resultTEj);
						if(resultTEj.status == Result::Success) anySuccesses = true; else anyFailures = true;
					}
					if(io.TM) {
						resultsTM.push_back(resultTMj);
						if(resultTMj.status == Result::Success) anySuccesses = true; else anyFailures = true;
					}
				}
			}

			// Print results to output file
			writeOutputFileResults(outputFile, outputFilePosition, resultsTE, resultsTM, io);

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