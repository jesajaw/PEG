/*
Copyright (C) 2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de) and 2012 Mark Boots (mark.boots@usask.ca).

This program was originally implemented as a part of the Parallel Efficiency of Gratings project PEG and got reworked in 2026. PEG is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation.
See <http://www.gnu.org/licenses/> for details.

This reworked version contains substantial modifications by Jesaja Weintritt (2026) and has not been independently verified against the original. It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; use at your own risk and verify results independently.
*/


#include "mainSupport.h"
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <exception>
#include <sstream>
#include <getopt.h>
#include <Eigen/Core>

/// Converts the whole string \c s to a finite double, otherwise throws \c errMsg.
static double parseDouble(const char* s, const char* errMsg) {
	char* end = nullptr;
	errno = 0;
	double v = strtod(s, &end);
	if(end == s || *end != '\0' || errno == ERANGE || !std::isfinite(v)) throw errMsg;
	return v;
}

/// Converts the whole string \c s to an integer or throw errMsg
static int parseInt(const char* s, const char* errMsg) {
	char* end = nullptr;
	errno = 0;
	long v = strtol(s, &end, 10);
	if(end == s || *end != '\0' || errno == ERANGE || v < INT_MIN || v > INT_MAX) throw errMsg;
	return static_cast<int>(v);
}

void CommandLineOptions::init() {
	mode = InvalidMode;

	/// "not given" markers; isValid() reports every required option that still has one of them.
	min = max = increment = incidenceAngle = includedAngle = wavelength = DBL_MAX;

	toOrder = N = INT_MAX;

	profile = Grating::InvalidProfile;
	period = DBL_MAX;

	eV = false;
	printDebugOutput = false;
	threads = 1;	/// by default, just one thread.
	measureTiming = false;
	integrationTolerance = 1e-5;	/// default: 1e-5 if not provided.
	coatingThickness = 0;	/// default: 0 (no coating) if not provided.

	rmsRoughnessNm = 0;

	showLegal = false;

	/// Polarization calculation options
	TE = false;
	TM = false;
}

/// Sets options based on command-line input arguments. Returns isValid().
bool CommandLineOptions::parseFromCommandLine(int argc, char** argv) {
	/// Parse Command-line options
	try {
		/// loop over all input options
		while(1) {
			/// define the options available
			static struct option long_options[] = {
				{"mode", required_argument, 0, 1},
				{"min", required_argument, 0, 2},
				{"max", required_argument, 0, 3},
				{"increment", required_argument, 0, 4},
				{"incidenceAngle", required_argument, 0, 5},
				{"includedAngle", required_argument, 0, 6},
				{"toOrder", required_argument, 0, 7},
				{"wavelength", required_argument, 0, 8},
				{"outputFile", required_argument, 0, 9},
				{"progressFile", required_argument, 0, 10},
				{"gratingType", required_argument, 0, 11},
				{"gratingGeometry", required_argument, 0, 12},
				{"gratingPeriod", required_argument, 0, 13},
				{"gratingMaterial", required_argument, 0, 14},
				{"N", required_argument, 0, 15},
				{"eV", no_argument, 0, 16},
				{"printDebugOutput", no_argument, 0, 17},
				{"threads", required_argument, 0, 18},
				{"measureTiming", no_argument, 0, 19},
				{"integrationTolerance", required_argument, 0, 20},
				{"coatingMaterial", required_argument, 0, 21},
				{"coatingThickness", required_argument, 0, 22},
				{"showLegal", no_argument, 0, 23},
				{"rmsRoughnessNm", required_argument, 0, 24},
				{"TE", no_argument, 0, 25},
				{"TM", no_argument, 0, 26},
				{"TETM", no_argument, 0, 27},
				{"TMTE", no_argument, 0, 27},
				{0, 0, 0, 0}
			};
				
			int option_index = 0;
			int c = getopt_long (argc, argv, "", long_options, &option_index);

			/// Detect the end of the options.
			if (c == -1)
				break;

			switch (c) {
				case 1: /// mode
					if(strcmp(optarg, "constantIncidence") == 0) mode = ConstantIncidence;
					else if(strcmp(optarg, "constantIncludedAngle") == 0) mode = ConstantIncludedAngle;
					else if(strcmp(optarg, "constantWavelength") == 0) mode = ConstantWavelength;
					else throw "The argument to --mode must be one of: constantIncidence, constantIncludedAngle, or constantWavelength.";
					break;
				case 2: /// min
					min = parseDouble(optarg, "--min must be a number.");
					break;
				case 3: /// max
					max = parseDouble(optarg, "--max must be a number.");
					break;
				case 4: /// increment
					increment = parseDouble(optarg, "--increment must be a number.");
					break;
				case 5: /// incidence angle
					incidenceAngle = parseDouble(optarg, "--incidenceAngle must be a number.");
					break;
				case 6: /// includedAngle
					includedAngle = parseDouble(optarg, "--includedAngle must be a number.");
					break;
				case 7: /// toOrder
					toOrder = parseInt(optarg, "--toOrder must be an integer.");
					break;
				case 8: /// wavelength
					wavelength = parseDouble(optarg,"--wavelength must be a number.");
					break;
				case 9: /// outputFile
					outputFile = optarg;
					break;
				case 10: /// progressFile
					progressFile = optarg;
					break;
				case 11: /// gratingType
					if(strcmp(optarg, "rectangular") == 0) profile = Grating::RectangularProfile;
					else if(strcmp(optarg, "blazed") == 0) profile = Grating::BlazedProfile;
					else if(strcmp(optarg, "sinusoidal") == 0) profile = Grating::SinusoidalProfile;
					else if(strcmp(optarg, "trapezoidal") == 0) profile = Grating::TrapezoidalProfile;
					else if(strcmp(optarg, "custom") == 0) profile = Grating::CustomProfile;
					else throw "The argument to --gratingType must be one of: rectangular, blazed, sinusoidal, trapezoidal, or custom.";
					break;
				case 12: {/// gratingGeometry
					char* saveptr;
					char* token;
					token = strtok_r(optarg, ",", &saveptr);
					while(token) {
						geometry.push_back(parseDouble(token, "--gratingGeometry values must be numbers."));
						token = strtok_r(0, ",", &saveptr);
					}
					break;
					}
				case 13: /// gratingPeriod
					period = parseDouble(optarg, "--gratingPeriod must be a number.");
					break;
				case 14: /// gratingMaterial
					material = optarg;
					break;
				case 15: /// N
					N = parseInt(optarg, "--N must be an integer.");
					break;
				case 16: /// interpret as eV instead of wavelength
					eV = true;
					break;
				case 17: /// print debug output to stdout
					printDebugOutput = true;
					break;
				case 18: /// number of threads to use for fine parallelization
					threads = parseInt(optarg, "--threads must be an integer.");
					break;
				case 19: /// measure timing
					measureTiming = true;
					break;
				case 20: /// integrationTolerance
					integrationTolerance = parseDouble(optarg, "--integrationTolerance must be a number.");
					break;
				case 21: /// coatingMaterial
					coating = optarg;
					break;
				case 22: /// coatingThickness
					coatingThickness = parseDouble(optarg, "--coatingThickness must be a number.");
					break;
				case 23: /// showLegal
					showLegal = true;
					break;
				case 24: /// rmsRoughnessNm
					rmsRoughnessNm = parseDouble(optarg, "--rmsRoughnessNm must be a number.");
					break;
				case 25: /// compute TE
					TE = true;
					break;
				case 26: /// compute TM
					TM = true;
					break;
				case 27: /// TE & TM
					TE = true;
					TM = true;
					break;
				default:	/// getopt_long returned '?': unknown option or missing argument
					throw "Unknown option or missing option argument.";
			}
		} /// end of loop over input options.
				
	}
	
	catch(const char* errMsg) {
		firstErrorMessage_ = errMsg;
		return false;
	}
	
	return isValid();
}

/// Returns true if all required input options have been provided correctly.
bool CommandLineOptions::isValid() {
	try {
		// (1) everything required is given
		if(mode == InvalidMode) throw "The operating mode --mode must be provided: constantIncidence, constantIncludedAngle, or constantWavelength";
		if(min == DBL_MAX) throw "The minimum range --min must be provided.";
		if(max == DBL_MAX) throw "The maximum range --max must be provided.";
		if(increment == DBL_MAX) throw "The range increment --increment must be provided.";

		if(mode == ConstantIncidence && incidenceAngle == DBL_MAX) throw "In constant incidence mode, the incidence angle --incidenceAngle must be provided.";
		if(mode == ConstantIncludedAngle && includedAngle == DBL_MAX) throw "In constant included angle mode, the included angle --includedAngle must be provided.";
		if(mode == ConstantIncludedAngle && toOrder == INT_MAX) throw "In constant included angle mode, the operating order --toOrder must be provided.";
		if(mode == ConstantWavelength && wavelength == DBL_MAX) throw "In constant wavelength mode, the wavelength --wavelength must be provided.";

		if(outputFile.empty()) throw "The output data file --outputFile must be provided.";

		if(profile == Grating::InvalidProfile) throw "The grating profile --gratingType must be provided.";
		if(material.empty()) throw "The grating material --gratingMaterial must be provided.";
		if(period == DBL_MAX) throw "The grating period --gratingPeriod must be provided.";

		if(N == INT_MAX) throw "The truncation index --N must be provided.";
		if(!TE && !TM) throw "At least one of --TE or --TM must be specified.";

		// (2) all numbers are finite
		if(!std::isfinite(min) || !std::isfinite(max) || !std::isfinite(increment))
			throw "--min, --max, and --increment must all be finite numbers.";
		if(incidenceAngle != DBL_MAX && !std::isfinite(incidenceAngle)) throw "--incidenceAngle must be a finite number.";
		if(includedAngle != DBL_MAX && !std::isfinite(includedAngle)) throw "--includedAngle must be a finite number.";
		if(wavelength != DBL_MAX && !std::isfinite(wavelength)) throw "--wavelength must be a finite number.";
		if(!std::isfinite(period)) throw "The grating period --gratingPeriod must be a finite number.";
		for(double gv : geometry) if(!std::isfinite(gv)) throw "Every --gratingGeometry argument must be a finite number.";

		// (3) value ranges
		// the sweep
		if(!(increment > 0)) throw "The range increment --increment must be a positive number.";
		if(max < min) throw "The range maximum --max must not be smaller than the minimum --min.";
		if(!(std::floor((max - min)/increment + 1e-9) < INT_MAX)) throw "The sweep has too many points; increase --increment or reduce the range --min..--max.";

		// the mode: angles in degrees, wavelengths/energies positive
		if(mode == ConstantIncidence && !(incidenceAngle >= 0 && incidenceAngle < 90))
			throw "The incidence angle --incidenceAngle must be in the range [0,90) degrees.";
		if(mode == ConstantIncludedAngle && !(includedAngle > 0 && includedAngle < 180))
			throw "The included angle --includedAngle must be in the range (0,180) degrees.";
		if(mode == ConstantWavelength && !(wavelength > 0)) throw "The wavelength --wavelength must be a positive number.";
		if(mode == ConstantWavelength && !(min >= 0 && max < 90))
			throw "In constant wavelength mode, --min and --max (the incidence angle sweep) must fall in [0,90) degrees.";
		if(mode != ConstantWavelength && !(min > 0))
			throw (eV ? "In this mode, --min and --max (the energy sweep) must be positive."
					: "In this mode, --min and --max (the wavelength sweep) must be positive.");

		// the grating, the numerics
		if(!(period > 0)) throw "The grating period --gratingPeriod must be a positive number.";
		if(coatingThickness < 0) throw "The coating thickness --coatingThickness must not be negative.";
		if(coatingThickness != 0 && coating.empty()) throw "A coating thickness was specified, but this requires a --coatingMaterial.";
		if(!coating.empty() && coatingThickness == 0) throw "A coating material was specified, but this requires a non-zero --coatingThickness.";
		if(rmsRoughnessNm < 0) throw "The RMS roughness must be in nm, larger than or equal to 0.";

		if(N < 0) throw "The truncation index --N must not be negative.";
		if(threads < 1) throw "The number of --threads to use for fine parallelization must be a positive number, at least 1.";
		if(!(integrationTolerance > 0)) throw "The integration tolerance --integrationTolerance must be a positive number.";

		// (4) the profile: number of --gratingGeometry values, then their ranges
		if(profile == Grating::RectangularProfile && geometry.size() != 2) throw "The rectangular profile requires two arguments to --gratingGeometry <depth>,<valleyWidth>.";
		if(profile == Grating::BlazedProfile && geometry.size() != 2) throw "The blazed profile requires two arguments to --gratingGeometry <blazeAngleDeg>,<antiBlazeAngleDeg>.";
		if(profile == Grating::SinusoidalProfile && geometry.size() != 1) throw "The sinusoidal profile requires one argument to --gratingGeometry <depth>.";
		if(profile == Grating::TrapezoidalProfile && geometry.size() != 4) throw "The trapezoidal profile requires four arguments to --gratingGeometry <depth>,<valleyWidth>,<blazeAngle>,<antiBlazeAngle>.";
		if(profile == Grating::CustomProfile && (geometry.size() < 7 || geometry.size()%2 != 1)) throw "The custom (point-wise) profile requires arguments to --gratingGeometry: a maximum height (um), followed by a sequence of (x,y) points along the profile going from (0,0) to (1,0). They will be scaled so that (0,0)->(1,1) maps to (0,0)->(period,maxHeight).";

		if(profile == Grating::CustomProfile && coatingThickness != 0) throw "The custom profile does not yet support coatings.";

		if(profile == Grating::RectangularProfile) {
			if(!(geometry[0] > 0)) throw "The rectangular profile depth (first --gratingGeometry value) must be positive.";
			if(!(geometry[1] >= 0 && geometry[1] < period)) throw "The rectangular profile valleyWidth (second --gratingGeometry value) must be in [0,period).";
		}
		if(profile == Grating::BlazedProfile) {
			if(!(geometry[0] > 0 && geometry[0] < 90)) throw "The blazed profile blazeAngle (first --gratingGeometry value) must be in (0,90) degrees.";
			if(!(geometry[1] > 0 && geometry[1] < 90)) throw "The blazed profile antiBlazeAngle (second --gratingGeometry value) must be in (0,90) degrees.";
		}
		if(profile == Grating::SinusoidalProfile) {
			if(!(geometry[0] > 0)) throw "The sinusoidal profile depth (the --gratingGeometry value) must be positive.";
		}
		if(profile == Grating::TrapezoidalProfile) {
			if(!(geometry[0] > 0)) throw "The trapezoidal profile depth (first --gratingGeometry value) must be positive.";
			if(!(geometry[1] >= 0 && geometry[1] < period)) throw "The trapezoidal profile valleyWidth (second --gratingGeometry value) must be in [0,period).";
			if(!(geometry[2] > 0 && geometry[2] < 90)) throw "The trapezoidal profile blazeAngle (third --gratingGeometry value) must be in (0,90) degrees.";
			if(!(geometry[3] > 0 && geometry[3] < 90)) throw "The trapezoidal profile antiBlazeAngle (fourth --gratingGeometry value) must be in (0,90) degrees.";
		}
		if(profile == Grating::CustomProfile) {
			if(!(geometry[0] > 0)) throw "The custom profile maximum height (first --gratingGeometry value) must be positive.";
			// The solver slices the profile from the first point on, and reads the point before the first one that lies above the slice: a profile that does not start at (0,0) would read outside of its point list.
			size_t last = geometry.size() - 2;	// index of the x-value of the last point
			if(geometry[1] != 0 || geometry[2] != 0 || geometry[last] != 1 || geometry[last+1] != 0)
				throw "The custom profile must start at (0,0) and end at (1,0).";
			for(size_t k=1; k<geometry.size(); ++k)
				if(!(geometry[k] >= 0 && geometry[k] <= 1)) throw "The custom profile points (x,y) must all lie in [0,1].";
		}
	}
	
	catch(const char* errMsg) {
		firstErrorMessage_ = errMsg;
		return false;
	}
	
	firstErrorMessage_ = "No errors.";
	return true;
}

// Sweep: from the options to one calculated point
int numSweepSteps(const CommandLineOptions& io) {
	return int(std::floor((io.max - io.min) / io.increment + 1e-9)) + 1;
}

bool checkMaterialData(const CommandLineOptions& io, std::string& error) {
	/// wavelengths (um) at both ends of the sweep
	double wlA, wlB;
	if(io.mode == CommandLineOptions::ConstantWavelength)
		wlA = wlB = io.eV ? M_HC / io.wavelength : io.wavelength;
	else {
		double first = io.min;
		double last  = io.min + io.increment * (numSweepSteps(io) - 1);
		wlA = io.eV ? M_HC / first : first;
		wlB = io.eV ? M_HC / last  : last;
	}

	// refractiveIndex() returns (0,0) for a missing, empty, unsorted table or a wavelength outside of the table.
	for(const std::string* mat : {&io.material, &io.coating}) {
		if(mat->empty()) continue;	// no coating
		for(double wl : {wlA, wlB})
			if(Grating::refractiveIndex(wl, *mat) == std::complex<double>(0,0)) {
				std::ostringstream msg;
				msg << "No refractive-index data for '" << *mat << "' at " << wl << " um (table missing, unsorted, or wavelength outside its range).";
				error = msg.str();
				return false;
			}
	}
	return true;
}

std::unique_ptr<Grating> makeGrating(const CommandLineOptions& io) {
	///\todo use array for coating & thickness -> multilayer
	switch(io.profile) {
		case Grating::RectangularProfile:
			return std::make_unique<RectangularGrating>(io.period, io.geometry[0], io.geometry[1], io.material, io.coating, io.coatingThickness);
		case Grating::BlazedProfile:
			return std::make_unique<BlazedGrating>(io.period, io.geometry[0], io.geometry[1], io.material, io.coating, io.coatingThickness);
		case Grating::SinusoidalProfile:
			return std::make_unique<SinusoidalGrating>(io.period, io.geometry[0], io.material, io.coating, io.coatingThickness);
		case Grating::TrapezoidalProfile:
			return std::make_unique<TrapezoidalGrating>(io.period, io.geometry[0], io.geometry[1], io.geometry[2], io.geometry[3], io.material, io.coating, io.coatingThickness);
		case Grating::CustomProfile:
			// no coating: isValid() refuses it for this profile.
			///\todo , io.coating, io.coatingThickness
			return std::make_unique<CustomProfileGrating>(io.period, io.geometry, io.material);
		default:
			return nullptr;
	}
}

void applyThreadSettings(const CommandLineOptions& io) {
	Eigen::setNbThreads(io.threads);
}

bool sweepPoint(const CommandLineOptions& io, int step, double& wavelength, double& incidenceDeg) {
	// value of the swept quantity: wavelength/energy or incidence angle in degrees
	double currentValue = io.min + io.increment*step;

	// wavelength (um): fixed (--wavelength) in constant-wavelength mode, swept in the other modes
	wavelength = (io.mode == CommandLineOptions::ConstantWavelength) ? io.wavelength : currentValue;
	// interpret the input as eV instead, and convert to the actual wavelength.
	// Formula: wavelength = hc / eV, with hc = 1.23984172 eV * um.
	if(io.eV)
		wavelength = M_HC / wavelength;

	// incidence angle: depends on the mode and possibly on the wavelength.
	switch(io.mode) {
		case CommandLineOptions::ConstantIncidence:		// fixed --incidenceAngle
			incidenceDeg = io.incidenceAngle;
			break;
		case CommandLineOptions::ConstantIncludedAngle: {	// monochromator mode
			double ciaRad = io.includedAngle * M_PI / 180;
			// satisfies alpha + beta = cia, and the grating equation toOrder*wavelength/d = sin(beta) - sin(alpha).
			// The asin() gives NaN if the order does not exist at this wavelength.
			incidenceDeg = (asin(-io.toOrder*wavelength/2/io.period/cos(ciaRad/2)) + ciaRad/2) * 180 / M_PI;
			break;
			}
		case CommandLineOptions::ConstantWavelength:		// the swept value itself is the incidence angle
			incidenceDeg = currentValue;
			break;
		default:
			// never happens: input validation assures a valid mode.
			incidenceDeg = 0;
			break;
	}

	return std::isfinite(incidenceDeg) && std::fabs(incidenceDeg) < 90;
}

PointResults computeSweepPoint(const Grating& grating, const CommandLineOptions& io, const MathOptions& mathOptions, int step, bool printDebugOutput) {
	PointResults r;
	double wavelength, incidenceDeg;
	bool usable = sweepPoint(io, step, wavelength, incidenceDeg);

	/// results start as 'InactiveCalculation' (= not requested); a requested polarization overwrites it.
	Result bad(Result::OtherFailure);
	bad.wavelength = wavelength;
	bad.incidenceDeg = incidenceDeg;

	if(!usable) {
		if(io.TE) r.te = bad;
		if(io.TM) r.tm = bad;
		return r;
	}

	// One complete solve per polarization: the efficiencies of all orders -N..+N. An exception must not end the sweep
	if(io.TE) {
		try {
			r.te = grating.getEffTE(incidenceDeg, wavelength, io.rmsRoughnessNm, mathOptions, printDebugOutput, io.threads, io.measureTiming);
		}
		catch(const std::exception& e) {
			std::cerr << "Step " << step << " (TE) failed with an exception: " << e.what() << std::endl;
			r.te = bad;
		}
	}
	if(io.TM) {
		try {
			r.tm = grating.getEffTM(incidenceDeg, wavelength, io.rmsRoughnessNm, mathOptions, printDebugOutput, io.threads, io.measureTiming);
		}
		catch(const std::exception& e) {
			std::cerr << "Step " << step << " (TM) failed with an exception: " << e.what() << std::endl;
			r.tm = bad;
		}
	}
	return r;
}

void printLegalNotice(std::ostream& os, bool full) {
	if(full) {
		os <<
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
		os <<
		"PEG Copyright (C)\n"
		"2026 Jesaja Weintritt (jesaja.weintritt@stud.eah-jena.de)\n"
		"2012 Mark Boots (mark.boots@usask.ca)\n"
		"This free software comes with ABSOLUTELY NO WARRANTY; you are welcome to redistribute it under certain conditions; run with --showLegal for details.\n";
	}
}

//// This helper function writes the header to the output file stream
void writeOutputFileHeader(std::ostream& of, const CommandLineOptions& io) {
	of << "# Input" << std::endl;
	
	switch(io.mode) {
		case CommandLineOptions::ConstantIncidence:
		of << "mode=constantIncidence" << std::endl;
		of << "incidenceAngle=" << io.incidenceAngle << std::endl;
		break;
		case CommandLineOptions::ConstantIncludedAngle:
		of << "mode=constantIncludedAngle" << std::endl;
		of << "includedAngle=" << io.includedAngle << std::endl;
		of << "toOrder=" << io.toOrder << std::endl;
		break;
		case CommandLineOptions::ConstantWavelength:
		of << "mode=constantWavelength" << std::endl;
		of << "wavelength/eV=" << io.wavelength << std::endl;
		break;
		default:
		of << "mode=invalid" << std::endl;
		break;
	}
	
	if(io.eV)
		of << "units=eV" << std::endl;
	else
		of << "units=um" << std::endl;
		
	of << "min=" << io.min << std::endl;
	of << "max=" << io.max << std::endl;
	of << "increment=" << io.increment << std::endl;
	
	switch(io.profile) {
		case Grating::RectangularProfile:
		of << "gratingType=rectangular" << std::endl;
		break;
		case Grating::BlazedProfile:
		of << "gratingType=blazed" << std::endl;
		break;
		case Grating::SinusoidalProfile:
		of << "gratingType=sinusoidal" << std::endl;
		break;
		case Grating::TrapezoidalProfile:
		of << "gratingType=trapezoidal" << std::endl;
		break;
		case Grating::CustomProfile:
		of << "gratingType=custom" << std::endl;
		break;
		default:
		of << "gratingType=invalid" << std::endl;
		break;
	}
	
	of << "gratingPeriod=" << io.period << std::endl;
	
	of << "gratingGeometry=";
	for(int i=0,cc=io.geometry.size(); i<cc; ++i) {
		if(i!=0)
			of << ",";
		of << io.geometry.at(i);
	}
	of << std::endl;
	
	of << "gratingMaterial=" << io.material << std::endl;
	if(io.coatingThickness > 0) {
		of << "coatingMaterial=" << io.coating << std::endl;
		of << "coatingThickness=" << io.coatingThickness << std::endl;
	}
	else {
		of << "coatingMaterial=[none]" << std::endl;
		of << "coatingThickness=0" << std::endl;
	}
	
	of << "N=" << io.N << std::endl;
	of << "integrationTolerance=" << io.integrationTolerance << std::endl;
}

/// This helper function appends a single efficiency result to the output file stream
void writeOutputFileResult(std::ostream& of, const Result& result, const CommandLineOptions& io) {
	/// x-column depends on mode: incidence angle in ConstantWavelength mode,
	/// otherwise wavelength (converted to eV if io.eV is set).
	if(io.mode == CommandLineOptions::ConstantWavelength) {
		of << result.incidenceDeg << "\t";
	}
	else {
		of << (io.eV ? M_HC / result.wavelength : result.wavelength) << "\t";
	}
	
	switch(result.status) {
	case Result::InvalidGratingFailure:
		of << "Error:InvalidGratingFailure" << std::endl;
		break;
	case Result::ConvergenceFailure:
		of << "Error:ConvergenceFailure" << std::endl;
		break;
	case Result::InsufficientCoefficientsFailure:
		of << "Error:InsufficientCoefficientsFailure" << std::endl;
		break;
	case Result::AlgebraFailure:
		of << "Error:AlgebraError" << std::endl;
		break;
	case Result::OtherFailure:
		of << "Error:OtherFailure" << std::endl;
		break;
	case Result::InactiveCalculation:
		of << "Error:InactiveCalculation" << std::endl;
		break;
	case Result::MissingRefractiveDataFailure:
		of << "Error:MissingRefractiveIndexData" << std::endl;
		break;
	case Result::Success:
		for(int i=0, cc=result.eff.size(); i<cc; ++i) {
			if(i!=0)
				of << "\t";
			of << result.eff.at(i);
		}
		of << std::endl;
		break;
	}
}

void writeOutputFileResults(std::ostream& of, std::streampos position, const std::vector<Result>& resultsTE, const std::vector<Result>& resultsTM, const CommandLineOptions& io) {
	of.seekp(position);
	of << "# Output" << std::endl;
	if(io.TE) {
		of << "# TE" << std::endl;
		for(const auto& result : resultsTE)
			writeOutputFileResult(of, result, io);
	}
	if(io.TM) {
		of << "# TM" << std::endl;
		for(const auto& result : resultsTM)
			writeOutputFileResult(of, result, io);
	}
	of.flush();
}

/// This helper function appends the progress to the given output stream
void writeOutputFileProgress(std::ostream& of, int completedSteps, int totalSteps, bool anySuccesses, bool anyFailures) {
	of << "# Progress" << std::endl;
	/// not done yet:
	if(completedSteps < totalSteps) {
		of << "status=inProgress" << std::endl;
	}
	/// or all done:
	else {
		if(anySuccesses && anyFailures)
			of << "status=someFailed" << std::endl;
		else if(anyFailures)
			of << "status=allFailed" << std::endl;
		else
			of << "status=succeeded" << std::endl;
	}
	
	of << "completedSteps=" << completedSteps << std::endl;
	of << "totalSteps=" << totalSteps << std::endl;
}