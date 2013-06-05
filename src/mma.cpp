/***************************************************************************
                                   mma.cpp
                             -------------------
    update               : 2005/06/30
    copyright            : (C) 2002-2005 by Michaël Roy
    email                : michaelroy@users.sourceforge.net
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#define _MMA_GUI_


//
// This the main program
// There is the choice between the graphical and the command line interface
//



//--
//
// Graphical Interface
//
//--
#ifdef _MMA_GUI_


// Grahical user interface class
#include "GraphicalInterface.h"

// Main program
int main(int argc, char **argv)
{
	// Create the main GUI
	GraphicalInterface* gui = new GraphicalInterface;
	// Show the GUI
	gui->Show(argc,argv);
	// FLTK main loop
	return Fl::run();
}


//--
//
// Command line interface
//
//--
#else


// Local references
#include "MultiresolutionMesh.h"
#include "MultiresolutionProcessing.h"
#include "Stopwatch.h"

// Global references
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

// Standard C++ stuff
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::setw;

// Main banner
static const std::string banner_string =
"                mmac 1.0\n"
" Mesh Multiresolution Analysis Software\n"
"     Copyright (C) 2002-2005 Michael Roy\n";

// Usage banner
static const std::string usage_string =
"Usage: wavemesh [options] <file> <level>\n"
"   <file>  file to be analysed\n"
"   <level> number of resolution levels\n"
"\n"
"   Options:\n"
"\n"
"     -filter           Enable detail coefficient filtering\n"
"\n"
"     -hard             Hard threshold filtering\n"
"\n"
"     -soft             Soft threshold filtering [default]\n"
"\n"
"     -output <file>    Output result file\n"
"\n"
"     -colorize         Write all levels of detail with color-coded details\n"
"\n"
"     -threshold <num>  Filter threshold value [default value = maximum]\n"
"\n"
"     -statistics       Print detail coefficient statistics\n"
"\n"
"     -levelstatistics  Print detail coefficient statistics per level\n"
"\n"
"     -simpli <t> <c>   Simplify the model down with <t> threshold and <c> colorize\n"
"\n"
"     -segment <ll> <hl> <mdm>  Colorize mesh\n"
"                                 <ll>  Low level\n"
"                                 <hl>  High level\n"
"                                 <mdm> Minimum detail magnitude\n";

// Parameter variables
static std::string file_name;
static bool filter_enabled(false);
static double threshold(DBL_MAX);
static bool soft_threshold_enabled(true);
static bool colorize_enabled(false);
static bool segment_enabled(false);
static bool statistics_enabled(false);
static bool level_statistics_enabled(false);
static std::string output_filename("");
static int level_number(50);
static double simplification_threshold(0);
static double simplification_colorize(false);


//--
//
// Parameters
//
//--
// Parse the command line parameters
static bool Parameters( int argc, char *argv[] )
{
	// No parameters
	// Print out the usage banner
	if( argc < 2 ) {
		cout<<endl<<usage_string<<endl;
	    	exit( EXIT_FAILURE );
	}

	// Iterate through the parameter list
	for( int loop=1; loop<argc; loop++ ) {
		// Options
		if( !strcmp(argv[loop], "-h") || !strcmp(argv[loop], "-help") || !strcmp(argv[loop], "--help") ) {
			cout<<endl<<usage_string<<endl;
	    		exit( EXIT_SUCCESS );
		}
		// Filter
		else if( !strcmp(argv[loop], "-filter") ) {
			filter_enabled = true;
		}
		// Hard thresolding
		else if( !strcmp(argv[loop], "-hard") ) {
			soft_threshold_enabled = false;
		}
		// Soft thresholding
		else if( !strcmp(argv[loop], "-soft") ) {
			soft_threshold_enabled = true;
		}
		// Coloring
		else if( !strcmp(argv[loop], "-colorize") ) {
			colorize_enabled = true;
		}
		// Threshold
		else if( !strcmp(argv[loop], "-threshold") ) {
			if( loop > (argc-3) ) return false;
			threshold = atof( argv[++loop] );
		}
		// Output mesh
		else if( !strcmp(argv[loop], "-output") ) {
			if( loop > (argc-3) ) return false;
			output_filename = argv[++loop];
		}
		// Statistics
		else if( !strcmp(argv[loop], "-statistics") ) {
			statistics_enabled = true;
		}
		// Level Statistics
		else if( !strcmp(argv[loop], "-levelstatistics") ) {
			level_statistics_enabled = true;
		}
		// Simplification
		else if( !strcmp(argv[loop], "-simpli") ) {
			if( loop > (argc-4) ) return false;
			simplification_threshold = atof( argv[++loop] );
			simplification_colorize = (bool)atoi( argv[++loop] );
		}
		// Input file names
		else if( loop == (argc - 2) ) {
			file_name = argv[loop++];
			level_number = atoi( argv[loop] );
		}
		// Error
		else {
			return false;
		}
	}

	// Error tests
	if( file_name.empty() ) return false;

	// Let's go !
	return true;
}

//--
//
// main
//
//--
// Main program
// Command line interface
int main(int argc, char *argv[])
{
	// Get a timer
	Stopwatch timer;

	// Print out the software banner
	cout<<endl<<banner_string<<endl;

	// Analyse the command line parameters
	if( Parameters(argc, argv) == false ) {
		cout<<"Error in command line"<<endl;
		cout<<"use -h option to see manual"<<endl<<endl;
		return EXIT_FAILURE;
	}

	// The mesh
	MultiresolutionMesh mesh;

	// The processor
	MultiresolutionProcessing processor;

	//
	// Read file
	//
	cout<<"+ Read file...          "<<flush;
	timer.Start();
	if( !mesh.ReadFile(file_name) ) {
		cout<<"    Failed\n"<<endl;
		return EXIT_FAILURE;
	}
	timer.Stop();
	cout<<setw(6)<<timer.Intermediate()<<endl;

	//
	// Multiresolution decomposition
	//
	cout<<"+ Analysis...           "<<flush;
	timer.Start();
	if( !mesh.Analysis(level_number) ) {
		cout<<"    Failed\n"<<endl;
		return EXIT_FAILURE;
	}
	timer.Stop();
	cout<<setw(6)<<timer.Intermediate()<<endl;

	//
	// Filtering
	//
	if( filter_enabled ) {
		cout<<"+ Filtering...          "<<flush;
		timer.Start();
//		mesh.Filter(low_threshold, high_threshold, soft_threshold_enabled, low_level_bound, high_level_bound);
	//	DetailShrinkage(mesh, threshold, soft_threshold_enabled, low_level_bound, high_level_bound);
	//	AnisotropicFiltering(mesh);
		timer.Stop();
		cout<<setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Multiresolution recontruction
	//
	cout<<"+ Synthesis...          "<<flush;
	timer.Start();
	if( !mesh.Synthesis() ) {
		cout<<"    Failed\n"<<endl;
		return EXIT_FAILURE;
	}
	timer.Stop();
	cout<<setw(6)<<timer.Intermediate()<<endl;

	//
	// Colorization
	//
	if( colorize_enabled ) {
		cout<<"+ Colorize mesh...      "<<std::flush;
		timer.Start();
//		if( Colorize(mesh) == false ) {
//			cout<<"    Failed\n"<<endl;
//			return EXIT_FAILURE;
//		}
		timer.Stop();
		cout<<std::setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Thresholding
	//
	if( segment_enabled ) {
		cout<<"+ Segmentation...       "<<flush;
		timer.Start();
//		DetailSegmentation(mesh, detail_magnitude_color, low_level_color, high_level_color);
		timer.Stop();
		cout<<setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Simplification
	//
	if( simplification_threshold > 0.0 ) {
		// Simplication
		cout<<"+ Simplification...     "<<std::flush;
		timer.Start();
//		if( !Simplification(mesh, simplification_threshold, simplification_colorize) ) {
//			cout<<"    Failed\n"<<endl;
//			return EXIT_FAILURE;
//		}
		timer.Stop();
		cout<<std::setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Statistics
	//
	if( statistics_enabled ) {
		cout<<"+ Statistics...         "<<flush;
		timer.Start();
//		if( !Statistics(mesh) ) {
//			cout<<"    Failed\n"<<endl;
//			return EXIT_FAILURE;
//		}
		timer.Stop();
		cout<<setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Level Statistics
	//
	if( level_statistics_enabled ) {
		cout<<"+ Level statistics...   "<<flush;
		timer.Start();
//		if( !LevelStatistics(mesh) ) {
//			cout<<"    Failed\n"<<endl;
//			return EXIT_FAILURE;
//		}
		timer.Stop();
		cout<<setw(6)<<timer.Intermediate()<<endl;
	}

	//
	// Output file
	//
	if( !output_filename.empty() ) {
		cout<<"+ Compute normals...    "<<std::flush;
		timer.Start();
		mesh.ComputeFaceNormals();
		mesh.ComputeVertexNormals();
		timer.Stop();
		cout<<std::setw(6)<<timer.Intermediate()<<endl;
		cout<<"+ Write output file...  "<<std::flush;
		timer.Start();
		if( !mesh.WriteFile(output_filename, VRML_1_FILE) ) {
			cout<<"    Failed\n"<<endl;
			return EXIT_FAILURE;
		}
		timer.Stop();
		cout<<std::setw(6)<<timer.Intermediate()<<endl;
	}

	// Print out the total time
	cout<<endl<<"Mesh Multiresolution Analysis done in "<<timer.Total()<<endl<<endl;

	// The end !
	return EXIT_SUCCESS;
}


#endif

