/***************************************************************************
                         MultiresolutionProcessing.h
                             -------------------
    update               : 2004/05/17
    copyright            : (C) 2003-2004 by Michaël Roy
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

#ifndef _MULTIRESOLUTIONPROCESSING_
#define _MULTIRESOLUTIONPROCESSING_

#include "MultiresolutionMesh.h"
#include "PairContraction.h"
#include "VectorT.h"
#include "Viewport.h"

#include <map>
#include <vector>


//--
//
// Enumerations
//
//--

// Segmentation type
enum SegmentationType
{
	SEGMENTATION_NONE = 0,
	SEGMENTATION_DILATION,
	SEGMENTATION_EROSION,
	SEGMENTATION_OPENING,
	SEGMENTATION_CLOSING
};

// Detail type
enum DetailType
{
	DETAIL_GEOMETRIC = 0,
	DETAIL_NORMAL,
	DETAIL_COLOR,
	DETAIL_TEXTURE
};



//--
//
// MultiresolutionProcessing
//
//--
class MultiresolutionProcessing
{

	//--
	//
	// Member variables
	//
	//--
	protected :

		MultiresolutionMesh* mesh;

	//	std::vector< std::vector<Vector3d> > original_geometric_details;
		std::vector< std::map<int, Vector3d> > original_geometric_details;
		std::vector< std::map<int, Vector3d> > original_color_details;
		std::vector<Vector3d> original_vertices;
		std::vector<Vector3i> original_faces;
		std::vector<Vector3d> original_colors;
		PairContractionList simplification_contractions;

	//--
	//
	// Member functions
	//
	//--
	public :

		// Default constructor
		MultiresolutionProcessing();

		// Constructor with a multiresolution mesh
		MultiresolutionProcessing(MultiresolutionMesh* m);

		// Destructor
		~MultiresolutionProcessing();

		// Set the multiresolution mesh to process
		bool Model(MultiresolutionMesh* m);

		// Reset the multiresolution mesh
		bool Reset();

		//
		// Multiresolution processings
		//
		bool Filter( const std::vector<double>& fb );
		bool Denoise( double threshold, int level_number, bool color=false );
		bool Threshold( double threshold, int level_number, SegmentationType type=SEGMENTATION_NONE, int iteration=1 );
		bool Color( DetailType type, int level_number, bool smooth=false );
		bool Simplify(double threshold, int level_number, DetailType type, bool view_dependent = false, Viewport vp = Viewport());
		bool Statistics();
		void Test();
		
	private :
	
		double GeometricDetailNoiseVariance();
		double GeometricDetailAngleNoiseVariance();
		double LocalGeometricDetailMedian(int i, int l);
		double GeometricDetailVarianceN1(int i, int l);
		double GeometricDetailVarianceN2(int i, int l);
		double GeometricDetailAngleVarianceN1(int i, int l);
};

#endif // _MULTIRESOLUTIONPROCESSING_
