/***************************************************************************
                              ResolutionLevel.h
                             -------------------
    update               : 2003/10/24
    copyright            : (C) 2003 by Michaël Roy
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

#ifndef _RESOLUTIONLEVEL_
#define _RESOLUTIONLEVEL_

#include "PairContraction.h"
#include "VectorT.h"
#include <list>
#include <vector>
#include <map>

//--
//
// SubdivisionWeight
//
//--
struct SubdivisionWeight
{
	//
	// Member Data
	//
	int vertex;
	double weight;

	//
	// Constructor
	//
	SubdivisionWeight(int v = -1, double w = 0) : vertex(v), weight(w) {}
};

typedef std::list<SubdivisionWeight> SubdivisionWeightList;

// Pair contraction list
typedef std::list<PairContraction> PairContractionList;

// Pair contraction iterator
typedef std::list<PairContraction>::iterator PairContractionIterator;

// Pair contraction constant iterator
typedef std::list<PairContraction>::const_iterator ConstPairContractionIterator;


//--
//
// DetailType
//
//--
enum VertexType
{
	UNDEFINED_VERTEX,
	EVEN_VERTEX,
	ODD_VERTEX
};


//--
//
// ResolutionLevel
//
//--
class ResolutionLevel
{
	//
	// Member Data
	//
	public :

		std::map<int,Vector3d> geometric_details;
		std::map<int,Vector3d> normal_details;
		std::map<int,Vector3d> color_details;
		std::map<int,Vector3d> texture_details;

		std::vector<VertexType> vertex_types;

		std::vector<SubdivisionWeightList> subdivision_weights;

		std::list<PairContraction> pair_contractions;

	//	std::map<Vector2i, double> weights;

	//
	// Member Functions
	//
	public :

		//
		// Constructor / Destructor
		//
		ResolutionLevel() {}
		~ResolutionLevel() {}

		//
		void Resize( int size );
};



#endif // _RESOLUTIONLEVEL_

