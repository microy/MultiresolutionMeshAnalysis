/***************************************************************************
                              PairContraction.h
                             -------------------
    update               : 2002-09-04
    copyright            : (C) 2002 by Michaël Roy
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

#ifndef _PAIRCONTRACTION_
#define _PAIRCONTRACTION_

#include <assert.h>
#include <float.h>
#include <vector>

//
// PairType
// 
enum PairType
{
	//
	// Different type of vertex pair contraction
	// 
	NORMAL_PAIR,
	BORDER_PAIR,
	MULTIPLE_FACE_PAIR,
	SINGLE_VERTEX_PAIR
};

//
// PairContraction
//
class PairContraction
{
 //
	// This class represent a pair of vertices going to be collapsed.
	// That is represented by a candidate vertex, a target vertex,
	// a collapse cost, common faces shared by candidate and target vertex,
	// and a pair contraction type
	//
	public :

		//
		// Constructor
		//
		inline PairContraction( const double& c = DBL_MAX, const int& v = -1, const int& t = -1, const PairType& pt = NORMAL_PAIR )
  			: cost(c), candidate(v), target(t), type( pt) {
		}

		//
		// Cost Interface
		//
		inline double& Cost() {
			return cost;
		}
		
		inline const double& Cost() const {
			return cost;
		}
		
		//
		// Candidate Interface
		// 
		inline int& Candidate() {
			return candidate;
		}
		
		inline const int& Candidate() const {
			return candidate;
		}

		//
		// Target Interface
		// 
		inline int& Target() {
			return target;
		}
		
		inline const int& Target() const {
			return target;
		}
		
		//
		// Type interface
		//
		inline PairType& Type() {
			return type;
		}

		inline const PairType& Type() const {
			return type;
		}

		//
		// Common Face Interface
		//
		inline int CommonFaceNumber() const {
			return (int)common_faces.size();
		}

		inline std::vector<int>& CommonFaces() {
			return common_faces;
		}

		inline const std::vector<int>& CommonFaces() const {
			return common_faces;
		}

		inline int& CommonFace(int i) {
			assert( (i>=0) && (i<CommonFaceNumber()) );
			return common_faces[i];
		}

		inline const int& CommonFace(int i) const {
			assert( (i>=0) && (i<CommonFaceNumber()) );
			return common_faces[i];
		}

		inline void AddCommonFace(int i) {
			common_faces.push_back(i);
		}

		inline void ClearCommonFaces() {
			common_faces.clear();
		}

	protected :

		double           cost;
		int              candidate;
		int              target;
		PairType         type;
		std::vector<int> common_faces;

};

//
// Comparison Operators
//
inline bool operator<(const PairContraction& a, const PairContraction& b)
{
	if( a.Cost() < b.Cost() ) return true;
	if( a.Cost() == b.Cost() ) return a.Candidate() != b.Candidate();
	return false;
}
		
inline bool operator==(const PairContraction& a, const PairContraction& b)
{
	return a.Candidate() == b.Candidate();
}

#endif // _PAIRCONTRACTION_

