/***************************************************************************
                              ProgressiveMesh.h
                             -------------------
    update               : 2004/06/18
    copyright            : (C) 2002-2004 by Michaël Roy
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

#ifndef _PROGRESSIVEMESH_
#define _PROGRESSIVEMESH_

#include "NeighborMesh.h"
#include "PairContraction.h"
#include "Quadric3.h"

#include <algorithm>
#include <list>
#include <set>
#include <vector>


//--
//
// ProgressiveMesh
//
//--
class ProgressiveMesh : public NeighborMesh
{

	//
	// Member Data
	//
	protected :

		//
		// Pair Contraction
		//
		std::set<PairContraction> pair_contractions;
		std::vector<std::set<PairContraction>::iterator > pair_contraction_indices;

		// Quadrics
		std::vector<Quadric3> vertex_quadrics;

		double invalid_contraction_penalty;
		double border_penalty;

		std::vector<int> map;
		std::vector<int> permutation;

		std::vector<bool> valid_faces;
		std::vector<bool> valid_vertices;

	//
	// Member Functions
	// 
	public :

		//
		// Constructor / Destructor
		//
		ProgressiveMesh() : vertex_quadrics(0) {
		}
		
		~ProgressiveMesh() {
		}

		//
		// Vertex Map Interface
		//
		inline int& Map(int i) {
			assert( (i>=0) && (i<(int)map.size()) );
			return map[i];
		}
		
		inline const int& Map(int i) const {
			assert((i>=0) && (i<(int)map.size()) );
  			return map[i];
		}

		// Write mesh in a file
		bool WriteFile( const std::string& file_name, const FileFormat& file_format=VRML_1_FILE ) const;

 		//
		// Initializations
		// 
		void BeginProgressiveDecomposition();
		void EndProgressiveDecomposition();

		//
		// Face and Vertex validities
		//
		void SetAllValid();

		bool IsProgressiveFaceValid(int f) const {
			assert( (f>=0) && (f<(int)valid_faces.size()) );
			return valid_faces[f];
		}

		bool IsProgressiveVertexValid(int v) const {
			assert( (v>=0) && (v<(int)valid_vertices.size()) );
			return valid_vertices[v];
		}

		int ValidVertexNumber() const {
			int count = 0;
			for( int i=0; i<VertexNumber(); i++ ) if( valid_vertices[i] ) count++;
			return count;
		}

		int ValidFaceNumber() const {
			int count = 0;
			for( int i=0; i<FaceNumber(); i++ ) if( valid_faces[i] ) count++;
			return count;
		}

		// Compute face normals
		void ComputeProgressiveFaceNormals();

		// Compute vertex normals
		// Assume that face normals are computed
		void ComputeProgressiveVertexNormals();


//	protected :

		//
		// Information
		//
		void CollectQuadrics();

		//
		// Cost Computation
		//
		void ComputeAllEdgeCollapseCosts();
  		double ComputeEdgeCollapseCost( int va, int vb);
		void ComputeEdgeCostAtVertex( int v );

		//
		// Checks
		//
		bool IsValidPairContraction(int va, int vb) const;
		inline bool IsValidPairContraction(const PairContraction& pair) const {
			return IsValidPairContraction(pair.Candidate(), pair.Target());
		}

		//
		// Simplification
		//
		void RemoveFace( int f );
		void ReplaceVertex( int f, int vold, int vnew );
		void RemoveVertex( int v );
		void Collapse( const PairContraction& pair );
		void UpdateNormals(const PairContraction& pair);
		void UpdateQuadrics(const PairContraction& pair);
		void UpdateEdgeCosts(const PairContraction& pair);

		//
		// Subdivision
		//
		void InsertVertex(int v);
		void InsertFace(int f);
		void Expand( const PairContraction& pair );

		//
		// Border Edge Handle
		//
		void PrepareBorderEdge(int va, int vb, const Vector3d& normal);
		void PrepareDiscontinuities();

   		//
		// Vertex Quadric Interface
		//
		int VertexQuadricNumber() const {
			return (int)vertex_quadrics.size();
		}

		Quadric3& VertexQuadric(int i) {
			assert( (i>=0) && (i<(int)vertex_quadrics.size()) );
			return vertex_quadrics[i];
		}

		//
		// Pair Contraction Interface
		//
		void AddPairContraction( int v, const PairContraction& p ) {
			assert( (v>=0) && (v<(int)pair_contraction_indices.size()) );
			assert( p.Candidate() == v );
			pair_contraction_indices[v] = pair_contractions.insert( p ).first;
		}

		void RemovePairContraction( int v ) {
			assert( (v>=0) && (v<(int)pair_contraction_indices.size()) );
			if( pair_contraction_indices[v] != pair_contractions.end() )
				pair_contractions.erase( pair_contraction_indices[v] );
			pair_contraction_indices[v] = pair_contractions.end();
		}
     	
		const PairContraction& MinimumPairContraction() const {
  			return *(pair_contractions.begin());
		}
};

#endif // _PROGRESSIVEMESH_

