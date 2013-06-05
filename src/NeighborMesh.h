/***************************************************************************
                               NeighborMesh.h
                             -------------------
    update               : 2004/01/21
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

#ifndef _NEIGHBORMESH_
#define _NEIGHBORMESH_

#include "Mesh.h"
#include <assert.h>
#include <algorithm>
#include <list>
#include <vector>

//--
//
// Type definitions
//
//--

// Neighbor list
typedef std::list<int> NeighborList;

// Neighbor iterator
typedef std::list<int>::iterator NeighborIterator;

// Constant neighbor iterator
typedef std::list<int>::const_iterator ConstNeighborIterator;


//--
//
// NeighborMesh
//
//--
class NeighborMesh : public Mesh
{
	//--
	//
	// Member Data
	//
	//--
	protected :

		// Neighbor faces of every vertex
		std::vector<std::list<int> > neighbor_faces;

		// Neighbor vertices of every vertex
		std::vector<std::list<int> > neighbor_vertices;
		
	//--
	//
	// Member Functions
	//
	//--
	public :

		//--
		//
		// Constructor / Destructor
		//
		//--
		inline NeighborMesh() {
		}

		inline ~NeighborMesh() {
		}

		//--
		//
		// Neighborhood management
		//
		//--

		// Collect neighbor faces and vertices
		void CollectNeighbors();

		// Check neighborhood of vertex #v
		// Return true if neighborhood of v is ok
		bool CheckNeighborhood(int v, bool verbose=false) const;

		// Is border vertex
		bool IsBorderVertex(int v, bool verbose=false) const;
	
		//--
		//
		// Neighbor Vertex Interface
		//
		//--

		// Vertex neighborhood number
		inline int VertexNeighborhoodNumber() const {
			return (int)neighbor_vertices.size();
		}

		// Vertex neighborhood array
		inline std::vector<std::list<int> >& VertexNeighborhoods() {
			return neighbor_vertices;
		}

		// Vertex neighborhood array (constant)
		inline const std::vector<std::list<int> >& VertexNeighborhoods() const {
			return neighbor_vertices;
		}

		// Neighbor vertex number of vertex #i
		inline int NeighborVertexNumber(int i) const {
			assert( (i>=0) && (i<VertexNeighborhoodNumber()) );
			return (int)neighbor_vertices[i].size();
		}

		// Neighbor vertex list of vertex #i
		inline std::list<int>& NeighborVertices(int i) {
			assert( (i>=0) && (i<VertexNeighborhoodNumber()) );
			return neighbor_vertices[i];
		}

		// Neighbor vertex list of vertex #i (constant)
		inline const std::list<int>& NeighborVertices(int i) const {
			assert( (i>=0) && (i<VertexNeighborhoodNumber()) );
			return neighbor_vertices[i];
		}

		// Find vertex #i in vertex neighborhood of vertex #v
		inline bool FindNeighborVertex( int v, int i ) const {
			assert( (v>=0) && (v<VertexNeighborhoodNumber()) );
			return std::find(neighbor_vertices[v].begin(), neighbor_vertices[v].end(), i) != neighbor_vertices[v].end();
		}

		// Add vertex #i in vertex neighborhood of vertex #v
		inline void AddNeighborVertex( int v, int i ) {
			assert( (v>=0) && (v<VertexNeighborhoodNumber()) );
			if( std::find(neighbor_vertices[v].begin(), neighbor_vertices[v].end(), i) == neighbor_vertices[v].end() )
				neighbor_vertices[v].push_back(i);
		}

		// Remove vertex #i in vertex neighborhood of vertex #v
		inline void RemoveNeighborVertex( int v, int i ) {
			assert( (v>=0) && (v<VertexNeighborhoodNumber()) );
			neighbor_vertices[v].remove(i);
		}

		// Clear vertex neighborhood of vertex #v
		inline void ClearNeighborVertices(int v) {
			assert( (v>=0) && (v<VertexNeighborhoodNumber()) );
			neighbor_vertices[v].clear();
		}

		//--
		//
		// Neighbor Face Interface
		//
		//--
		
		// Face neighborhood number
		inline int FaceNeighborhoodNumber() const {
			return (int)neighbor_faces.size();
		}

		// Face neighborhood array
		inline std::vector<std::list<int> >& FaceNeighborhoods() {
			return neighbor_faces;
		}

		// Face neighborhood array (constant)
		inline const std::vector<std::list<int> >& FaceNeighborhoods() const {
			return neighbor_faces;
		}

		// Neighbor face number of vertex #i
		inline int NeighborFaceNumber(int i) const {
			assert( (i>=0) && (i<FaceNeighborhoodNumber()) );
			return (int)neighbor_faces[i].size();
		}

		// Neighbor face list of vertex #i
		inline std::list<int>& NeighborFaces(int i) {
			assert( (i>=0) && (i<FaceNeighborhoodNumber()) );
			return neighbor_faces[i];
		}
		
		// Neighbor face list of vertex #i (constant)
		inline const std::list<int>& NeighborFaces(int i) const {
			assert( (i>=0) && (i<FaceNeighborhoodNumber()) );
			return neighbor_faces[i];
		}

		// Find face #i in face neighborhood of vertex #v
		inline bool FindNeighborFace(int v, int i) const {
			assert( (v>=0) && (v<FaceNeighborhoodNumber()) );
			return std::find(neighbor_faces[v].begin(), neighbor_faces[v].end(), i) != neighbor_faces[v].end();
		}

		// Add face #i in face neighborhood of vertex #v
		inline void AddNeighborFace(int v, int i) {
			assert( (v>=0) && (v<FaceNeighborhoodNumber()) );
			if( std::find(neighbor_faces[v].begin(), neighbor_faces[v].end(), i) == neighbor_faces[v].end() )
				neighbor_faces[v].push_back(i);
		}
        
		// Remove face #i in face neighborhood of vertex #v
		inline void RemoveNeighborFace(int v, int i) {
			assert( (v>=0) && (v<FaceNeighborhoodNumber()) );
			neighbor_faces[v].remove(i);
		}

		// Clear face neighborhood of vertex #v
		inline void ClearNeighborFaces(int v) {
			assert( (v>=0) && (v<FaceNeighborhoodNumber()) );
			neighbor_faces[v].clear();
		}
};

#endif // _NEIGHBORMESH_

