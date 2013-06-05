/***************************************************************************
                              NeighborMesh.cpp
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

#include "NeighborMesh.h"
#include <iostream>


//--
//
// CollectNeighbors
//
//--
// Collect neighbors
void NeighborMesh::CollectNeighbors()
{
	// Resize vertex neighborhood array
	neighbor_vertices.resize( VertexNumber() );

	// Resize face neighborhood array
	neighbor_faces.resize( VertexNumber() );
	
	// For every face
	for( int i=0; i<FaceNumber(); i++ )
	{
		// Shortcuts to vertex indices of the current face
		const int& a = Face( i, 0 );
		const int& b = Face( i, 1 );
		const int& c = Face( i, 2 );
		
		// Add neighbor face
		AddNeighborFace( a, i );
		AddNeighborFace( b, i );
		AddNeighborFace( c, i );

		// Add neighbor vertices
		AddNeighborVertex( a, b );
		AddNeighborVertex( a, c );
		
		AddNeighborVertex( b, a );
		AddNeighborVertex( b, c );
		
		AddNeighborVertex( c, a );
		AddNeighborVertex( c, b );
	}
}

//--
//
// CheckNeighborhood
//
//--
// Verify if neighbor links are correct
bool NeighborMesh::CheckNeighborhood(int v, bool verbose) const
{
	// Use Standard C++ Library output stream
	using std::cerr;
	using std::endl;
	
	// No neighbors ?
	if( (NeighborVertexNumber(v) == 0) && (NeighborFaceNumber(v) == 0) )
	{
		return true;
	}

	// Wrong neighbor number ?
	if( (NeighborVertexNumber(v) == 0) || (NeighborFaceNumber(v) == 0) )
	{
		if( verbose )
		{
			cerr<<"\n\nNeighborhood check: wrong neighbor number\n"<<endl;
		}
		return false;
	}

	ConstNeighborIterator itv;
	ConstNeighborIterator itf;
	
	//
	// Check neighbor vertices
	// 
	itv = NeighborVertices(v).begin();
	while( itv != NeighborVertices(v).end() )
	{
		bool vertex_found(false);
		itf = NeighborFaces(v).begin();
		while( itf != NeighborFaces(v).end() )
		{
			if( FaceHasVertex(*itf, *itv) )
			{
				vertex_found = true;
				break;
			}
			
			++itf;
		}
		
		if( vertex_found == false )
		{
			if( verbose )
			{
				cerr<<"\n\nNeighborhood check: neighbor vertex "<<*itv<<" not found\n"<<endl;
			}
			return false;
		}
		
		++itv;
	}
	
	itv = NeighborVertices(v).begin();
	while( itv != NeighborVertices(v).end() )
	{
		if( FindNeighborVertex(*itv, v) == false )
		{
			if( verbose )
			{
				cerr<<"\n\nNeighborhood check: vertex not found in neighbor vertex "<<*itv<<"\n"<<endl;
			}
			return false;
		}
		++itv;
	}
	
	//
	// Check neighbor faces
	// 
	itf = NeighborFaces(v).begin();
	while( itf != NeighborFaces(v).end() )
	{
		const int& va = Face(*itf, 0);
		const int& vb = Face(*itf, 1);
		const int& vc = Face(*itf, 2);
		
		if( (va==vb) || (va==vc) || (vb==vc) )
		{
			if( verbose )
			{
				cerr<<"\n\nNeighborhood check: invalid face "<<*itf<<"\n"<<endl;
			}
			return false;
		}
		
		if( (va<0) || (vb<0) || (vc<0) )
		{
			if( verbose )
			{
				cerr<<"\n\nNeighborhood check: invalid face "<<*itf<<"\n"<<endl;
			}
			return false;
		}
		
		bool vertex_1_found(false);
		bool vertex_2_found(false);

		if( va == v )
		{
			vertex_1_found = FindNeighborVertex(v, vb);
			vertex_2_found = FindNeighborVertex(v, vc);
		}
		else if( vb == v )
		{
			vertex_1_found = FindNeighborVertex(v, va);
			vertex_2_found = FindNeighborVertex(v, vc);
		}
		else if( vc == v )
		{
			vertex_1_found = FindNeighborVertex(v, va);
			vertex_2_found = FindNeighborVertex(v, vb);
		}
		
		if( (vertex_1_found==false) || (vertex_2_found==false) )
		{
			if( verbose )
			{
				cerr<<"\n\nNeighborhood check: vertex not found in face "<<*itf<<"\n"<<endl;
			}
			return false;
		}

		++itf;
	}

	return true;
}

//--
//
// IsBorderVertex
//
//--
// Return true if vertex is on border edge
// Assume neighborhood of the given vertex is valid
bool NeighborMesh::IsBorderVertex(int v, bool verbose) const
{
	// Use Standard C++ Library output stream
	using std::cerr;
	using std::endl;

	assert( (v>=0) && (v<VertexNumber()) );
	assert( VertexNeighborhoodNumber() == VertexNumber() );
	assert( FaceNeighborhoodNumber() == VertexNumber() );

	// Check neighbor vertex number
//	if( neighbor_vertices[v].size() < 3 ) {
//		if( verbose ) cerr<<"Vertex is on the border: less than 3 neighbor vertices"<<endl;
//		return true;
//	}

	// Check neighbor face number
//	if( neighbor_faces[v].size() < 3 ) {
//		if( verbose ) cerr<<"Vertex is on the border: less than 3 neighbor faces"<<endl;
//		return true;
//	}


	ConstNeighborIterator itv = neighbor_vertices[v].begin();
	while( itv != neighbor_vertices[v].end() ) {
		int common_face = 0;
		ConstNeighborIterator itf = neighbor_faces[*itv].begin();
		while( itf != neighbor_faces[*itv].end() ) {
			if( FaceHasVertex(*itf, v) == true ) common_face++;
			++itf;
		}
		if( common_face < 2 ) {
			if( verbose ) cerr<<"Vertex is on the border: less than 2 common faces"<<endl;
			return true;
		}
		++itv;
	}
	return false;
}

