/***************************************************************************
                             ProgressiveMesh.cpp
                             -------------------
    update               : 2004/05/18
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

#include "ProgressiveMesh.h"
#include <algorithm>
#include <float.h>


//
// BeginProgressiveReconstruction
// 
void ProgressiveMesh::BeginProgressiveDecomposition()
{
	pair_contractions.clear();
	pair_contraction_indices.resize( VertexNumber(), pair_contractions.end() );

	valid_faces.resize(FaceNumber(), true);
	valid_vertices.resize(VertexNumber(), true);

	invalid_contraction_penalty = 1e9;
	border_penalty = 1e3;

	// Collect information
	CollectNeighbors();
	CollectQuadrics();
	ComputeAllEdgeCollapseCosts();
}

//
// EndProgressiveDecomposition
// 
void ProgressiveMesh::EndProgressiveDecomposition()
{
	pair_contraction_indices.clear();
	pair_contractions.clear();
	vertex_quadrics.clear();
}


void ProgressiveMesh::SetAllValid()
{
	valid_faces.resize(FaceNumber(), true);
	valid_vertices.resize(VertexNumber(), true);
}

//--
//
// IsValidPairContraction
//
//--
bool ProgressiveMesh::IsValidPairContraction(int va, int vb) const
{
	assert(va>=0);
	assert(vb>=0);

	// Test to see if after the pair contraction,
	// one vertex has a valence of 3
	// This case screws up the predict operator
//	ConstNeighborIterator it = NeighborVertices(va).begin();
//	while( it != NeighborVertices(va).end() ) {
//		if(NeighborVertexNumber(*it)==4 && FindNeighborVertex(vb, *it) && !IsBorderVertex(*it)) {
//			return false;
//		}
//		++it;
//	}			

	// Check for invalid faces after pair contraction
	ConstNeighborIterator it = NeighborFaces(va).begin();
	while( it != NeighborFaces(va).end() ) {

		// Process only faces that are not common to va and vb
		if( !FaceHasVertex( *it, vb ) )	{

			// Find the two other vertices making the face
			int o1, o2;
			if( Face(*it,0) == va )	{
				o1 = Face(*it,1);
				o2 = Face(*it,2);
			}
			else if( Face(*it,1) == va ) {
				o1 = Face(*it,2);
				o2 = Face(*it,0);
			}
			else {
				o1 = Face(*it,0);
				o2 = Face(*it,1);
			}

			// Check degenerate face
			Vector3d u = Vertex(o1) - Vertex(vb);
			Vector3d v = Vertex(o2) - Vertex(vb);
			double angle = (u | v) / (u.Length()*v.Length());
			if( (angle<-0.99) || (angle>0.99) ) return false;

			// Face normal before contraction
			Vector3d current_normal = ComputeFaceNormal( *it );

			// Compute normal of the new face after contraction
			Vector3d new_normal = ComputeFaceNormal( vb, o1, o2 );

			// Compute dot product (cos angle) between normals of
			// the previous face and the new face
			// If the dot product is less than 0 (angle<pi/2)
			// the new face is inverted
			if( (current_normal | new_normal) < 0.0 ) return false;
		}
		++it;
	}

	// Valid pair contraction (va->vb)
	return true;
}


//--
//
// ComputeEdgeCollapseCost
//
//--
double ProgressiveMesh::ComputeEdgeCollapseCost( int va, int vb )
{
	// Sanity check
	assert(va>=0);
	assert(vb>=0);

	// Compute the contraction cost
//	Quadric3 quad( VertexQuadric(va) );
//	quad += VertexQuadric(vb);
//	double error = quad.Evaluate( Vertex(vb) );
	
	// TEST: edge length as the error metric
	double error = (Vertex(va) - Vertex(vb)).Length();
	
	// TEST: random vertex selection (i.e. random error value)
//	static bool first_time = true;
//	if( first_time == true ) {
		// Initialize the random seed
//		srand( time(NULL) );
//		first_time = false;
//	}
//	error = (double)rand() / (double)RAND_MAX;

	// Border test
	if( IsBorderVertex(va) && !IsBorderVertex(vb) )	return error + invalid_contraction_penalty;
//	if( !IsBorderVertex(va) && IsBorderVertex(vb) )	return error + invalid_contraction_penalty;

	// Check for invalid faces after pair contraction
	if( !IsValidPairContraction(va, vb) ) return error + invalid_contraction_penalty;

	// Pair conctraction error
	return error;
}


//
// ComputeEdgeCostAtVertex
//
void ProgressiveMesh::ComputeEdgeCostAtVertex( int v )
{
	// Security check
	assert( CheckNeighborhood(v) );

	if( NeighborFaceNumber(v) == 0 )
	{
		// v doesn't have neighbors so it costs nothing to collapse
		AddPairContraction( v, PairContraction(-0.01, v, -1, SINGLE_VERTEX_PAIR) );
		return;
	}


	double tmp_cost;
	PairContraction pair( DBL_MAX, v );
	
	// search all neighboring edges for "least cost" edge
	std::list<int>::iterator it = NeighborVertices(v).begin();
	while( it != NeighborVertices(v).end() )
	{
		// Security checks
		assert( NeighborFaceNumber(*it) > 0 );
		assert( FindNeighborVertex(*it, v) == true );
		// Compute edge v-it collapse cost 
		tmp_cost = ComputeEdgeCollapseCost(v, *it );

		// Record lowest collapse cost
		if( tmp_cost < pair.Cost() ) {
			// candidate for edge collapse
			pair.Target() = *it;
			// cost of the collapse
			pair.Cost() = tmp_cost;
		}
		++it;
	}
	
	// Record faces on edge va-vb
	it = NeighborFaces(pair.Candidate()).begin();
	while( it != NeighborFaces(pair.Candidate()).end() )
	{
		if( FaceHasVertex( *it, pair.Target() ) ) {
			pair.AddCommonFace(*it);
		}
		++it;
	}

	// Only valid pair must be added
	assert( pair.CommonFaceNumber() > 0 );

	if( pair.CommonFaceNumber() == 1 ) pair.Type() = BORDER_PAIR;
	else if( pair.CommonFaceNumber() > 2 ) pair.Type() = MULTIPLE_FACE_PAIR;

	// Add pair to the pair constraction set
	AddPairContraction( v, pair );

}


void ProgressiveMesh::ComputeAllEdgeCollapseCosts()
{
	// For all the edges, compute the difference it would make
	// to the model if it was collapsed.  The least of these
	// per vertex is cached in each vertex object.
	for(int i=0; i<VertexNumber(); i++ )
	{
		ComputeEdgeCostAtVertex(i);
	}
}

void ProgressiveMesh::RemoveFace( int f )
{
	RemoveNeighborFace( Face(f, 0), f );
	RemoveNeighborFace( Face(f, 1), f );
	RemoveNeighborFace( Face(f, 2), f );
}

void ProgressiveMesh::InsertFace( int f )
{
	const int& va = Face( f, 0 );
	const int& vb = Face( f, 1 );
	const int& vc = Face( f, 2 );
	
	AddNeighborFace( va, f );
	AddNeighborFace( vb, f );
	AddNeighborFace( vc, f );

	AddNeighborVertex( va, vb );
	AddNeighborVertex( va, vc );
	
	AddNeighborVertex( vb, va );
	AddNeighborVertex( vb, vc );
		
	AddNeighborVertex( vc, va );
	AddNeighborVertex( vc, vb );
}

void ProgressiveMesh::ReplaceVertex( int f, int vold, int vnew )
{
	// Check
	assert( FaceHasVertex(f, vold) );
	
	if( Face(f, 0) == vold )
	{
		Face(f, 0) = vnew;
	}
	else if( Face(f, 1) == vold )
	{
		Face(f, 1) = vnew;
	}
	else
	{
		Face(f, 2) = vnew;
	}
	AddNeighborFace(vnew, f);
}

void ProgressiveMesh::RemoveVertex( int v )
{
	std::list<int>::iterator it = NeighborVertices(v).begin();
	while( it != NeighborVertices(v).end() )
	{
		RemoveNeighborVertex( *it, v );
		++it;
	}
}

void ProgressiveMesh::InsertVertex( int v )
{
	std::list<int>::iterator it = NeighborVertices(v).begin();
	while( it != NeighborVertices(v).end() )
	{
		AddNeighborVertex( *it, v );
		++it;
	}
}

//--
//
// Collapse
//
//--
// Collapse a pair of vertices
void ProgressiveMesh::Collapse( const PairContraction& pair )
{
	//
	// Collapse a vertex pair. This is an half-edge collapse:
	// the vb vertex is the result of the contraction.
	//
	
	const int& va = pair.Candidate();
	const int& vb = pair.Target();
	std::list<int>::iterator it;

	// Sanity checks
	assert( CheckNeighborhood(va) );
	assert( (vb<0) ? true : CheckNeighborhood(vb) );
	assert( (vb<0) ? true : FindNeighborVertex(va, vb) );
	assert( (vb<0) ? true : FindNeighborVertex(vb, va) );

	// Remove vertex va
	RemoveVertex( va );
	valid_vertices[va] = false;

	// Remove a single vertex
	if( vb < 0 )
	{
		assert( NeighborFaceNumber(va) == 0 );
		assert( NeighborVertexNumber(va) == 0 );
		return;
	}

	// Delete faces on edge va-vb
	for( int i=0; i<pair.CommonFaceNumber(); i++ )
	{
		RemoveFace( pair.CommonFace(i) );
		valid_faces[pair.CommonFace(i)] = false;
	}

	// Update neighbor face links
	it = NeighborFaces(va).begin();
	while( it != NeighborFaces(va).end() )
	{
		ReplaceVertex(*it, va, vb);
		assert( FaceHasVertex(*it, vb) );
		++it;
	}
	
	// If vb has no more neighbor face
	if( NeighborFaceNumber(vb) == 0 )
	{
		// Remove vb
		RemoveVertex(vb);
		// Update previous neighbor vertices
		it = NeighborVertices(vb).begin();
		while( it != NeighborVertices(vb).end() )
		{
			++it;
		}
		// Update vb
		ClearNeighborVertices(vb);
		return;
	}
	
	// Look for single vertices
	std::list<int> alone_vertices;
	it = NeighborVertices(va).begin();
	while( it != NeighborVertices(va).end() )
	{
		// If it has no neighbor face
		// this is a single vertex
		if( NeighborFaceNumber(*it) == 0 )
		{
			alone_vertices.push_back(*it);
		}
		++it;
	}

	// Remove single vertices
	it = alone_vertices.begin();
	while( it != alone_vertices.end() )
	{
		RemoveVertex(*it);
		ClearNeighborVertices(*it);
		++it;
	}

	// Update neighbor vertex links
	it = NeighborVertices(va).begin();
	while( it != NeighborVertices(va).end() )
	{
		if( *it != vb )
		{
			// Remove this neighbor (sanity)
			RemoveNeighborVertex(*it, vb);
			RemoveNeighborVertex(vb, *it);
			std::list<int>::iterator itf = NeighborFaces(*it).begin();
			while( itf != NeighborFaces(*it).end() )
			{
				// Check if vb has a common face with *it
				if( FaceHasVertex( *itf, vb ) )
				{
					// Link vertices *it & vb
					AddNeighborVertex(*it, vb);
					AddNeighborVertex(vb, *it);
					break;
				}
				++itf;
			}
		}
		++it;
	}
	
	// Check vb neighborhood
	assert( CheckNeighborhood(vb) );
}

//--
//
// Expand
//
//--
void ProgressiveMesh::Expand( const PairContraction& pair )
{
	//
	// This function inverts a pair contraction
	// + Insert vertex va
	// + Update links of vb
	//

	const int& va = pair.Candidate();
	const int& vb = pair.Target();

	valid_vertices[va] = true;

	if( vb < 0 )
	{
		InsertVertex( va );
		valid_vertices[va] = true;
		return;
	}

	assert( CheckNeighborhood(vb) );

	// Update remaining faces to have vertex_a instead of vertex_b
	std::list<int>::iterator it = NeighborFaces(va).begin();
	while( it != NeighborFaces(va).end() )
	{
		ReplaceVertex(*it, vb, va);
		RemoveNeighborFace(vb, *it);
		++it;
	}

	// Add faces on edge va-vb
	for( int i=0; i<pair.CommonFaceNumber(); i++ )
	{
		InsertFace( pair.CommonFace(i) );
		valid_faces[ pair.CommonFace(i) ] = true;
	}

	it = NeighborVertices(va).begin();
	while( it != NeighborVertices(va).end() )
	{
		if( *it != vb )
		{
			RemoveNeighborVertex( *it, vb );
			RemoveNeighborVertex( vb, *it );
			std::list<int>::iterator itf = NeighborFaces(*it).begin();
			while( itf != NeighborFaces(*it).end() )
			{
				// Check if vb has a common face with *it
				if( FaceHasVertex( *itf, vb ) )
				{
					// Link vertices *it & vb
					AddNeighborVertex(*it, vb);
					AddNeighborVertex(vb, *it);
					break;
				}
				++itf;
			}
		}
		++it;
	}

	InsertVertex( va );
	valid_vertices[va] = true;

	assert( CheckNeighborhood(va) );
	assert( CheckNeighborhood(vb) );
}

//--
//
// UpdateNormals
//
//--
// Update face and vertex normals on the neighborhood of a given pair after its contraction
void ProgressiveMesh::UpdateNormals( const PairContraction& pair )
{
	const int& v = pair.Target();

	// Single vertex removed
	if( v < 0 ) return;
	
	ConstNeighborIterator itv, itf;

	VertexNormal(v) = 0.0;
	itf = NeighborFaces(v).begin();
	while( itf != NeighborFaces(v).end() )
	{
		ComputeFaceNormal(*itf);
		VertexNormal(v) += FaceNormal(*itf);
		++itf;
	}
	VertexNormal(v).Normalize();
	
	itv = NeighborVertices(v).begin();
	while( itv != NeighborVertices(v).end() )
	{
		VertexNormal(*itv) = 0.0;
		itf = NeighborFaces(*itv).begin();
		while( itf != NeighborFaces(*itv).end() )
		{
			VertexNormal(*itv) += FaceNormal(*itf);
			++itf;
		}
		VertexNormal(*itv).Normalize();
		++itv;
	}
}



//--
//
// UpdateQuadrics
//
//--
// Update quadrics on the neighborhood of a given pair after its contraction
void ProgressiveMesh::UpdateQuadrics( const PairContraction& pair )
{
	const int& va = pair.Candidate();
	const int& vb = pair.Target();

	// Single vertex removed
	if( vb < 0 ) return;
	
	// Add candidate vertex quadric to target vertex quadric
	VertexQuadric(vb) += VertexQuadric(va);
}

//--
//
// UpdateEdgeCosts
//
//--
// Update edge costs
void ProgressiveMesh::UpdateEdgeCosts( const PairContraction& pair )
{
	const int& va = pair.Candidate();
	const int& vb = pair.Target();

	// Compute neighbor vertex cost
	ConstNeighborIterator it = NeighborVertices(va).begin();
	while( it != NeighborVertices(va).end() )
	{
		RemovePairContraction(*it);
		ComputeEdgeCostAtVertex(*it);
		++it;
	}

	// Single vertex removed
	if( vb < 0 ) return;
	
	// Compute neighbor vertex cost
	it = NeighborVertices(vb).begin();
	while( it != NeighborVertices(vb).end() )
	{
		RemovePairContraction(*it);
		ComputeEdgeCostAtVertex(*it);
		++it;
	}
		
	// Compute target vertex cost
	RemovePairContraction(vb);
	ComputeEdgeCostAtVertex(vb);
}

//
// CollectQuadrics
//
void ProgressiveMesh::CollectQuadrics()
{
	face_normals.resize( FaceNumber() );
	vertex_quadrics.resize( VertexNumber() );
	
	/*
	// Weight the quadric by the edge length
	for( int i=0; i<VertexNumber(); i++ ) {
		VertexQuadric(i) += Quadric3(1.0, 0.0, 0.0, -Vertex(i)[0],
		                                  1.0, 0.0, -Vertex(i)[1],
		                                       1.0, -Vertex(i)[2],
		                                             Vertex(i).Length());
	}
 	*/
 	
	for( int i=0; i<FaceNumber(); i++ ) {
		const int& va = Face(i,0);
		const int& vb = Face(i,1);
		const int& vc = Face(i,2);

		FaceNormal(i) = ComputeRawFaceNormal( va, vb, vc );
//	 	double area = 0.5 * FaceNormal(i).Length();
		FaceNormal(i).Normalize();

 		Quadric3 quad = Quadric3( FaceNormal(i), -(Vertex(va)|FaceNormal(i)) );

		// Area weighted quadric
//		quad *= quad.Area();
		
		VertexQuadric(va) += quad;
		VertexQuadric(vb) += quad;
		VertexQuadric(vc) += quad;
	}

	//
	// Border
	//
	PrepareDiscontinuities();
	
}


//--
//
// PrepareBorderEdge
//
//--
// Constrain the vertex quadric of the border edge (va-vb)
// Construct a perpendicular plane to the surface on the border edge,
// and add the plane quadric to each vertex quadric of the border edge.
void ProgressiveMesh::PrepareBorderEdge(int va, int vb, const Vector3d& normal)
{
	Vector3d dir = Vertex(vb) - Vertex(va);

	Vector3d v = (dir ^ normal).Normalize();
	Quadric3 q( v, -(v|Vertex(va)) );

	q *= border_penalty;
	
	// Area weighted quadric
	q.Area() = dir.SquareLength();
//	q3 *= q3.Area();
	 
	VertexQuadric(va) += q;
	VertexQuadric(vb) += q;
}


//--
//
// PrepareDiscontinuities
//
//--
// Find and constrain border edges.
void ProgressiveMesh::PrepareDiscontinuities()
{
	for( int i=0; i<VertexNumber(); i++ ) {
		ConstNeighborIterator itv = NeighborVertices(i).begin();
		while( itv != NeighborVertices(i).end() ) {
			if( i < *itv ) {
				std::list<int> sides;
				ConstNeighborIterator itf = NeighborFaces(i).begin();
				while( itf != NeighborFaces(i).end() ) {
					if( FaceHasVertex( *itf, *itv ) ) {
						sides.push_back(*itf);
					}
					++itf;
				}
				if( sides.size() == 1 ) {
					PrepareBorderEdge(i, *itv, FaceNormal(*sides.begin()));
				}
			}
			++itv;
		}
	}
}


//--
//
// ComputeProgressiveFaceNormals
//
//--
// Compute unit normal of every faces
void ProgressiveMesh::ComputeProgressiveFaceNormals()
{
	// Resize face normal array
	if( face_normals.size() != faces.size() ) {
		face_normals.resize( FaceNumber() );
	}

	// For every valid face
	for( int i=0; i<FaceNumber(); i++ ) {
		if( IsProgressiveFaceValid(i) == false ) continue;
		// Compute unit face normal
		FaceNormal(i) = ComputeFaceNormal(i);
	}
}


//--
//
// ComputeProgressiveVertexNormals
//
//--
// Compute unit normal of every vertex
void ProgressiveMesh::ComputeProgressiveVertexNormals()
{
	int i;

	// Assume that face normals are computed
	assert( FaceNormalNumber() == FaceNumber() );

	// Resize and initialize vertex normal array
	vertex_normals.assign( VertexNumber(), Vector3d(0,0,0) );
	
	// For every face
	for( i=0 ; i<FaceNumber() ; i++ ) {
		if( IsProgressiveFaceValid(i) == false ) continue;
		// Add face normal to vertex normal
		VertexNormal(i,0) += FaceNormal(i);
		VertexNormal(i,1) += FaceNormal(i);
		VertexNormal(i,2) += FaceNormal(i);
	}

	// For every vertex
	for( i=0 ; i<VertexNumber() ; i++) {
		if( IsProgressiveVertexValid(i) == false ) continue;
		// Normalize vertex normal
		VertexNormal(i).Normalize();
	}
}


//--
//
// WriteFile
//
//--
// Write progressive mesh into a file
bool ProgressiveMesh::WriteFile( const std::string& file_name, const FileFormat& file_format ) const
{
	Mesh tmp_mesh;
	std::vector<int> vertex_lut(VertexNumber());
	int count = 0;

	// Add valid vertices
	for( int i=0; i<VertexNumber(); i++ ) {
		if( valid_vertices[i] == false ) {
			vertex_lut[i] = -1;
			continue;
		}
		tmp_mesh.AddVertex( Vertex(i) );
		vertex_lut[i] = count++;
	}

	// Add valid colors
	if( ColorNumber() == VertexNumber() ) {
		for( int i=0; i<VertexNumber(); i++ ) {
			if( valid_vertices[i] == false ) continue;
			tmp_mesh.AddColor( Color(i) );
		}
	}

	// Add valid textures
	if( (TextureNumber()!=0) && (TextureName().empty()==false) ) {
		tmp_mesh.TextureName() = TextureName();
		for( int i=0; i<VertexNumber(); i++ ) {
			if( valid_vertices[i] == false ) continue;
			tmp_mesh.AddTexture( Texture(i) );
		}
	}

	// Remesh with valid vertices and valid faces
	for( int i=0; i<FaceNumber(); i++ )
	{
		if( valid_faces[i] == false ) continue;
		if( valid_vertices[Face(i,0)] == false ) continue;
		if( valid_vertices[Face(i,1)] == false ) continue;
		if( valid_vertices[Face(i,2)] == false ) continue;
		tmp_mesh.AddFace( Vector3i(vertex_lut[Face(i,0)], vertex_lut[Face(i,1)], vertex_lut[Face(i,2)]) );
	}

	// Compute normals
	tmp_mesh.ComputeFaceNormals();
	tmp_mesh.ComputeVertexNormals();

	// Write the progressive mesh
	return tmp_mesh.WriteFile(file_name, file_format);
}

