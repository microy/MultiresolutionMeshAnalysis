/***************************************************************************
                                Curvature.cpp
                             -------------------
    update               : 2003-06-11
    copyright            : (C) 2002-2003 by Michaël Roy
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

#include "Curvature.h"
#include "VectorT.h"
#include <assert.h>
#include <math.h>


// Region Area
static double RegionArea(const Vector3d& vo, const Vector3d& va, const Vector3d& vb);


//--
//
// ComputeCurvatureNormals
//
//--
// Compute vertex normals according to the mean curvature
// Assume that the neighbors are collected
void ComputeCurvatureNormals(NeighborMesh* mesh)
{
	// Sanity checks
	assert( mesh );
	assert( mesh->VertexNeighborhoodNumber() == mesh->VertexNumber() );
	assert( mesh->FaceNeighborhoodNumber() == mesh->VertexNumber() );
	
	// Allocate vertex normal array
	mesh->VertexNormals().resize( mesh->VertexNumber() );
	
	// Temporary variables
	ConstNeighborIterator it1, it2;
	double coef;
	Vector3d average_normal;
	
	// Compute curvature normals for each vertex
	for( int i=0; i<mesh->VertexNumber(); i++ )
	{
		// Compute average of face normals in the neighborhood
		average_normal = 0;
		it1 = mesh->NeighborFaces(i).begin();
		while( it1 != mesh->NeighborFaces(i).end() )
		{
			average_normal += mesh->ComputeRawFaceNormal(*it1);
			++it1;
		}
		// Compute curvature normal if the vertex is not on a border edge
		mesh->VertexNormal(i) = 0;
		if( mesh->IsBorderVertex(i) == false )
		{
			it1 = mesh->NeighborVertices(i).begin();
			while( it1 != mesh->NeighborVertices(i).end() )
			{
				coef = 0;
				it2 = mesh->NeighborVertices(*it1).begin();
				while( it2 != mesh->NeighborVertices(*it1).end() )
				{
					if( mesh->FindNeighborVertex(i, *it2) )
					{
						coef += Cotan( mesh->Vertex(*it2), mesh->Vertex(i), mesh->Vertex(*it1) );
					}
					++it2;
				}
				mesh->VertexNormal(i) += (coef * (mesh->Vertex(i) - mesh->Vertex(*it1)));
				++it1;
			}
		}
		// Curvature normal can be 0 for flat plane or saddle point
		if( mesh->VertexNormal(i) == 0.0 )
		{
			mesh->VertexNormal(i) = average_normal;
		}
		// Check curvature normal orientation
		else if( (mesh->VertexNormal(i) | average_normal) < 0.0 )
		{
			mesh->VertexNormal(i) *= -1.0;
		}
		// Normalize the final vector
		mesh->VertexNormal(i).Normalize();
	}
}

//--
//
// ComputeMeanCurvatureNormal
//
//--
// Compute vertex mean curvature normal
Vector3d MeanCurvatureNormal(const NeighborMesh* mesh, int v)
{
	// Sanety checks
	assert( mesh );
	assert( (v>=0) && (v<mesh->VertexNumber()) );
	assert( mesh->VertexNeighborhoodNumber() == mesh->VertexNumber() );
	assert( mesh->FaceNeighborhoodNumber() == mesh->VertexNumber() );

	ConstNeighborIterator it1, it2;
	double area = 0;
	Vector3d curvature(0,0,0);
	double coef;

	it1 = mesh->NeighborFaces(v).begin();
	while( it1 != mesh->NeighborFaces(v).end() )
	{
		if( mesh->Face(*it1, 0) == v )
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 1), mesh->Vertex(*it1, 2));
		}
		else if( mesh->Face(*it1, 1) == v )
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 2), mesh->Vertex(*it1, 0));
		}
		else
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 0), mesh->Vertex(*it1, 1));
		}
		++it1;
	}

	// Sanity check
	assert( area != 0.0 );

	// Compute mean curvature normal
	it1 = mesh->NeighborVertices(v).begin();
	while( it1 != mesh->NeighborVertices(v).end() )
	{
		coef = 0;
		it2 = mesh->NeighborVertices(*it1).begin();
		while( it2 != mesh->NeighborVertices(*it1).end() )
		{
			if( mesh->FindNeighborVertex(v, *it2) )
			{
				coef += Cotan( mesh->Vertex(*it2), mesh->Vertex(v), mesh->Vertex(*it1) );
			}
			++it2;
		}
		curvature += coef * (mesh->Vertex(v)-mesh->Vertex(*it1));
		++it1;
	}

	area = 1.0 / area;
	area *= 0.5;
	return curvature * area;
}

//--
//
// ComputeGaussianNormal
//
//--
// Compute vertex gaussian curvature
double GaussianCurvature(const NeighborMesh* mesh, int v)
{
	assert( mesh );
	assert( (v>=0) && (v<mesh->VertexNumber()) );
	assert( mesh->VertexNeighborhoodNumber() == mesh->VertexNumber() );
	assert( mesh->FaceNeighborhoodNumber() == mesh->VertexNumber() );

	ConstNeighborIterator it1;
	double area = 0;
	double angle_sum = 0;
	
	it1 = mesh->NeighborFaces(v).begin();
	while( it1 != mesh->NeighborFaces(v).end() )
	{
		if( mesh->Face(*it1, 0) == v )
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 1), mesh->Vertex(*it1, 2));
			angle_sum += AngleFromCotan(mesh->Vertex(v), mesh->Vertex(*it1, 1), mesh->Vertex(*it1, 2));
		}
		else if( mesh->Face(*it1, 1) == v )
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 2), mesh->Vertex(*it1, 0));
			angle_sum += AngleFromCotan(mesh->Vertex(v), mesh->Vertex(*it1, 2), mesh->Vertex(*it1, 0));
		}
		else
		{
			area += RegionArea(mesh->Vertex(v), mesh->Vertex(*it1, 0), mesh->Vertex(*it1, 1));
			angle_sum += AngleFromCotan(mesh->Vertex(v), mesh->Vertex(*it1, 0), mesh->Vertex(*it1, 1));
		}
		++it1;
	}

	// Sanity check
	assert( area != 0.0 );

	return (2.0*M_PI - angle_sum) / area;
}

//--
//
// RegionArea
//
//--
double RegionArea(const Vector3d& vo, const Vector3d& va, const Vector3d& vb)
{
	double face_area = ((va-vo)^(vb-vo)).Length() * 0.5;

	if( face_area == 0.0 ) return 0.0;

	if( ObtuseAngle(vo, va, vb) )
	{
		return face_area * 0.5;
	}
	if( ObtuseAngle(va, vb, vo) || ObtuseAngle(vb, vo, va) )
	{
		return face_area * 0.25;
	}
	
	return (Cotan(va, vo, vb)*(vo-vb).SquareLength() + Cotan(vb, vo, va)*(vo-va).SquareLength()) * 0.125;
}

