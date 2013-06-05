/***************************************************************************
                        MultiresolutionProcessing.cpp
                             -------------------
    update               : 2004/09/18
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

#include "Frustum.h"
#include "MultiresolutionProcessing.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <valarray>
#include <cfloat>
#include <cmath>


//--
//
// MultiresolutionProcessing
//
//--
MultiresolutionProcessing::MultiresolutionProcessing(): mesh(0)
{
}


//--
//
// MultiresolutionProcessing
//
//--
MultiresolutionProcessing::MultiresolutionProcessing(MultiresolutionMesh* m)
{
	Model(m);
}


//--
//
// ~MultiresolutionProcessing
//
//--
MultiresolutionProcessing::~MultiresolutionProcessing()
{
}


//--
//
// Reset
//
//--
bool MultiresolutionProcessing::Model(MultiresolutionMesh* m)
{
	mesh = m;
	if( !mesh ) {
		original_vertices.clear();
		original_faces.clear();
		original_colors.clear();
		original_geometric_details.clear();
		original_color_details.clear();
		return false;
	}

	original_vertices = mesh->Vertices();
	original_faces = mesh->Faces();
	original_geometric_details.resize( mesh->LevelNumber() );
	for( int i=0; i<mesh->LevelNumber(); i++ ) {
		original_geometric_details[i] = mesh->Level(i).geometric_details;
	}

	if( mesh->ColorNumber() == mesh->VertexNumber() ) {
		original_colors = mesh->Colors();
		original_color_details.resize( mesh->LevelNumber() );
		for( int i=0; i<mesh->LevelNumber(); i++ ) {
			original_color_details[i] = mesh->Level(i).color_details;
		}
	}
	else {
		original_colors.clear();
		original_color_details.clear();
	}

//	Statistics();
	return true;
}


//--
//
// Reset
//
//--
bool MultiresolutionProcessing::Reset()
{
	// Check mesh
	if( !mesh ) return false;

	// Restore original geometric details
	for( int i=0; i<mesh->LevelNumber(); i++ ) {
		mesh->Level(i).geometric_details = original_geometric_details[i];
	}

	if( simplification_contractions.size() != 0 ) {
		// Insert odd vertices
		ConstPairContractionIterator it = simplification_contractions.begin();
		while( it != simplification_contractions.end() ) {
			// Expand candidate of current pair contraction
			mesh->Expand( *it );
			// Next pair contraction
			++it;
		}
		simplification_contractions.clear();
	}

	// Restore original vertices
	// Necessary to avoid weird behavior during reconstruction (yes, I need to check that !)
	mesh->Vertices() = original_vertices;

	// Restore original colors
	if( (int)original_colors.size() == mesh->VertexNumber() ) {
		mesh->Colors() = original_colors;
		for( int i=0; i<mesh->LevelNumber(); i++ ) {
			mesh->Level(i).color_details = original_color_details[i];
		}
	}
	else {
		mesh->ClearColors();
	}

	// Go to the base level
	mesh->GoToBaseLevel();

	// Rebuild the model
	mesh->Synthesis();

	// Update the normals
	mesh->ComputeFaceNormals();
	mesh->ComputeVertexNormals();

	return true;
}


//--
//
// Filter
//
//--
bool MultiresolutionProcessing::Filter(const std::vector<double>& fb)
{
	// Check mesh
	if( !mesh ) return false;

	// Go to the base level
	mesh->GoToBaseLevel();

	for( int l=0; l<3; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[5];
		}
	}

	for( int l=3; l<6; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[4];
		}
	}

	for( int l=6; l<9; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[3];
		}
	}

	for( int l=9; l<12; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[2];
		}
	}

	for( int l=12; l<16; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[1];
		}
	}

	for( int l=16; l<19; l++ ) {
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Shortcut to the detail of the current vertex
			level.geometric_details[i] = original_geometric_details[l][i] * fb[0];
		}
	}

	// Rebuild the model
	mesh->Synthesis();

	// Update the normals
	mesh->ComputeFaceNormals();
	mesh->ComputeVertexNormals();
	return true;
}


//--
//
// GeometricDetailNoiseVariance
//
//-
// Estimate the noise variance using the median of the detail coefficient lengths
//--
double MultiresolutionProcessing::GeometricDetailNoiseVariance()
{
	double min = DBL_MAX;
	double max = DBL_MIN;
	double temp;
	int total = 0;
	double noise_variance = 0;
	std::vector<int> histo(1000,0);
		
	// Find min/max details and count the number of details
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( mesh->Level(0).vertex_types[i] == UNDEFINED_VERTEX ) continue;
		temp = mesh->Level(0).geometric_details[i].Length();
		total++;
		if( temp > max ) max = temp;
		if( temp < min ) min = temp;
	}
	
	// Compute the histogram
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( mesh->Level(0).vertex_types[i] == UNDEFINED_VERTEX ) continue;
		temp = mesh->Level(0).geometric_details[i].Length();
		histo[(int)(temp/max*999.0)]++;
	}
	
	// Find the variance
	for( int i=0, count=0; i<1000; i++ ) {
		count += histo[i];
		if( count > (total/2) ) {
			noise_variance = (double)i/999.0 * max;
			break;
		}
	}

	return noise_variance;	

}


//--
//
// GeometricDetailAngleNoiseVariance
//
//-
// Estimate the variance of the noise on the angle between the detail coefficient and the vertex normal
//--
double MultiresolutionProcessing::GeometricDetailAngleNoiseVariance()
{
	double min = DBL_MAX;
	double max = DBL_MIN;
	double temp;
	int total = 0;
	double noise_variance = 0;
	std::vector<int> histo(1000,0);
	
	// Find min/max details and count the number of details
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( mesh->Level(0).vertex_types[i] == UNDEFINED_VERTEX ) continue;
		temp = fabs(mesh->Level(0).geometric_details[i] | mesh->VertexNormal(i));
		total++;
		if( temp > max ) max = temp;
		if( temp < min ) min = temp;
	}
	
	// Compute the histogram
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( mesh->Level(0).vertex_types[i] == UNDEFINED_VERTEX ) continue;
		temp = fabs(mesh->Level(0).geometric_details[i] | mesh->VertexNormal(i));
		histo[(int)(temp/max*999.0)]++;
	}
	
	// Find the variance
	for( int i=0, count=0; i<1000; i++ ) {
		count += histo[i];
		if( count > (total/2) ) {
			noise_variance = (double)i/999.0 * max;
			break;
		}
	}
	
	return noise_variance;
}

template<class Type>
Type sqr(const Type& x)
{
	return x*x;
}



//--
//
// GeometricDetailVarianceN1
//
//--
double MultiresolutionProcessing::GeometricDetailVarianceN1(int i, int l)
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->Level(l);
	
	// Local detail coefficient variance estimation
	double mean = sqr(level.geometric_details[i].Length());
	int count = 1;
	// Shortcut to the local neighborhood
	const NeighborList& neighborhood = mesh->NeighborVertices(i);
	// Vertex iterator in the neighborhood
	ConstNeighborIterator it;

	// Compute the mean
	for( it = neighborhood.begin(); it != neighborhood.end(); ++it ) {
		if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
		mean += sqr(level.geometric_details[*it].Length());
		count++;
	}
	mean /= (double)count; // Mean of the local detail coefficients
	return mean;
}	


//--
//
// GeometricDetailVarianceN2
//
//--
double MultiresolutionProcessing::GeometricDetailVarianceN2(int i, int l)
{

	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->Level(l);

	// neighbor vertex list
	NeighborList n2 = mesh->NeighborVertices(i);

	// Get the 2-ring neighborhood
	for( ConstNeighborIterator it = mesh->NeighborVertices(i).begin(); it != mesh->NeighborVertices(i).end(); ++it ) {
		n2.insert( n2.end(), mesh->NeighborVertices(*it).begin(), mesh->NeighborVertices(*it).end() );
	}
	
	// Remove redundant entries
	n2.sort();
	n2.unique();
	
	// Local detail coefficient variance estimation
	double mean = 0;
	for( ConstNeighborIterator it = n2.begin(); it != n2.end(); ++it ) {
		if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
		mean += level.geometric_details[*it].Length();
	}
	return mean /= (double)n2.size(); // Mean of the local detail coefficients
}	
		
//--
//
// GeometricDetailAngleVarianceN1
//
//--
double MultiresolutionProcessing::GeometricDetailAngleVarianceN1(int i, int l)
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->Level(l);
	
	// Local detail coefficient variance estimation
	double signal_variance = 0;
	double mean = 0;
	int count = 1;
	for( ConstNeighborIterator it = mesh->NeighborVertices(i).begin(); it != mesh->NeighborVertices(i).end(); ++it ) {
		if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
		mean += fabs(level.geometric_details[i] | level.geometric_details[*it]);
		count++;
	}
	mean /= (double)count; // Mean of the local detail coefficients
	count--; // Variance is normalize by N-1 (see Numerical Recipes)
	for( ConstNeighborIterator it = mesh->NeighborVertices(i).begin(); it != mesh->NeighborVertices(i).end(); ++it ) {
		if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
		signal_variance += sqr( fabs(level.geometric_details[i] | level.geometric_details[*it]) - mean );
	}
	return signal_variance /= (double)count;
}


//========================================================
//
// Get2ringNeighborVertices
//
//========================================================
NeighborList Get2ringNeighborVertices( NeighborMesh* mesh, int v )
{
	// Get 1-ring neighboring faces
	NeighborList neighbors = mesh->NeighborVertices(v);
	// Get 2-ring neighboring faces
	NeighborIterator it = mesh->NeighborVertices(v).begin();
	while( it != mesh->NeighborVertices(v).end() ) {
		neighbors.insert( neighbors.end(), mesh->NeighborVertices(*it).begin(), mesh->NeighborVertices(*it).end() );
		++it;
	}
	// Sort the neighbor list and remove redundant entries
	neighbors.sort();
	neighbors.unique();
	// Return the 2-ring neighbor list
	return neighbors;
}


//========================================================
//
// AverageEdgeLength
//
//========================================================
double AverageEdgeLength( const NeighborMesh* mesh )
{
	// Compute average edge length
	int edge_number = 0;
	double average_edge_length = 0.0;
	ConstNeighborIterator it;
	// Compute average edge length
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		// Get the neighborhood
		const NeighborList& neighborhood = mesh->NeighborVertices(i);
		for( it = neighborhood.begin(); it != neighborhood.end(); ++it ) {
			average_edge_length += (mesh->Vertex(i)-mesh->Vertex(*it)).Length();
			edge_number++;
		}
	}
	
	// Double the edge number because they are counted twice
	edge_number *= 2;
	
	// Return the average edge length
	return average_edge_length /= (double)edge_number;
}



//========================================================
//
// AverageOffsetStandardDeviation
//
//========================================================
double AverageOffsetStandardDeviation( NeighborMesh* mesh )
{
	// Compute average standard deviation of the offset in the neighborhood
	double standard_deviation = 0;
	ConstNeighborIterator it;
	
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
	
		// Get the neighborhood
	//	const NeighborList& neighborhood = mesh->NeighborVertices(i);
		NeighborList neighborhood = Get2ringNeighborVertices(mesh, i);

		// First pass to get the mean
		double mean_offset = 0;
		for( it = neighborhood.begin(); it != neighborhood.end(); ++it ) {
			mean_offset += mesh->VertexNormal(i) | (mesh->Vertex(i) - mesh->Vertex(*it));
		}
		mean_offset /= (double)mesh->NeighborVertexNumber(i);
	
		// Second pass to get the standard offset variance
		double offset_variance = 0;
		for( it = neighborhood.begin(); it != neighborhood.end(); ++it ) {
			double temp = (mesh->VertexNormal(i) | (mesh->Vertex(i) - mesh->Vertex(*it))) - mean_offset;
			offset_variance += temp*temp;
		}
		offset_variance /= (double)mesh->NeighborVertexNumber(i) + 1.0;
		standard_deviation += sqrt(offset_variance);
		
	}
	
	return standard_deviation /= (double)mesh->VertexNumber();
}


//--
//
// Denoise
//
//-
// Denoise the model using hard or soft thresholding of the detail coefficients
//--
bool MultiresolutionProcessing::Denoise( double threshold, int level_number, bool color )
{
	// Sanety checks
	if( mesh == 0 ) return false;
	if( threshold < 0 ) threshold = DBL_MAX;
	if( level_number < 0 ) level_number = 0;
	else if( level_number >= mesh->LevelNumber() ) level_number = mesh->LevelNumber()-1;

	// Go to the base level
//	mesh->GoToBaseLevel();
	int l = level_number;
	mesh->GoToLevel(l+1);
	
	double noise_variance = GeometricDetailNoiseVariance()/0.6745;
	std::cout<<noise_variance<<std::endl;
	noise_variance = threshold;

	
	// Apply the detail shrinkage denoising method
	for( int l=0; l<=level_number; l++ ) {
	
		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->Level(l);
		
		// For each vertex
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
		
			// Valid vertex ?
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;


			// Shortcut to the detail coefficients
			std::map<int,Vector3d>* details = 0;
			// Color or geometric details
			if( color ) details = &(level.color_details);
			else details = &(level.geometric_details);
			
			// Length of the detail
			double length = (*details)[i].Length();
			
			// Compute signal variance
//			double signal_variance = GeometricDetailVarianceN1( i, l );
			double signal_variance = (*details)[i].SquareLength();
			int count = 1;
			const NeighborList& neighborhood = mesh->NeighborVertices(i);
			for( ConstNeighborIterator it = neighborhood.begin(); it != neighborhood.end(); ++it ) {
				if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
				signal_variance += (*details)[*it].SquareLength();
				count++;
			}
			signal_variance /= (double)count + 1.0;

			// Remove noise variance from signal variance
			signal_variance = std::max(0.0, signal_variance - noise_variance*noise_variance);

			// Compute threshold according to BayesShrink rule
			if( signal_variance == 0 ) threshold = DBL_MAX;
			else threshold = noise_variance*noise_variance / sqr(signal_variance);

			// Detail length higher than threshold
			if( length > threshold ) {
				// Detail vector direction
				Vector3d direction = (*details)[i];
				direction.Normalize();
				// Soft thresholding
				(*details)[i] = (*details)[i] - direction * threshold;
			}
			// Detail length lower than threshold
			else {
				// Annihilate detail
				(*details)[i] = 0.0;
			}
		}
	}

	// Rebuild the model
	mesh->Synthesis();

	// Update the normals
	mesh->ComputeFaceNormals();
	mesh->ComputeVertexNormals();
	
	// Get a fresh, denoised model !
	return true;
}


//--
//
// LocalGeometricDetailMedian
//
//-
// Compute the median of the geometric detail coefficient lengths
//--
double MultiresolutionProcessing::LocalGeometricDetailMedian(int i, int l)
{
	static std::list<double> details;
	
	ResolutionLevel level = mesh->Level(l);
	
	// Store the details in a list
	details.clear();
	details.push_back( level.geometric_details[i].Length() );
	for( ConstNeighborIterator it = mesh->NeighborVertices(i).begin(); it != mesh->NeighborVertices(i).end(); ++it ) {
		if( level.vertex_types[*it] == UNDEFINED_VERTEX ) continue;
		details.push_back( level.geometric_details[*it].Length() );
	}
	
	// Sort the detail list
	details.sort();

	// Find the variance (the middle value of the list)
	int num = (int)((double)details.size()/2.0);
	int count = 0;
	for( std::list<double>::iterator it = details.begin(); it != details.end(); ++it ) {
		count++;
		if( count >= num ) return *it;
	}
	
	return 0;
}


//--
//
// Threshold
//
//--
bool MultiresolutionProcessing::Threshold( double threshold, int level_number, SegmentationType type, int iteration )
{
	// Check mesh
	if( !mesh ) return false;

	// If the mesh doesn't have any colors, add some
	if( mesh->ColorNumber() != mesh->VertexNumber() ) mesh->Colors().resize( mesh->VertexNumber() );

	// Check level number
	if( level_number < 0 ) level_number = 0;
	else if( level_number >= mesh->LevelNumber() ) level_number = mesh->LevelNumber()-1;

	// Go to the desired level number
	mesh->GoToLevel(level_number);

	std::valarray<bool> morpho1(mesh->VertexNumber());
	std::valarray<bool> morpho2(mesh->VertexNumber());

	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->Level(level_number);

	// Shortcut to the details
	std::map<int, Vector3d>* details = 0;

	switch( type ) {

		default :
		case DETAIL_GEOMETRIC :
			details = &(level.geometric_details);
			break;

		case DETAIL_NORMAL :
			details = &(level.normal_details);
			break;

		case DETAIL_COLOR :
			details = &(level.color_details);
			break;

		case DETAIL_TEXTURE :
			details = &(level.texture_details);
			return false;
	}



	// Thresholding
	for( int i=0; i<mesh->VertexNumber(); i++ ) {

		if( level.vertex_types[i] == UNDEFINED_VERTEX ) {
			morpho1[i] = false;
			continue;
		}
		
		//
		// Detail coefficient length thresholding
		//
		double length = (*details)[i].Length();
		if( length >= threshold ) morpho1[i] = true;
		else morpho1[i] = false;
//		if( GeometricDetailAngleVarianceN1(i, level_number) >= threshold ) morpho1[i] = true;
//		else morpho1[i] = false;
	}

	// Segmentation
	switch( type ) {

		//
		// No segmentation
		//
		default :
		case SEGMENTATION_NONE :
			morpho2 = morpho1;
			break;

		//
		// Dilation
		//
		case SEGMENTATION_DILATION :
			for( int it=0; it<iteration; it++ ) {
				morpho2 = false;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == false ) continue;
					morpho2[i] = true;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = true;
						++it;
					}
				}
			}
			break;

		//
		// Erosion
		//
		case SEGMENTATION_EROSION :
			for( int it=0; it<iteration; it++ ) {
				morpho2 = true;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == true ) continue;
					morpho2[i] = false;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = false;
						++it;
					}
				}
			}
			break;

		//
		// Opening
		//
		case SEGMENTATION_OPENING :
			// Dilation
			for( int it=0; it<iteration; it++ ) {
				morpho2 = false;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == false ) continue;
					morpho2[i] = true;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = true;
						++it;
					}
				}
				morpho1 = morpho2;
			}
			// Erosion
			for( int it=0; it<iteration; it++ ) {
				morpho2 = true;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == true ) continue;
					morpho2[i] = false;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = false;
						++it;
					}
				}
			}
			break;

		//
		// Closing
		//
		case SEGMENTATION_CLOSING :
			// Erosion
			for( int it=0; it<iteration; it++ ) {
				morpho2 = true;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == true ) continue;
					morpho2[i] = false;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = false;
						++it;
					}
				}
				morpho1 = morpho2;
			}
			// Dilation
			for( int it=0; it<iteration; it++ ) {
				morpho2 = false;
				for( int i=0; i<mesh->VertexNumber(); i++ ) {
					if( morpho1[i] == false ) continue;
					morpho2[i] = true;
					ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
					while( it != mesh->NeighborVertices(i).end() ) {
						morpho2[*it] = true;
						++it;
					}
				}
			}
			break;
	}

	// Coloring
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( morpho2[i] == false ) mesh->Color(i) = Vector3d(0.0,0.5,1.0);
		else mesh->Color(i) = Vector3d(1.0,0.5,0.0);
	}

	// Update the normals
	mesh->ComputeFaceNormals();
	mesh->ComputeVertexNormals();

	return true;
}


//--
//
// Color
//
//--
// Color mesh according to detail lengths
bool MultiresolutionProcessing::Color( DetailType type, int level_number, bool smooth )
{
	if( !mesh ) return false;

	// If the mesh doesn't have any colors, add some
	if( mesh->ColorNumber() != mesh->VertexNumber() ) mesh->Colors().resize( mesh->VertexNumber() );

	// Check level number
	if( level_number < 0 ) level_number = 0;
	else if( level_number >= mesh->LevelNumber() ) level_number = mesh->LevelNumber() - 1;

	// Go to the desired level number
	mesh->GoToLevel(level_number);

	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->Level(level_number);

	// Shortcut to the detail array
	std::map<int, Vector3d>* details = 0;

	switch( type ) {

		default :
		case DETAIL_GEOMETRIC :
			details = &(level.geometric_details);
			break;

		case DETAIL_NORMAL :
			details = &(level.normal_details);
			break;

		case DETAIL_COLOR :
			details = &(level.color_details);
			break;

		case DETAIL_TEXTURE :
			details = &(level.texture_details);
			break;
	}


	// Get min/max detail length
	double min = DBL_MAX;
	double max = DBL_MIN;
	double temp;
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
		// Geometic detail length
		temp = (*details)[i].Length();
//		temp = GeometricDetailVarianceN1(i, level_number);	
//		temp = GeometricDetailVarianceN2(i, level_number);
		if( temp > max ) max = temp;
		if( temp < min ) min = temp;
	}

	// Colorize the vertex according to the detail length
	max -= min;
	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( mesh->IsProgressiveVertexValid(i) == false ) continue;

		if( level.vertex_types[i] != UNDEFINED_VERTEX ) {
			temp = (*details)[i].Length();
//			temp = GeometricDetailVarianceN1(i, level_number);
//			temp = GeometricDetailVarianceN2(i, level_number);
			mesh->Color(i) = Double2Color( (temp-min)/max );
		}
		else {
			mesh->Color(i) = Vector3d(0,0,0);
		}
	}

	// Apply a median filter to smooth the resulting coloring
	if( smooth == true ) {
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			// Valid vertex ?
			if( mesh->IsProgressiveVertexValid(i) == false ) continue;
			// Valid detail
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			// Replace the color of the current vertex by
			// the average of the colors in the neighborhood
			Vector3d average = mesh->Color(i);
			int count = 1;
			ConstNeighborIterator itv = mesh->NeighborVertices(i).begin();
			while( itv != mesh->NeighborVertices(i).end() ) {
				if( level.vertex_types[*itv] != UNDEFINED_VERTEX ) {
					average += mesh->Color(*itv);
					count++;
				}
				++itv;
			}
			mesh->Color(i) = average / (double)count;
		}
	}

	// Update the normals
	mesh->ComputeFaceNormals();
	mesh->ComputeVertexNormals();

	return true;
}


//--
//
// Simplify
//
//--
// Adaptative simplification
bool MultiresolutionProcessing::Simplify(double threshold, int level_number, DetailType type, bool view_dependent, Viewport vp)
{
	// Check mesh
	if( !mesh ) return false;

	// Check level number
	if( level_number < 0 ) level_number = 0;
	else if( level_number >= mesh->LevelNumber() ) level_number = mesh->LevelNumber() - 1;

	// Check treshold
	if( threshold<0 ) threshold = DBL_MAX;

	// Go the initial level
	if( simplification_contractions.size() != 0 ) {
		// Insert odd vertices
		ConstPairContractionIterator it = simplification_contractions.begin();
		while( it != simplification_contractions.end() ) {
			// Expand candidate of current pair contraction
			mesh->Expand( *it );
			// Next pair contraction
			++it;
		}
		simplification_contractions.clear();
	}
	else {
		// Go to the initial level
		mesh->GoToInitialLevel();
	}

	Vector3d dir(0,0,-1);

	if( view_dependent ) {
		// Viewing vector
		dir = vp.z_axis;
	}

	// Field of vision angle / 2
	double fov_2 = vp.fov / 2.0;

	Frustum frustum;
	frustum.CalculateFrustum();

	// Simplify
	for( int l=0; l<level_number; l++ ) {

		// Get level
		ResolutionLevel& level = mesh->Level( l );

		// Shortcut to the details of this level
		std::map<int, Vector3d>* details = 0;
		switch( type ) {
			default :
			case DETAIL_GEOMETRIC : details = &(level.geometric_details); break;
			case DETAIL_NORMAL : details = &(level.normal_details); break;
			case DETAIL_COLOR : details = &(level.color_details); break;
			case DETAIL_TEXTURE : details = &(level.texture_details); break;
		}

		// Go through all the pair contraction of this level
		for( ConstPairContractionIterator itp=level.pair_contractions.begin(); itp!=level.pair_contractions.end(); ++itp ) {

			// Valid odd vertex ?
			if( level.vertex_types[itp->Candidate()] != ODD_VERTEX ) continue;

			// Valid pair contraction ?
			if( !mesh->IsValidPairContraction(*itp) ) continue;

			// Visibility test
			if( view_dependent &&
			// Back face culling
			((mesh->VertexNormal(itp->Candidate())|dir)>-0.1/*sin(fov_2)*/ ||
			// Frustum clipping
			!frustum.PointInFrustum(mesh->Vertex(itp->Candidate())))
			) {
				// Create a new pair contraction
				PairContraction pair(*itp);
				// Initialize the common faces
				pair.ClearCommonFaces();
				// Record faces on edge va-vb
				ConstNeighborIterator it = mesh->NeighborFaces(pair.Candidate()).begin();
				while( it != mesh->NeighborFaces(pair.Candidate()).end() ) {
					if( mesh->FaceHasVertex( *it, pair.Target() ) ) {
						pair.AddCommonFace(*it);
					}
					++it;
				}
				// Collapse
				mesh->Collapse( pair );
				mesh->UpdateNormals(pair);
				// Store the contraction
				simplification_contractions.push_front(pair);
			}

			// Visible point
			// Detail relevance test
			else if( /*((mesh->VertexNormal(itp->Candidate())|dir)<-0.2) &&*/ ((*details)[itp->Candidate()].Length()<threshold) ) {
				// Create a new pair contraction
				PairContraction pair(*itp);
				// Initialize the common faces
				pair.ClearCommonFaces();
				// Record faces on edge va-vb
				ConstNeighborIterator it = mesh->NeighborFaces(pair.Candidate()).begin();
				while( it != mesh->NeighborFaces(pair.Candidate()).end() ) {
					if( mesh->FaceHasVertex( *it, pair.Target() ) ) {
						pair.AddCommonFace(*it);
					}
					++it;
				}
				// Collapse
				mesh->Collapse( pair );
				mesh->UpdateNormals(pair);
				// Store the contraction
				simplification_contractions.push_front(pair);
			}
		}
	}

	return true;
}


//--
//
// Statistics
//
//--
// Compute and write in files statistics about the multiresolution model
// essentially about the details
bool MultiresolutionProcessing::Statistics()
{
	// Check mesh
	if( !mesh ) return false;

	std::ofstream file_levels("levels.log");
	using std::endl;
	
	// Go the finest level
	mesh->GoToInitialLevel();

	// Stat
	for( int l=0; l<mesh->LevelNumber(); l++ ) {

		// Temp variables
		double min = DBL_MAX;
		double max = DBL_MIN;
		double mean = 0;
		double standard_deviation = 0, absolute_deviation = 0;
		double temp;
		int count = 0;
		std::vector<int> histo(1000,0);
		std::vector<double> detail_vector;

		// Get level
		ResolutionLevel& level = mesh->Level( l );
	
		// Compute min, max, mean and histogram
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;
			temp = level.geometric_details[i].Length();
			mean += temp;
			detail_vector.push_back(temp);
			count++;
			if( temp > max ) max = temp;
			if( temp < min ) min = temp;
		}
		mean /= count;

		// Compute histogram and deviations
		for( int i=0; i<(int)detail_vector.size(); i++ ) {
			histo[(int)(detail_vector[i]/max*999.0)]++;
			standard_deviation += sqr(detail_vector[i] - mean);
			absolute_deviation += fabs(detail_vector[i] - mean);
		}
		standard_deviation /= count-1;
		standard_deviation = sqrt( standard_deviation );
		absolute_deviation /= count;
		
		// Compute median
		sort( detail_vector.begin(), detail_vector.end() );
		double median = detail_vector[ detail_vector.size()/2 ];
	
		// Store histogram
		char filename[100];
		sprintf(filename, "histogram-%02d.log", l);
		std::ofstream file_histo(filename);
		for( int i=0; i<1000; i++ ) {
			file_histo<<(double)i/999.0*max<<" "<<histo[i]<<endl;
		}
		file_histo.close();
	
		// Print level statistics
		file_levels<<l<<'\t'<<mesh->ValidVertexNumber()<<'\t'<<mesh->ValidFaceNumber()<<'\t'<<min<<'\t'<<max<<'\t'<<mean<<'\t'<<median<<'\t'<<standard_deviation<<'\t'<<absolute_deviation<<endl;

		// Go to next coarse level
		if( l < mesh->LevelNumber() ) mesh->DecomposeCurrentLevel();
	}

	// Close level stats file
	file_levels.close();
	
	// Go back to finest level
	mesh->GoToInitialLevel();
	
	// Finished
	return true;
}



//--
//
// Test
//
//--
void MultiresolutionProcessing::Test()
{
	// Check mesh
	if( !mesh ) return;

	static bool first_time = true;
	static std::valarray<bool> noise(mesh->VertexNumber());

	if( first_time ) {

		// If the mesh doesn't have any colors, add some
		if( mesh->ColorNumber() != mesh->VertexNumber() ) mesh->Colors().resize( mesh->VertexNumber() );

		mesh->GoToInitialLevel();

		// Shortcut to the current resolution level
		ResolutionLevel& level = mesh->CurrentLevel();

		// Thresholding
		for( int i=0; i<mesh->VertexNumber(); i++ ) {

			if( level.vertex_types[i] == UNDEFINED_VERTEX ) {
				noise[i] = false;
				continue;
			}
		
			//
			// Detail coefficient length thresholding
			//
			if( level.geometric_details[i].Length() >= 0.2 ) {
				int count=0;
				ConstNeighborIterator itv = mesh->NeighborVertices(i).begin();
				while( itv != mesh->NeighborVertices(i).end() ) {
					if( level.normal_details[*itv].Length() >= 0.01 ) {
						count++;
					}
					++itv;
				}
				if( count > (mesh->NeighborVertexNumber(i)/2) ) {
					noise[i] = true;
				}
			}
			else {
				noise[i] = false;
			}
		}

		// Coloring
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( noise[i] == false ) mesh->Color(i) = Vector3d(0.7,0.7,0.7);
			else mesh->Color(i) = Vector3d(1.0,0.0,0.0);
		}
		
		first_time = false;
	}
	else {
		// Denoising
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( noise[i] == true )  {
				mesh->Level(0).geometric_details[i]=0.0;
				mesh->Color(i) = Vector3d(0.7,0.7,0.7);
			}
		}
		mesh->GoToLevel(2);
		mesh->Synthesis();
		// Update the normals
		mesh->ComputeFaceNormals();
		mesh->ComputeVertexNormals();
	}

}




