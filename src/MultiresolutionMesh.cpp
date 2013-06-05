/***************************************************************************
                           MultiresolutionMesh.cpp
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

#include "MultiresolutionMesh.h"
#include "Curvature.h"

#include <algorithm>
#include <float.h>
#include <fstream>
#include <iomanip>
#include <valarray>


//#include "DebugLog.h"


static double detail_length_tolerance = 1e15;


//#define _TEST_PLAN_
//#define _OUTPUT_LEVELS_

//#ifdef _HSI_COLOR_

inline void Rgb2Hsi(const Vector3d& rgb_color, Vector3d& hsi_color)
{
	// Shortcuts
	const double& r = rgb_color[0];
	const double& g = rgb_color[1];
	const double& b = rgb_color[2];
	double& h = hsi_color[0];
	double& s = hsi_color[1];
	double& i = hsi_color[2];

	assert( r>=0.0 && r<=1.0 );
	assert( g>=0.0 && g<=1.0 );
	assert( b>=0.0 && b<=1.0 );

	// Intensity
	i = (r+g+b)/3.0;

	// Saturation
	if(i == 0) s = 1;
	else if((r == g) && (g == b)) s = 0;
	else s = 1.0 - 3.0*std::min(r, std::min(g, b))/(r+g+b);

	// Hue
	if((s==0) || (i==0)) h = 0;
	else {
		h = acos( 0.5 * ((r-g)+(r-b)) / sqrt((r-g)*(r-g)+(r-b)*(g-b)) );
		if(b/i > g/i) h = 2.0*M_PI - h;
	}
	// Normalize hue
	h /= 2.0*M_PI;
}

inline Vector3d Rgb2Hsi(const Vector3d& rgb_color)
{
	Vector3d hsi_color;
	Rgb2Hsi(rgb_color, hsi_color);
	return hsi_color;
}


inline void Hsi2Rgb(const Vector3d& hsi_color, Vector3d& rgb_color)
{
	// Shortcuts
	double h = hsi_color[0];
	const double& s = hsi_color[1];
	const double& i = hsi_color[2];
	double& r = rgb_color[0];
	double& g = rgb_color[1];
	double& b = rgb_color[2];

	assert( h>=0.0 && h<=1.0 );
	assert( s>=0.0 && s<=1.0 );
	assert( i>=0.0 && i<=1.0 );

	static const double TWO_PI_3 = 2.0 * M_PI / 3.0;
	static const double FOUR_PI_3 = 4.0 * M_PI / 3.0;
	static const double PI_3 = M_PI / 3.0;

	h *= 2.0 * M_PI;

	if( s == 0.0 ) {
		r = g = b = i;
	}
	else if( h < TWO_PI_3 ) {
		b = i * (1.0 - s);
		r = i * (1.0 + s * cos(h) / cos(PI_3 - h));
		g = 3.0 * i - (r + b);
	}
	else if( h < FOUR_PI_3 ) {
		h -= TWO_PI_3;
		r = i * (1.0 - s);
		g = i * (1.0 + s * cos(h) / cos(PI_3 - h));
		b = 3.0 * i - (r + g);
	}
	else {
		h -= FOUR_PI_3;
		g = i * (1.0 - s);
		b = i * (1.0 + s * cos(h) / cos(PI_3 - h));
		r = 3.0 * i - (g + b);
	}
}

inline double HsiDifference( const Vector3d& c1, const Vector3d& c2 )
{
	// Shortcuts
	double h1 = c1[0];
	const double& s1 = c1[1];
	const double& i1 = c1[2];
	double h2 = c2[0];
	const double& s2 = c2[1];
	const double& i2 = c2[2];

	assert( h1>=0.0 && h1<=1.0 );
	assert( s1>=0.0 && s1<=1.0 );
	assert( i1>=0.0 && i1<=1.0 );
	assert( h2>=0.0 && h2<=1.0 );
	assert( s2>=0.0 && s2<=1.0 );
	assert( i2>=0.0 && i2<=1.0 );

	h1 *= 2.0*M_PI;
	h2 *= 2.0*M_PI;

	if( c1 == c2 ) return 0;

	double delta_i = fabs(i1 - i2);
	double theta;
	double delta_s;
	if( i1==0  || i2==0 ) {
		theta = 0;
		delta_s = 0;
	}
	else {
		theta = fabs(h1 - h2);
		if( theta > M_PI ) theta = 2.0*M_PI - theta;
		delta_s = sqrt( s1*s1 + s2*s2 - 2.0*s1*s2*cos(theta) );
	}
	return sqrt( delta_i*delta_i + delta_s*delta_s );
}

//#endif


//-- Test --
//#include "Stopwatch.h"
//std::ofstream test("test1.log", std::ios::app); using std::endl; using std::flush;
//std::ofstream test("levelsize.log", std::ios::app); using std::endl; using std::flush;
//Stopwatch timer;
//int memory(0), tmp_mem;
//--



//#define _RGB_2_GRAYSCALE_



//--
//
// Analysis
//
//--
bool MultiresolutionMesh::Analysis(int)
{
	//
	// HSI
	//
	#ifdef _HSI_COLOR_
		std::vector<Vector3d> color_backup;
		if( ColorNumber() == VertexNumber() )
		{
			color_backup = colors;
			Vector3d tmp_hsi;
			for( int i=0; i<ColorNumber(); i++ )
			{
				Rgb2Hsi( Color(i), tmp_hsi );
				Color(i) = tmp_hsi;
			}
		}
	#endif


	//
	// RGB 2 Grayscale
	//
	#ifdef _RGB_2_GRAYSCALE_
		std::vector<Vector3d> color_backup;
		if( ColorNumber() == VertexNumber() ) {
			color_backup = colors;
			double tmp;
			for( int i=0; i<ColorNumber(); i++ ) {
				tmp = (Color(i)[0]+Color(i)[1]+Color(i)[2])/3.0;
				Color(i) = tmp;
			}
		}
	#endif
	

	// Initialize member variables
	current_level_number = -1;
	levels.clear();
	mean_curvature.clear();
	gaussian_curvature.clear();
	texture.Reset();
	vertex_map.clear();
	edge_list.clear();

	// Initialize normals
	ComputeFaceNormals();
	ComputeVertexNormals();

	// Initialize progressive decomposition
	BeginProgressiveDecomposition();
	

	// Initialize curvatures
//	mean_curvature.resize( VertexNumber() );
//	gaussian_curvature.resize( VertexNumber() );
//	ComputeMeanCurvature();
//	ComputeGaussianCurvature();

	// Initialize texture
	if( (TextureNumber()!=0) && (TextureName().empty()==false) ) {
		use_texture = texture.ReadFile( TextureName() );
	}
	else {
		use_texture = false;
	}

	if( ColorNumber() == VertexNumber() ) {
		use_color = true;
	}
	else {
		use_color = false;
	}

	use_normal = true;

//	permutation.assign(VertexNumber(), -1);
//	map.assign( VertexNumber(), -1);

	// Initialize vertex map
	// This table indicates on which vertex each vertex has collapsed to
	vertex_map.resize( VertexNumber() );
	for( int i=0; i<VertexNumber(); i++ ) {
		vertex_map[i] = i;
	}

//	std::cout<<bounding_box.Center()<<std::endl;
//	for( int i=0; i<VertexNumber(); i++ ) {
//		if( 
//			Vertex(i)[0]>-0.01 &&
//			Vertex(i)[2]>-0.01 &&
//			Vertex(i)[0]<0.01 &&
//			Vertex(i)[2]<0.01
//		)
//			std::cout<<i<<" "<<Vertex(i)<<std::endl;
//	}


	ConstPairContractionIterator itp;
	ConstNeighborIterator itn;
	std::vector<bool> locked_vertices(VertexNumber(), false);
//	levels.resize( base_level );
	current_level_number = 0;

	//-- Test --
	//int vnum = VertexNumber();
	//timer.Reset();
	//timer.Start();
	//test<<current_level_number<<'\t'<<vnum<<endl;
	//--

	#ifdef _OUTPUT_LEVELS_
		// Write the finest level statistics
		std::ofstream level_stats_file("levels.log");
		level_stats_file<<current_level_number<<'\t'<<ValidVertexNumber()<<'\t'<<ValidFaceNumber()<<std::endl;
		// Write the finest level mesh (the initial mesh)
		char filename[255];
		sprintf(filename,"level%02d.wrl",current_level_number);
		WriteFile(filename);
	#endif


	//
	// Create the levels of details
	//
	while( 1 ) {

		// Add a new level
		levels.push_back( ResolutionLevel() );

		//
		// Select the odd vertices
		// Remove the odd vertices
		// Lock the neighbors
		// Odd vertices are a set of independent vertices selected to be removed
		// The remaining (even) vertices compose the next coarse level
		//
		while( !pair_contractions.empty() ) {

			#ifdef _TEST_PLAN_
				// Lock one predefined vertex
				locked_vertices[1844] = true; // plan-20000
			//	locked_vertices[30670] = true; // plan-irregulier
			#endif

			// Get pair contraction with minimum cost
			PairContraction pair = MinimumPairContraction();
			
			// Remove pair from contraction candidate list
			RemovePairContraction( pair.Candidate() );
			
			// Do not perform bad contraction until the level of detail is very low
			if( pair.Cost()>=invalid_contraction_penalty && current_level_number<15 ) {
				locked_vertices[pair.Candidate()] = true;
				continue;
			}
			
			// Test to see if after the pair contraction,
			// one vertex has a valence of 3
			// This case screws up the predict operator
			itn = NeighborVertices(pair.Candidate()).begin();
			while( itn != NeighborVertices(pair.Candidate()).end() ) {
				if(current_level_number<15 && NeighborVertexNumber(*itn)==4 && FindNeighborVertex(pair.Target(), *itn) && !IsBorderVertex(*itn)) {
					locked_vertices[pair.Candidate()] = true;
				}
				++itn;
			}
			
			if( locked_vertices[pair.Candidate()] ) continue;
			
			// Map the candidate to the target
			vertex_map[pair.Candidate()] = pair.Target();
			
			// Lock the selected vertex and its neighbors
		//	valid_vertices[pair.Candidate()] = false;
		
			// Collapse candidate
			Collapse( pair );
			
			// Update normals
			UpdateNormals( pair );
			
			// Update quadrics
			UpdateQuadrics(pair);
			
			// Update pair cost of surronding vertices
			if( pair.Target() >= 0 ) {	
				// Compute neighbor vertex cost
				itn = NeighborVertices(pair.Target()).begin();
				while( itn != NeighborVertices(pair.Target()).end() ) {
					if( locked_vertices[*itn] == false ) {
						RemovePairContraction(*itn);
						ComputeEdgeCostAtVertex(*itn);
					}
					++itn;
				}
			}
			
			// Lock the neighborhood of the selected vertex
			itn = NeighborVertices(pair.Candidate()).begin();
			while( itn != NeighborVertices(pair.Candidate()).end() ) {
				// If the vertex is not already locked
				if( locked_vertices[*itn] == false ) {
					// Remove its pair contraction estimation
					RemovePairContraction(*itn);
					// Lock it
					locked_vertices[*itn] = true;
				}
				// Next neighbor
				++itn;
			}
			
			// Record pair contraction
			levels[current_level_number].pair_contractions.push_front( pair );
		}

//		bool left = false;
		int left_vertex_number = 0;

		// Compute new edge collapse pairs
		for(int i=0; i<VertexNumber(); i++ ) {
			if( valid_vertices[i] == false ) continue;
			left_vertex_number++;
			locked_vertices[i] = false;
			ComputeEdgeCostAtVertex(i);
		}

		// Update the normals
		ComputeProgressiveFaceNormals();
		ComputeProgressiveVertexNormals();

		// Next level
		current_level_number++;
		
		#ifdef _OUTPUT_LEVELS_
			// Write the current level statistics
			level_stats_file<<current_level_number<<'\t'<<ValidVertexNumber()<<'\t'<<ValidFaceNumber()<<std::endl;
			// Write the current level mesh
			sprintf(filename,"level%02d.wrl",current_level_number);
			WriteFile(filename);
		#endif

		// Is there any vertex left ?
//		if( left_vertex_number == 0 ) break;
	//	if( left_vertex_number < 100 ) break;
		if( left_vertex_number < 10 ) break;

	}

	// End the progressive decomposition
	EndProgressiveDecomposition();

	// Reset
	locked_vertices.clear();


	//
	// Go back the initial fine level
	// Reconstruct all levels
	//
	while( current_level_number > 0 ) {
		// Next level
		current_level_number--;

		// Insert odd vertices
		itp = levels[current_level_number].pair_contractions.begin();
		while( itp != levels[current_level_number].pair_contractions.end() ) {
			// Expand candidate of current pair contraction
			Expand( *itp );
			// Next pair contraction
			++itp;
		}
	}

	// Update the normals
	ComputeProgressiveFaceNormals();
	ComputeProgressiveVertexNormals();

	for( int i=0; i<(int)levels.size(); i++ ) {
		levels[i].Resize( VertexNumber() );
	}

	//
	// Compute the details
	//
	while( current_level_number < LevelNumber() ) {

		//
		// Predict odd vertices
		//
		itp = levels[current_level_number].pair_contractions.begin();
		while( itp != levels[current_level_number].pair_contractions.end() ) {
			// Compute wavelet coefficient
			PredictVertex( itp->Candidate(), ODD_VERTEX );
			
			// Invalidate the odd vertex
			valid_vertices[itp->Candidate()] = false;

			#ifdef _TEST_PLAN_
				// Null the details
				levels[current_level_number].geometric_details[itp->Candidate()] = 0;
			#endif
			
			// Next pair contraction
			++itp;
		}


		//
		// Predict even vertices
		//
		for(int i=0; i<VertexNumber(); i++ ) {
			// Valid vertex ?
			if( valid_vertices[i] == false ) continue;
			
			// Predict even vertex
			PredictVertex( i, EVEN_VERTEX );

			#ifdef _TEST_PLAN_
				// Null the details
				levels[current_level_number].geometric_details[i] = 0;
			#endif
		}

		//
		// Simplify
		// Removed odd vertices to go to next coarse level
		//
		itp = levels[current_level_number].pair_contractions.begin();
		while( itp != levels[current_level_number].pair_contractions.end() ) {
		
			// Collapse candidate
			Collapse( *itp );
			
			// Next Pair contrction
			++itp;
		}

		// Update the normals
		ComputeProgressiveFaceNormals();
		ComputeProgressiveVertexNormals();

		// Go to the next coarse level
		current_level_number++;
	}

	// Test plan
	#ifdef _TEST_PLAN_
		Vertex(1844)[1] = -1.0; // plan-2000
//		Vertex(30670)[1] = -1.0; // plan-irregulier
		char filename[255];
		sprintf(filename,"plan-base.wrl");
		WriteFile(filename);
	#endif


	//-- Test --
	//timer.Stop();
	//test<<VertexNumber()<<'\t'<<vnum<<'\t'<<((double)VertexNumber()/vnum - 1.0)*100.0<<'\t'<<timer.Total()<<'\t'<<memory<<endl;
	//--

	//
	// HSI
	//
/*
	#ifdef _HSI_COLOR_
	if( ColorNumber() == VertexNumber() )
	{
		colors = color_backup;
	}
	#endif
*/

	return true;
}


//--
//
// Synthesis
//
//--
bool MultiresolutionMesh::Synthesis()
{
	return GoToInitialLevel();
}


//--
//
// DecomposeCurrentLevel
//
//--
bool MultiresolutionMesh::DecomposeCurrentLevel()
{
	// Test if we are not already at the coarsest level
	if( current_level_number >= (int)levels.size() ) return false;

	// Removed odd vertices of the current level
	std::list<PairContraction>::iterator it = levels[current_level_number].pair_contractions.begin();
	while( it != levels[current_level_number].pair_contractions.end() ) {
		// Collapse candidate
		Collapse( *it );
		// Update normals
		UpdateNormals( *it );
		// Next Pair contrction
		++it;
	}

	// Go to the next coarse level
	current_level_number++;

	return true;
}


//--
//
// ReconstructCurrentLevel
//
//--
bool MultiresolutionMesh::ReconstructCurrentLevel()
{
	// Test if we are not already at the finest level
	if( current_level_number <= 0 ) return false;

	// Recontruct the next fine level
	current_level_number--;

	// Insert odd vertices
	std::list<PairContraction>::iterator it = levels[current_level_number].pair_contractions.begin();
	while( it != levels[current_level_number].pair_contractions.end() )
	{
		// Expand candidate of the current pair contraction
		Expand( *it );
		// Next pair contraction
		++it;
	}

	// Reconstruct odd vertices
	for( int i=0; i<VertexNumber(); i++ ) {
		if( levels[current_level_number].vertex_types[i] != ODD_VERTEX ) continue;
		Vertex(i) = ReconstructVertex(i);
	}
/*
	// Reconstruct even vertices
	std::vector<Vector3d> temp_vertices = vertices;
	for( int i=0; i<VertexNumber(); i++ ) {
		if( levels[current_level_number].vertex_types[i] != EVEN_VERTEX ) continue;
		temp_vertices[i] = ReconstructVertex(i);
	}
	vertices = temp_vertices;
*/	
	if( use_color ) {
		// Reconstruct odd vertices
		for( int i=0; i<VertexNumber(); i++ ) {
			if( levels[current_level_number].vertex_types[i] != ODD_VERTEX ) continue;
			Color(i) = ReconstructVertexColor(i);
		}
/*
		// Reconstruct even vertices
		temp_vertices = colors;
		for( int i=0; i<VertexNumber(); i++ ) {
			if( levels[current_level_number].vertex_types[i] != EVEN_VERTEX ) continue;
			temp_vertices[i] = ReconstructVertexColor(i);
		}
		colors = temp_vertices;
*/
	}

	return true;
}


//--
//
// GoToLevel
//
//--
bool MultiresolutionMesh::GoToLevel(int i)
{
	// Check level number
	if( i < 0 ) i = 0;
	else if( i > LevelNumber() ) i = LevelNumber();

	// Go to the desired level number
	if( i < CurrentLevelNumber() ) {
		while( i < CurrentLevelNumber() ) {
			ReconstructCurrentLevel();
		}
	}
	else if( i > CurrentLevelNumber() ) {
		while( i > CurrentLevelNumber() ) {
			DecomposeCurrentLevel();
		}
	}
	return true;
}


//--
//
// GoToBaseLevel
//
//--
bool MultiresolutionMesh::GoToBaseLevel()
{
	while( CurrentLevelNumber() < LevelNumber() ) {
		DecomposeCurrentLevel();
	}
	return true;
}


//--
//
// GoToInitialLevel
//
//--
bool MultiresolutionMesh::GoToInitialLevel()
{
	while( CurrentLevelNumber() > 0 ) {
		ReconstructCurrentLevel();
	}
	return true;
}


//--
//
// PredictVertex
//
//--
bool MultiresolutionMesh::PredictVertex(int v, VertexType vt)
{
	double coef;
	double weight(0);
	ConstNeighborIterator itv1, itv2;

	// Initialise detail validity
	levels[current_level_number].vertex_types[v] = UNDEFINED_VERTEX;

	//-- Test --
//	tmp_mem = 0;
	//--
	
	// Check border vertex
	if( IsBorderVertex(v) ) return false;

	// For every neighbor vertices
	itv1 = NeighborVertices(v).begin();
	while( itv1 != NeighborVertices(v).end() ) {
		// Check for topological noise
	//	if( NeighborVertexNumber(*itv1) < 4 ) return false;

	
		coef = 0;
		itv2 = NeighborVertices(*itv1).begin();
		while( itv2 != NeighborVertices(*itv1).end() ) {
			if( FindNeighborVertex(v, *itv2) ) {
				if( IsColinear( Vertex(*itv2), Vertex(v), Vertex(*itv1) ) ) return false;
				coef += Cotan( Vertex(*itv2), Vertex(v), Vertex(*itv1) );
			}
			++itv2;
		}
	
		// Laplacian relaxation
	//	coef = 1.0 / (double)NeighborVertexNumber(v);
	
		// Save the current coeff
		levels[current_level_number].subdivision_weights[v].push_back( SubdivisionWeight(*itv1,coef) );

		//-- Test --
		//tmp_mem++;
		//

		// Add current coef to normalization coef
		weight += coef;
		++itv1;
	}
	
	// Sanity check
	if( weight == 0 ) {
		return false;
	}

	// Normalize weights
	weight = 1.0 / weight;
	SubdivisionWeightList::iterator itw = levels[current_level_number].subdivision_weights[v].begin();
	while( itw != levels[current_level_number].subdivision_weights[v].end() ) {
		itw->weight *= weight;
		++itw;
	}

	// Compute predicted vertex position
	Vector3d relax(0,0,0);
	itw = levels[current_level_number].subdivision_weights[v].begin();
	while( itw != levels[current_level_number].subdivision_weights[v].end() ) {
		relax += itw->weight * Vertex(itw->vertex);
		++itw;
	}

	// Compute detail (vertex - subdivided_vertex)
	levels[current_level_number].geometric_details[v] = Vertex(v) - relax;

	//-- Test --
	//tmp_mem++;
	//

	// Sanety check (degenerated detail)
//	if( levels[current_level_number].geometric_details[v].Length() > detail_length_tolerance ) return false;

	// Validate detail
	levels[current_level_number].vertex_types[v] = vt;

	// Create predicted mesh fo geometrical attribute details
	std::swap( Vertex(v), relax );

	// Normal vector detail
	Vector3d relax_normal(0,0,0);
	itv1 = NeighborFaces(v).begin();
	while( itv1 != NeighborFaces(v).end() )
	{
		relax_normal += ComputeFaceNormal( *itv1 );
		++itv1;
	}
	relax_normal.Normalize();
//	levels[current_level_number].normal_details[v] = VertexNormal(v) - relax_normal;
	levels[current_level_number].normal_details[v] = Vector3d(acos(VertexNormal(v) | relax_normal), 0, 0);

	// Gaussian curvature
//	relax_curvature = GaussianCurvature(this, v);
//	levels[current_level_number].gaussian_curvature_details[v] = gaussian_curvature[v] - relax_curvature;

	// Go back to the initial mesh
	std::swap( Vertex(v), relax );

	// Color detail
	if( use_color )
	{
		relax = 0.0;
		itw = levels[current_level_number].subdivision_weights[v].begin();
		while( itw != levels[current_level_number].subdivision_weights[v].end() ) {
			relax += itw->weight * Color(itw->vertex);
			++itw;
		}
		relax.Clamp(0,1);
		#ifdef _HSI_COLOR_
			levels[current_level_number].color_details[v][0] = HsiDifference(Rgb2Hsi(Color(v)), Rgb2Hsi(relax));
			levels[current_level_number].color_details[v][1] = 0;
			levels[current_level_number].color_details[v][2] = 0;
		#else
			levels[current_level_number].color_details[v] = Color(v) - relax;
		#endif
	}

	// Texture detail
	if( use_texture )
	{
	//	Vector2d relax_texture_coord(0,0);
		relax = 0.0;
		itw = levels[current_level_number].subdivision_weights[v].begin();
		while( itw != levels[current_level_number].subdivision_weights[v].end() )
		{
		//	relax_texture_coord += itw->weight * Texture(itw->vertex);
			relax += itw->weight * Pixel(itw->vertex);
			++itw;
		}
		levels[current_level_number].texture_details[v] = Pixel(v) - relax;
	//	levels[current_level_number].texture_details[v] = Pixel(v) - Pixel(relax_texture_coord);
	//	levels[current_level_number].texture_details[v][0] = Texture(v)[0] - relax_texture_coord[0];
	//	levels[current_level_number].texture_details[v][1] = Texture(v)[1] - relax_texture_coord[1];
	//	levels[current_level_number].texture_details[v][2] = 0;
	}

	return true;
}


//--
//
// ReconstructVertex
//
//--
Vector3d MultiresolutionMesh::ReconstructVertex( int v )
{
	// Predict vertex position using subdivision weights
	Vector3d temp(0,0,0);
	SubdivisionWeightList::const_iterator it = levels[current_level_number].subdivision_weights[v].begin();
	while( it != levels[current_level_number].subdivision_weights[v].end() )
	{
		temp += it->weight * Vertex(it->vertex);
		++it;
	}
	// Add details
	return temp += levels[current_level_number].geometric_details[v];
}


//--
//
// ReconstructVertexColor
//
//--
Vector3d MultiresolutionMesh::ReconstructVertexColor( int v )
{
	// Predict vertex position using subdivision weights
	Vector3d temp(0,0,0);
	SubdivisionWeightList::const_iterator it = levels[current_level_number].subdivision_weights[v].begin();
	while( it != levels[current_level_number].subdivision_weights[v].end() )
	{
		temp += it->weight * Color(it->vertex);
		++it;
	}
	// Add details
	return temp += levels[current_level_number].color_details[v];
}


//--
//
// ComputeMeanCurvature
//
//--
void MultiresolutionMesh::ComputeMeanCurvature()
{
	for( int i=0; i<VertexNumber(); i++ )
	{
		if( IsProgressiveVertexValid(i) == false )
		{
			continue;
		}
		mean_curvature[i] = MeanCurvatureNormal(this, i).Length() / 2.0;
	}
}


//--
//
// ComputeGaussianCurvature
//
//--
void MultiresolutionMesh::ComputeGaussianCurvature()
{
	for( int i=0; i<VertexNumber(); i++ )
	{
		if( IsProgressiveVertexValid(i) == false )
		{
			continue;
		}
		gaussian_curvature[i] = GaussianCurvature(this, i);
	}
}


//--
//
// Pixel
//
//--
Vector3d MultiresolutionMesh::Pixel(const double& s, const double& t)
{
	const int& w = texture.Width();
	const int& h = texture.Height();
	double r = texture.Buffer((int)(s*((double)w-1.0))*3+(int)(t*((double)h-1.0))*w*3);
	double g = texture.Buffer((int)(s*((double)w-1.0))*3+(int)(t*((double)h-1.0))*w*3+1);
	double b = texture.Buffer((int)(s*((double)w-1.0))*3+(int)(t*((double)h-1.0))*w*3+2);
	return Vector3d( r/255.0, g/255.0, b/255.0 );
}


//--
//
// Pixel
//
//--
Vector3d MultiresolutionMesh::Pixel(const Vector2d& coord)
{
	return Pixel( coord[0], coord[1] );
}


//--
//
// Pixel
//
//--
Vector3d MultiresolutionMesh::Pixel(int v)
{
	return Pixel( Texture(v) );
}


//--
//
// BuildProgressiveEdgeList
//
//--
void MultiresolutionMesh::BuildProgressiveEdgeList()
{
	// Reset the edge list
	edge_list.clear();
	// For every face
	for( int i=0; i<FaceNumber(); i++ )
	{
		// Test if the face is valid
		if( !IsProgressiveFaceValid(i) ) continue;
		// Shortcuts to vertex indices of the current face
		const int& a = Face(i, 0);
		const int& b = Face(i, 1);
		const int& c = Face(i, 2);
		// Edge (a,b)
		if( a < b ) edge_list.push_back(Vector2i(a, b));
		else edge_list.push_back(Vector2i(b, a));
		// Edge (b,c)
		if( b < c ) edge_list.push_back(Vector2i(b, c));
		else edge_list.push_back(Vector2i(c, b));
		// Edge (c,a)
		if( c < a ) edge_list.push_back(Vector2i(c, a));
		else edge_list.push_back(Vector2i(a, c));
	}
	// Sort the edge list
	edge_list.sort();
	// Remove redundant edges
	edge_list.unique();
}


//--
//
// ComputeEdgeCoefficients
//
//--
void MultiresolutionMesh::ComputeEdgeCoefficients()
{
	double coef;
	ConstNeighborIterator itv;
	// Reset the edge coefficients of the current level
//	levels[current_level_number].weights.clear();
	// Iterate through the edge list
	std::list<Vector2i>::const_iterator ite = edge_list.begin();
	while( ite != edge_list.end() ) {
		// Shortcut to the two vertices composing the edge
		const int& va = (*ite)[0];
		const int& vb = (*ite)[1];
		// Init the coef
		coef = 0;
		// Look into the neighborhood of va to find common
		// side vertices of (va,vb)
		itv = NeighborVertices(va).begin();
		while( itv != NeighborVertices(va).end() ) {
			if( !FindNeighborVertex(vb, *itv) ) continue;
			if( IsColinear(Vertex(*itv), Vertex(va), Vertex(vb)) ) continue;
			coef += Cotan( Vertex(*itv), Vertex(va), Vertex(vb) );
			++itv;
		}
		// Save the coef in the current level
//		levels[current_level_number].weights[*ite] = coef;
		// Next edge
		++ite;
	}
}
