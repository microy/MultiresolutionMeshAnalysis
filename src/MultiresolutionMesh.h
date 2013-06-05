/***************************************************************************
                            MultiresolutionMesh.h
                             -------------------
    update               : 2003/11/11
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

#ifndef _MULTIRESOLUTIONMESH_
#define _MULTIRESOLUTIONMESH_

#include "ProgressiveMesh.h"
#include "ResolutionLevel.h"
#include "TextureImage.h"

#include <list>
#include <vector>



//#define HSI_COLOR


//--
//
// MultiresolutionMesh
//
//--
class MultiresolutionMesh : public ProgressiveMesh
{

	//--
	//
	// Member Data
	//
	//--
	protected :

		int current_level_number; // Current resolution level
		std::vector<ResolutionLevel> levels; // Array of resolution level information
		std::vector<double> mean_curvature;
		std::vector<double> gaussian_curvature;
		TextureImage texture; // The texture image of the model
		std::vector<int> vertex_map; // Vertex map
		std::list<Vector2i> edge_list; // Edge list
		bool use_color;
		bool use_texture;
		bool use_normal;

	//--
	//
	// Member Functions
	//
	//--
	public :

		//
		// Constructor / Destructor
		//
		inline MultiresolutionMesh():current_level_number(-1) {}
		inline ~MultiresolutionMesh() {}

		//
  		// Multiresolution Analysis
		//
		bool Analysis(int base_level = 0);
		bool Synthesis();

		// Level decomposition/reconstruction
		bool DecomposeCurrentLevel();
		bool ReconstructCurrentLevel();

		// Level management
		bool GoToLevel(int i);
		bool GoToBaseLevel();
		bool GoToInitialLevel();

		//
		// Level interface
		//
		inline int LevelNumber() const {
			return levels.size();
		}

		inline int CurrentLevelNumber() const {
			return current_level_number;
		}

		inline std::vector<ResolutionLevel>& Levels() {
			return levels;
		}

		inline const std::vector<ResolutionLevel>& Levels() const {
			return levels;
		}

		inline ResolutionLevel& Level(int i) {
			assert( (i>=0) && (i<(int)levels.size()) );
			return levels[i];
		}

		inline const ResolutionLevel& Level(int i) const {
			assert( (i>=0) && (i<(int)levels.size()) );
			return levels[i];
		}

		inline ResolutionLevel& CurrentLevel() {
			assert( (current_level_number>=0) && (current_level_number<(int)levels.size()) );
			return levels[current_level_number];
		}

		inline const ResolutionLevel& CurrentLevel() const {
			assert( (current_level_number>=0) && (current_level_number<(int)levels.size()) );
			return levels[current_level_number];
		}

		Vector3d Pixel(int v);

		//
		// Vertex Map Interface
		//
		inline int& VertexMap(int i) {
			assert( (i>=0) && (i<(int)vertex_map.size()) );
			return vertex_map[i];
		}

		inline const int& VertexMap(int i) const {
			assert((i>=0) && (i<(int)vertex_map.size()) );
  			return vertex_map[i];
		}



	protected :

		bool PredictVertex( int v, VertexType vt );
		Vector3d ReconstructVertex( int v );
		Vector3d ReconstructVertexColor( int v );

		void ComputeMeanCurvature();
		void ComputeGaussianCurvature();

		Vector3d Pixel(const double& s, const double& t);
		Vector3d Pixel(const Vector2d& coord);

		void BuildProgressiveEdgeList();
		void ComputeEdgeCoefficients();

};



#endif // _MULTIRESOLUTIONMESH_

