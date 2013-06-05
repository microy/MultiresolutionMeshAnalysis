/***************************************************************************
                                Visualizer.h
                             -------------------
    update               : 2009/02/19
    copyright            : (C) 2003-2009 by Michaël Roy
    email                : michael.roy@u-bourgogne.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef _VISUALIZER_
#define _VISUALIZER_

#include "MultiresolutionMesh.h"
#include "TextureImage.h"
#include "Viewport.h"

#include <GL/gl.h>
#include <GL/glu.h>

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>


//--
//
// Visualizer
//
//--
// Manage the display of a multiresolution mesh
// Inherit from Fl_Gl_Window (FLTK OpenGL window)
class Visualizer : public Fl_Gl_Window
{

	//--
	//
	// Enumerations
	//
	//--
	protected :

		// Mouse buttons
		enum MouseButton {
			LEFT_BUTTON = 1,
			MIDDLE_BUTTON = 2,
			RIGHT_BUTTON = 3
		};

		// Motion states
		enum MotionState {
			MOTION_NONE,
			MOTION_ROTATION,
			MOTION_TRANSLATION_XY,
			MOTION_TRANSLATION_Z
		};


	//--
	//
	// Member variables
	//
	//--
	protected :

		// Mesh related variables
		MultiresolutionMesh* mesh;

		// Texture image
		TextureImage texture;

		// Display options
		bool wireframe_enabled;
		bool solid_enabled;
		bool vertex_enabled;
		bool color_enabled;
		bool texture_enabled;
		bool smooth_enabled;
		bool normal_enabled;
		bool antialiasing_enabled;
		bool geometrical_details_enabled;
		bool normal_details_enabled;
		bool color_details_enabled;
		bool texture_details_enabled;
		bool simplification_enabled;
		bool color_bar_enabled;
		bool bounding_box_enabled;

		// Scaling variables
		double average_edge_length;
		double vector_scale;

		// Current mouse motion state
		MotionState motion_state;

		// OpenGL list number for a sphere
		int sphere_list_number;

		// Transformation related variables
		double view_angle;
		double aspect_ratio;
		Vector3d model_position;
		double scale_factor;
		Vector3d translations;
		Vector3d rotation_axis;
		double rotation_angle;
		Vector2i previous_mouse_position;
		Vector3d previous_trackball_position;
		double trackball_transform[4][4];
		double external_rotation_x;
		double external_rotation_y;
		double external_rotation_z;

		int frustum_list_number;

		
	//--
	//
	// Member functions
	//
	//--
	public :

		// Default constructor
		Visualizer(int position_x=100, int position_y=100, int window_width=400, int window_height=300, const char* label = "");

		// Destructor
		~Visualizer();

		// Set the model to visualize
		bool SetModel(MultiresolutionMesh* m);

		void ModelCorrections(Vector3d position, double scale) {
			model_position = position;
			scale_factor = scale;
		}

		// Update the model
		void UpdateModel(bool update_colors=true);

		// Overloaded function from Fl_Gl_Window inheritance
		void show();

		// Get the viewport
		Viewport GetViewport();

		// Option functions
		void EnableSolid(bool enable);
		void EnableWireframe(bool enable);
		void EnableVertices(bool enable);
		void EnableNormals(bool enable);
		void EnableColors(bool enable);
		void EnableTexture(bool enable);
		void EnableSmooth(bool enable);
		void EnableAntialiasing(bool enable);
		void EnableDetailsGeometrical(bool enable);
		void EnableDetailsNormal(bool enable);
		void EnableDetailsColor(bool enable);
		void EnableDetailsTexture(bool enable);
		void SetBackgroundColor(Vector3d color);
		void ToggleSimplification(bool enable);
		void ToggleColorBar(bool enable);
		void SetTransformation(double rotx, double roty, double rotz);

	protected :

		// Overloaded functions from Fl_Gl_Window inheritance
		void resize(int x, int y, int w, int h);
		int handle(int);
		void draw();

		// OpenGL initialization
		void InitializeOpenGL();

		// OpenGL viewport resizing
		void ResizeOpenGL( int width, int height );

		// Mouse callback function
		bool Mouse( int button, int x, int y );

		// Mouse motion callback function
		bool Motion( int x, int y );

		// Keyboard callback function
		bool Keyboard( int key );

		// Trackball function
		Vector3d TrackballMapping( int x, int y ) const;
		void InitializeTrackball();
		Vector3d ApplyTrackballTransformation(const Vector3d& p) const;

		// Display functions
		void DrawViewingFrustum();
		void DrawMesh();
		void DrawColorMesh();
		void DrawTextureMesh();
		void DrawSmoothMesh();
		void DrawSmoothColorMesh();
		void DrawSmoothTextureMesh();
		void DrawVertices();
		void DrawNormals();
		void DrawGeometricalDetails();
		void DrawNormalDetails();
		void DrawColorDetails();
		void DrawTextureDetails();
		void DrawColorBar();

};


#endif // _VISUALIZER_
