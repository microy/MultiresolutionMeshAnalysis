/***************************************************************************
                               Visualizer.cpp
                             -------------------
    update               : 2004/10/16
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

#include "Visualizer.h"
#include "GraphicalInterface.h"

#include <cmath>


//--
//
// Visualizer
//
//--
// Default constructor
Visualizer::Visualizer(int position_x, int position_y, int window_width, int window_height, const char* label)
	// Initialize base class
	: Fl_Gl_Window(position_x, position_y, window_width, window_height, label),
	// Initialize member variables
	mesh(0),
	wireframe_enabled(false), solid_enabled(true), vertex_enabled(false),
	color_enabled(false), texture_enabled(false), smooth_enabled(false),
	normal_enabled(false), antialiasing_enabled(false), geometrical_details_enabled(false),
	normal_details_enabled(false), color_details_enabled(false), texture_details_enabled(false),
	simplification_enabled(false), color_bar_enabled(false), bounding_box_enabled(false),
	average_edge_length(0.0), vector_scale(0.0), motion_state(MOTION_NONE),
	sphere_list_number(0), view_angle(45.0), model_position(0.0,0.0,0.0),
	scale_factor(1.0), translations(0.0,0.0,0.0), rotation_axis(0.0,0.0,0.0),
	rotation_angle(0.0), previous_mouse_position(0,0), previous_trackball_position(0.0,0.0,0.0),
	external_rotation_x(0.0), external_rotation_y(0.0), external_rotation_z(0.0)
{
	// Initialize OpenGL display mode
	// RGB colors + Double buffer + Depth test
	mode(FL_RGB|FL_DOUBLE|FL_DEPTH);

	// Initialize the trackball tranformation matrix
	InitializeTrackball();
}


//--
//
// ~Visualizer
//
//--
// Destructor
Visualizer::~Visualizer()
{
}


//--
//
// InitializeOpenGL
//
//--
// Initialize OpenGL
void Visualizer::InitializeOpenGL()
{
	// OpenGL color configuration
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClearDepth(1.0);

	// Back face culling
	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	// Depth test
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);

	// Normal vector normalization
	// Due to glScale usage
	glEnable(GL_NORMALIZE);

	// Anti-aliasing settings
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glEnable( GL_BLEND );

	// Model rendering settings
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	// Color
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// Light
	int LightPos[4] = {0,0,3,1};
	glLightiv(GL_LIGHT0,GL_POSITION,LightPos);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	// Model matrix
	glMatrixMode(GL_MODELVIEW);
}

//--
//
// SetModel
//
//--
bool Visualizer::SetModel(MultiresolutionMesh* m)
{
	// Save mesh address
	mesh = m;

	// Check mesh validity
	if( mesh == 0 ) return false;
	if( mesh->VertexNumber() == 0 ) return false;
	if( mesh->FaceNumber() == 0 ) return false;
	if(  mesh->LevelNumber()!=0 && mesh->CurrentLevelNumber()!=-1 ) mesh->Synthesis();

	int count = 0;
	average_edge_length = 0.0;
	// Compute average edge length
	for( int i=0; i<mesh->VertexNumber(); i++ )
	{
		ConstNeighborIterator it = mesh->NeighborVertices(i).begin();
		while( it != mesh->NeighborVertices(i).end() )
		{
			if( *it > i )
			{
				average_edge_length += (mesh->Vertex(i)-mesh->Vertex(*it)).Length();
				count++;
			}
			++it;
		}
	}
	vector_scale = average_edge_length /= (double)count;

	// Generate a sphere in a list
	sphere_list_number = glGenLists(1);
	glNewList(sphere_list_number, GL_COMPILE_AND_EXECUTE);
		GLUquadricObj* sphere = gluNewQuadric();
		gluSphere( sphere, average_edge_length/5.0, 6, 6 );
	glEndList();

	// Color management
	if(mesh->ColorNumber() == mesh->VertexNumber()) color_enabled = true;
	else color_enabled = false;

	// Texture management
	if(mesh->TextureNumber()!=0 && mesh->TextureName().empty()==false) {
		texture_enabled = true;
		texture.ReadFile( mesh->TextureName() );
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
		glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
		glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
		gluBuild2DMipmaps( GL_TEXTURE_2D, 3, texture.Width(), texture.Height(), GL_RGB, GL_UNSIGNED_BYTE, texture.Buffer());
	}
	else {
		texture_enabled = false;
	}

	// Compute normals
//	mesh->ComputeProgressiveFaceNormals();
//	mesh->ComputeProgressiveVertexNormals();

	mesh->ComputeBoundingBox();
	model_position = mesh->BoundingBox().Barycenter();
	Vector3d model_size = mesh->BoundingBox().Length();
	scale_factor = 1.0 / std::max( std::max(model_size[0],model_size[1]), model_size[2] );

	// Initialize transformations
	translations = 0;
	rotation_angle = 0;

	// Initialize the trackball tranformation matrix
	InitializeTrackball();

	// Update the view
	redraw();

	return true;
}


//--
//
// UpdateModel
//
//--
void Visualizer::UpdateModel(bool update_colors)
{
	// Check if we have a mesh
	if( !mesh ) return;

	// Color management
//	if(update_colors && mesh->ColorNumber()==mesh->VertexNumber()) color_enabled = true;
//	else color_enabled = false;


	// Update the display
	redraw();
}


//--
//
// DrawViewingFrustum
//
//--
void Visualizer::DrawViewingFrustum()
{
	glCallList(frustum_list_number);
}


//--
//
// DrawMesh
//
//--
void Visualizer::DrawMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->FaceNormal(f) );

			glVertex3dv( mesh->Vertex(f,0) );
			glVertex3dv( mesh->Vertex(f,1) );
			glVertex3dv( mesh->Vertex(f,2) );
		}

	glEnd();
}


//--
//
// DrawColorMesh
//
//--
void Visualizer::DrawColorMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->FaceNormal(f) );

			glColor3dv( mesh->Color(f,0) );
			glVertex3dv( mesh->Vertex(f,0) );

			glColor3dv( mesh->Color(f,1) );
			glVertex3dv( mesh->Vertex(f,1) );

			glColor3dv( mesh->Color(f,2) );
			glVertex3dv( mesh->Vertex(f,2) );
		}

	glEnd();
}


//--
//
// DrawTextureMesh
//
//--
void Visualizer::DrawTextureMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->FaceNormal(f) );

			glTexCoord2dv( mesh->Texture(f,0) );
			glVertex3dv( mesh->Vertex(f,0) );

			glTexCoord2dv( mesh->Texture(f,1) );
			glVertex3dv( mesh->Vertex(f,1) );

			glTexCoord2dv( mesh->Texture(f,2) );
			glVertex3dv( mesh->Vertex(f,2) );
		}

	glEnd();
}


//--
//
// DrawSmoothMesh
//
//--
void Visualizer::DrawSmoothMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->VertexNormal(f,0) );
			glVertex3dv( mesh->Vertex(f,0) );

			glNormal3dv( mesh->VertexNormal(f,1) );
			glVertex3dv( mesh->Vertex(f,1) );

			glNormal3dv( mesh->VertexNormal(f,2) );
			glVertex3dv( mesh->Vertex(f,2) );
		}

	glEnd();
}


//--
//
// DrawSmoothColorMesh
//
//--
void Visualizer::DrawSmoothColorMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->VertexNormal(f,0) );
			glColor3dv( mesh->Color(f,0) );
			glVertex3dv( mesh->Vertex(f,0) );

			glNormal3dv( mesh->VertexNormal(f,1) );
			glColor3dv( mesh->Color(f,1) );
			glVertex3dv( mesh->Vertex(f,1) );

			glNormal3dv( mesh->VertexNormal(f,2) );
			glColor3dv( mesh->Color(f,2) );
			glVertex3dv( mesh->Vertex(f,2) );
		}

	glEnd();
}


//--
//
// DrawSmoothTextureMesh
//
//--
void Visualizer::DrawSmoothTextureMesh()
{
	glBegin(GL_TRIANGLES);

		for( int f=0; f<mesh->FaceNumber(); f++ ) {
			if( mesh->IsProgressiveFaceValid(f) == false ) continue;

			glNormal3dv( mesh->VertexNormal(f,0) );
			glTexCoord2dv( mesh->Texture(f,0) );
			glVertex3dv( mesh->Vertex(f,0) );

			glNormal3dv( mesh->VertexNormal(f,1) );
			glTexCoord2dv( mesh->Texture(f,1) );
			glVertex3dv( mesh->Vertex(f,1) );

			glNormal3dv( mesh->VertexNormal(f,2) );
			glTexCoord2dv( mesh->Texture(f,2) );
			glVertex3dv( mesh->Vertex(f,2));
		}

	glEnd();
}


//--
//
// DrawVertices
//
//--
void Visualizer::DrawVertices()
{
	// If the mesh has been simplified
	if( simplification_enabled ) {
		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( !mesh->IsProgressiveVertexValid(i) ) continue;
			const Vector3d& v = mesh->Vertex(i);
			glPushMatrix();
				glColor3f( 0.8, 0.4, 0.4 );
				glTranslated( v[0], v[1], v[2] );
				glCallList( sphere_list_number );
			glPopMatrix();
		}
		return;
	}

	// Shortcut to the current resolution level
	const ResolutionLevel& level = mesh->CurrentLevel();

	for( int i=0; i<mesh->VertexNumber(); i++ ) {
		if( !mesh->IsProgressiveVertexValid(i) ) continue;
		const Vector3d& v = mesh->Vertex(i);
		switch( level.vertex_types[i] ) {
		//	case UNDEFINED_VERTEX :
		//		glPushMatrix();
		//			glColor3f( 0.1, 0.1, 0.1 );
		//			glTranslatef( v[0], v[1], v[2] );
		//			glCallList( sphere_list_number );
		//		glPopMatrix();
		//		break;

			case ODD_VERTEX :
				glPushMatrix();
					glColor3f( 0.8, 0.4, 0.4 );
					glTranslated( v[0], v[1], v[2] );
					glCallList( sphere_list_number );
				glPopMatrix();
				break;

			case EVEN_VERTEX :
				glPushMatrix();
					glColor3f( 0.4, 0.8, 0.4 );
					glTranslated( v[0], v[1], v[2] );
					glCallList( sphere_list_number );
				glPopMatrix();
				break;
				
			default:
				break;
		}
	}
}


//--
//
// DrawNormals
//
//--
void Visualizer::DrawNormals()
{
	Vector3d point;

	glBegin(GL_LINES);

		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( mesh->IsProgressiveVertexValid(i) == false )  continue;

			point = mesh->Vertex(i) + mesh->VertexNormal(i)*vector_scale*10.0;
			glVertex3dv(mesh->Vertex(i));
			glVertex3dv(point);
		}

	glEnd();
}


//--
//
// DrawGeometricalDetails
//
//--
void Visualizer::DrawGeometricalDetails()
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->CurrentLevel();
	Vector3d point;

	glBegin(GL_LINES);

		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( mesh->IsProgressiveVertexValid(i) == false ) continue;
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;

			point = mesh->Vertex(i) + level.geometric_details[i]*vector_scale*500.0;
			glVertex3dv(mesh->Vertex(i));
			glVertex3dv(point);
		}

	glEnd();
}


//--
//
// DrawNormalDetails
//
//--
void Visualizer::DrawNormalDetails()
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->CurrentLevel();
	Vector3d point;

	glBegin(GL_LINES);

		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( mesh->IsProgressiveVertexValid(i) == false ) continue;
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;

			point = mesh->Vertex(i) + level.normal_details[i] * vector_scale;
			glVertex3dv(mesh->Vertex(i));
			glVertex3dv(point);
		}

	glEnd();
}


//--
//
// DrawColorDetails
//
//--
void Visualizer::DrawColorDetails()
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->CurrentLevel();
	Vector3d point;

	glBegin(GL_LINES);

		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( mesh->IsProgressiveVertexValid(i) == false ) continue;
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;

			point = mesh->Vertex(i) + level.color_details[i] * vector_scale;
			glVertex3dv(mesh->Vertex(i));
			glVertex3dv(point);
		}

	glEnd();
}


//--
//
// DrawTextureDetails
//
//--
void Visualizer::DrawTextureDetails()
{
	// Shortcut to the current resolution level
	ResolutionLevel& level = mesh->CurrentLevel();
	Vector3d point;

	glBegin(GL_LINES);

		for( int i=0; i<mesh->VertexNumber(); i++ ) {
			if( mesh->IsProgressiveVertexValid(i) == false ) continue;
			if( level.vertex_types[i] == UNDEFINED_VERTEX ) continue;

			point = mesh->Vertex(i) + level.texture_details[i] * vector_scale;
			glVertex3dv(mesh->Vertex(i));
			glVertex3dv(point);
		}

	glEnd();
}


//--
//
// draw
//
//--
void Visualizer::draw()
{
	// Initialize display
	glLoadIdentity();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// No mesh
	if( !mesh ) return;

	// Translations
	glTranslated(translations[0], translations[1], translations[2]);

	// Trackball rotation
	glPushMatrix();
		glLoadIdentity();
		glRotated(rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2]);
		glMultMatrixd((GLdouble*)trackball_transform);
		glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble*)trackball_transform);
	glPopMatrix();
	glMultMatrixd((GLdouble*)trackball_transform);
	
	// External rotation (fixed with the GUI)
	glRotated(external_rotation_x, 1, 0, 0);
	glRotated(external_rotation_y, 0, 1, 0);
	glRotated(external_rotation_z, 0, 0, 1);

	// Model corrections
	glScaled(scale_factor,scale_factor,scale_factor);
	glTranslated(-model_position[0], -model_position[1], -model_position[2]);


	// Solid model
	if( solid_enabled ) {
		glColor3f(0.7,0.7,0.7);
		// Smooth shading
		if( smooth_enabled ) {
			if( texture_enabled ) DrawSmoothTextureMesh();
			else if( color_enabled ) DrawSmoothColorMesh();
			else DrawSmoothMesh();
		}
		// Flat shading
		else {
			if( texture_enabled ) DrawTextureMesh();
			else if( color_enabled ) DrawColorMesh();
			else DrawMesh();
		}
	}

	// Invisible model (for wireframe without hidden faces)
	else {
		glDisable(GL_LIGHTING);
		glDisable(GL_COLOR_MATERIAL);
		glColor3f(1,1,1);
		DrawSmoothMesh();
		glEnable(GL_LIGHTING);
		glEnable(GL_COLOR_MATERIAL);
	}

	// Wireframe model
	if( wireframe_enabled ) {
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		glPolygonOffset(-1.0,2.0);
		glEnable(GL_POLYGON_OFFSET_LINE);
		glColor3f(1,0,0);
		DrawMesh();
		glDisable(GL_POLYGON_OFFSET_LINE);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glEnable(GL_LIGHTING);
	}

	// Normals
	if( normal_enabled ) {
		glDisable(GL_LIGHTING);
		glColor3f(0.2,0.8,0.4);
		DrawNormals();
		glEnable(GL_LIGHTING);
	}

	// Vertices
	if( vertex_enabled ) {
		DrawVertices();
	}

	// Geometrical detail vectors
	if( geometrical_details_enabled ) {
		glDisable(GL_LIGHTING);
		glColor3f(0.8,0.2,0.4);
		DrawGeometricalDetails();
		glEnable(GL_LIGHTING);
	}

	// Normal detail vectors
	else if( normal_details_enabled ) {
		glDisable(GL_LIGHTING);
		glColor3f(0.8,0.2,0.4);
		DrawNormalDetails();
		glEnable(GL_LIGHTING);
	}

	// Color detail vectors
	else if( color_details_enabled ) {
		glDisable(GL_LIGHTING);
		glColor3f(0.8,0.2,0.4);
		DrawColorDetails();
		glEnable(GL_LIGHTING);
	}

	// Texture detail vectors
	else if( texture_details_enabled ) {
		glDisable(GL_LIGHTING);
		glColor3f(0.8,0.2,0.4);
		DrawTextureDetails();
		glEnable(GL_LIGHTING);
	}

	// DrawViewingFrustum
	if( simplification_enabled ) {
		DrawViewingFrustum();
	}

	// Bounding box
	if( bounding_box_enabled ) {
		// No face culling
		glDisable(GL_CULL_FACE);
		// Transparent color
		glColor4f(0.4, 0.4, 0.8, 0.2);
		// Draw bb
		mesh->BoundingBox().Draw();
		glEnable(GL_CULL_FACE);
	}

	// Color bar
	if( color_bar_enabled ) DrawColorBar();

	// Initialize rotation
	rotation_axis = Vector3d(0,0,0);
	rotation_angle = 0;

	// Do not swap buffers
	// FLTK takes care of it
}


//--
//
// show
//
//--
void Visualizer::show()
{
	// Show the window
	Fl_Gl_Window::show();

	// Activate OpenGL context for this window
	make_current();

	// Initialize OpenGL
	InitializeOpenGL();

	// Setup the OpenGL viewport
	ResizeOpenGL(w(), h());
}


//--
//
// resize
//
//--
void Visualizer::resize(int x, int y, int w, int h)
{
	// Resize the window
	Fl_Gl_Window::resize(x, y, w, h);

	// Resize the OpenGL viewport
	ResizeOpenGL(w, h);
}


//--
//
// ResizeOpenGL
//
//--
void Visualizer::ResizeOpenGL(int window_width, int window_height)
{
	aspect_ratio = (double)window_width/window_height;
	glViewport(0, 0, window_width, window_height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective((float)view_angle, (float)aspect_ratio, 0.1, 1000.0);
	gluLookAt(0,0,2,0,0,0,0,1,0);
	glMatrixMode(GL_MODELVIEW);
}


//--
//
// InitializeTrackball
//
//--
// Set the trackball transformation matrix to identity
void Visualizer::InitializeTrackball()
{
	for( int i=0; i<4; i++ ) {
		for( int j=0; j<4; j++ ) {
			if( i != j ) trackball_transform[i][j] = 0.0;
			else trackball_transform[i][j] = 1.0;
		}
	}
}


//--
//
// TrackballMapping
//
//--
// Adapted from Nate Robins programs
// http://www.xmission.com/~nate
Vector3d Visualizer::TrackballMapping( int x, int y ) const
{
	static Vector3d v(0,0,0);
	static const double pi_2 = M_PI / 2.0;

	v[0] = ( 2.0 * (double)x - (double)w() ) / (double)w();
	v[1] = ( (double)h() - 2.0 * (double)y ) / (double)h();
	v[2] = 0.0;

	double  d = v.Length();
	if( d > 1.0 ) d = 1.0;

	v[2] = cos(pi_2 * d);

	return v.Normalize();
}


//--
//
// handle
//
//--
int Visualizer::handle(int event)
{
	switch(event) {
		// Mouse down event
		case FL_PUSH :
			Mouse(Fl::event_button(), Fl::event_x(), Fl::event_y());
			break;

		// Mouse up event
		case FL_RELEASE :
			motion_state = MOTION_NONE;
			cursor(FL_CURSOR_DEFAULT);
			break;

		// Mouse moved while down event
		case FL_DRAG :
			Motion(Fl::event_x(), Fl::event_y());
			break;

		// Keyboard event
		// Return 1 if you understand/use the keyboard event, 0 otherwise...
		case FL_KEYBOARD :
		case FL_SHORTCUT :
			if( Keyboard( Fl::event_key() ) ) break;
			return 0;

		default:
			// pass other events to the base class...
			return Fl_Gl_Window::handle(event);
	}

	// Update the display
	redraw();

	// Ok, I handled the event
	return 1;
}


//--
//
// Mouse
//
//--
bool Visualizer::Mouse(int button, int x, int y)
{
	switch( button ) {

		// Left button
		case FL_LEFT_MOUSE :
			motion_state = MOTION_ROTATION;
			// Trackball Rotation
			previous_trackball_position = TrackballMapping( x, y );
			// Change window cursor
			cursor(FL_CURSOR_HAND);
			break;

		// Middle button
		case FL_MIDDLE_MOUSE :
			motion_state = MOTION_TRANSLATION_XY;
			previous_mouse_position = Vector2i( x, y );
			cursor(FL_CURSOR_MOVE);
			break;

		// Right button
		case FL_RIGHT_MOUSE :
			motion_state = MOTION_TRANSLATION_Z;
			previous_mouse_position = Vector2i( x, y );
			cursor(FL_CURSOR_NS);
			break;

		default :
			motion_state = MOTION_NONE;
			cursor(FL_CURSOR_DEFAULT);
			return false;
	}

	return true;
}


//--
//
// Motion
//
//--
bool Visualizer::Motion(int x, int y)
{
	static Vector3d current_position;

	switch( motion_state ) {

		// Trackball Rotation
		case MOTION_ROTATION :
			// Map the mouse position to a logical
			current_position = TrackballMapping( x, y );
			// Rotate about the axis that is perpendicular to the great circle connecting the mouse movements.
			rotation_axis = previous_trackball_position ^ current_position;
			rotation_angle = 90.0 * Length(current_position-previous_trackball_position) * 1.5;
			previous_trackball_position = current_position;
			break;

		case MOTION_TRANSLATION_XY :
			translations[0] -= ((float)previous_mouse_position[0]-x)*0.01;
			translations[1] += ((float)previous_mouse_position[1]-y)*0.01;
			previous_mouse_position = Vector2i( x, y );
			break;

		case MOTION_TRANSLATION_Z :
			translations[2] += ((float)previous_mouse_position[1]-y)*0.01;
			previous_mouse_position = Vector2i( x, y );
			break;

		default :
			return false;
	}

	return true;
}


//--
//
// Keyboard
//
//--
bool Visualizer::Keyboard(int key)
{
	if( !mesh ) return false;
	switch( key ) {
		case 'r' :
		case 'R' :
			translations = 0;
			external_rotation_x = 0;
			external_rotation_y = 0;
			external_rotation_z = 0;
			// Initialize the trackball tranformation matrix
			InitializeTrackball();
			break;

		case 'b' :
		case 'B' :
			bounding_box_enabled = !bounding_box_enabled;
			break;
			
		case 'g' :
		case 'G' :
			geometrical_details_enabled = !geometrical_details_enabled;
			break;
		
		case '+' :
		case '=' :
			vector_scale *= 2.0;
			break;
			
		case '-' :
			vector_scale /= 2.0;
			break;

		// Coaser resolution level
		case FL_Up :
		case FL_Page_Up :
			if(simplification_enabled || !mesh->LevelNumber()) break;
			if(mesh->CurrentLevelNumber() < mesh->LevelNumber()) {
				mesh->DecomposeCurrentLevel();
				mesh->ComputeProgressiveFaceNormals();
				mesh->ComputeProgressiveVertexNormals();
				((GraphicalInterface*)user_data())->UpdateStatus();
			}
			break;

		// Finer resolution level
		case FL_Down :
		case FL_Page_Down :
			if(simplification_enabled || !mesh->LevelNumber()) break;
			if(mesh->CurrentLevelNumber() > 0) {
				mesh->ReconstructCurrentLevel();
				mesh->ComputeProgressiveFaceNormals();
				mesh->ComputeProgressiveVertexNormals();
				((GraphicalInterface*)user_data())->UpdateStatus();
			}
			break;

		default :
			return false;
	}
	return true;
}


//--
//
// EnableSolid
//
//--
void Visualizer::EnableSolid(bool enable)
{
	solid_enabled = enable;
	redraw();
}


//--
//
// EnableWireframe
//
//--
void Visualizer::EnableWireframe(bool enable)
{
	wireframe_enabled = enable;
	redraw();
}


//--
//
// EnableVertices
//
//--
void Visualizer::EnableVertices(bool enable)
{
	vertex_enabled = enable;
	redraw();
}


//--
//
// EnableNormals
//
//--
void Visualizer::EnableNormals(bool enable)
{
	normal_enabled = enable;
	redraw();
}


//--
//
// EnableColors
//
//--
void Visualizer::EnableColors(bool enable)
{
	if( !mesh ) return;
	if( mesh->ColorNumber() == mesh->VertexNumber() ) {
		color_enabled = enable;
	}
	redraw();
}


//--
//
// EnableTexture
//
//--
void Visualizer::EnableTexture(bool enable)
{
	if( !mesh ) return;
	if( (mesh->TextureNumber()!=0) && (mesh->TextureName().empty()==false) ) {
		texture_enabled = enable;
	}
	redraw();
}


//--
//
// EnableSmooth
//
//--
void Visualizer::EnableSmooth(bool enable)
{
	smooth_enabled = enable;
	redraw();
}


//--
//
// EnableAntialiasing
//
//--
void Visualizer::EnableAntialiasing(bool enable)
{
	antialiasing_enabled = enable;
	// Enable anti-aliasing
	if( antialiasing_enabled == true ) {
		glEnable( GL_POINT_SMOOTH );
		glEnable( GL_LINE_SMOOTH );
	}
	// Disable anti-aliasing
	else {
		glDisable( GL_LINE_SMOOTH );
		glDisable( GL_POINT_SMOOTH );
	}
	redraw();
}


//--
//
// EnableDetailsGeometrical
//
//--
void Visualizer::EnableDetailsGeometrical(bool enable)
{
	geometrical_details_enabled = enable;
	redraw();
}


//--
//
// EnableDetailsNormal
//
//--
void Visualizer::EnableDetailsNormal(bool enable)
{
	normal_details_enabled = enable;
	redraw();
}


//--
//
// EnableDetailsColor
//
//--
void Visualizer::EnableDetailsColor(bool enable)
{
	color_details_enabled = enable;
	redraw();
}


//--
//
// EnableDetailsTexture
//
//--
void Visualizer::EnableDetailsTexture(bool enable)
{
	texture_details_enabled = enable;
	redraw();
}


//--
//
// SetBackgroundColor
//
//--
void Visualizer::SetBackgroundColor(Vector3d color)
{
	// Clamp the RGB values of the color
	color.Clamp(0.0, 1.0);
	// Change OpenGL background color
	glClearColor((float)color[0], (float)color[1], (float)color[2], 1.0);
	// Redraw the scene
	redraw();
}


//--
//
// ToggleSimplificationMode
//
//--
void Visualizer::ToggleSimplification(bool enable)
{
	simplification_enabled = enable;
	redraw();
}


//--
//
// ToggleColorBar
//
//--
void Visualizer::ToggleColorBar(bool enable)
{
	color_bar_enabled = enable;
	redraw();
}


//--
//
// ApplyTrackballTransformation
//
//--
Vector3d Visualizer::ApplyTrackballTransformation(const Vector3d& p) const
{
	// Transformed point
	Vector3d pt;

	// Apply trackball transformation matrix to point p in homogeneous coordinates
	pt[0] = trackball_transform[0][0]*p[0] + trackball_transform[0][1]*p[1] + trackball_transform[0][2]*p[2] + trackball_transform[0][3];
	pt[1] = trackball_transform[1][0]*p[0] + trackball_transform[1][1]*p[1] + trackball_transform[1][2]*p[2] + trackball_transform[1][3];
	pt[2] = trackball_transform[2][0]*p[0] + trackball_transform[2][1]*p[1] + trackball_transform[2][2]*p[2] + trackball_transform[2][3];

	// Normalize by the homogeneous coordinate to go back to 3d coordinates
	pt *= 1.0 / (trackball_transform[3][0]*p[0] + trackball_transform[3][1]*p[1] + trackball_transform[3][2]*p[2] + trackball_transform[3][3]);

	return pt;
}


//--
//
// GetViewport
//
//--
Viewport Visualizer::GetViewport()
{
	Viewport vp;

	// Field of vision
	vp.fov = view_angle * 2.0 * M_PI / 360.0; // Convert viewing angle in radians
	double fov_2 = vp.fov / 2.0;

	//
	// Note about the transformations
	//
	// Visualization:
	// general_translation + rotation
	//
	// Viewport in model space:
	// - general_translation - rotation


	// general_translation
	Vector3d camera = (Vector3d(0,0,2)-translations);

	// rotation
	vp.origin = ApplyTrackballTransformation( camera ); // Camera position
	vp.x_axis = ApplyTrackballTransformation( camera+Vector3d(sin(fov_2*aspect_ratio),0,0) ); // special point on x axis
	vp.y_axis = ApplyTrackballTransformation( camera+Vector3d(0,sin(fov_2),0) ); // special point on y axis
	vp.z_axis = ApplyTrackballTransformation( camera+Vector3d(0,0,-cos(fov_2)) ); // special point on z axis

	// Viewport axis
	vp.x_axis -= vp.origin;
	vp.z_axis -= vp.origin;
	vp.y_axis -= vp.origin;

	// Viewport corner points
	Vector3d a = vp.origin + vp.z_axis - vp.x_axis + vp.y_axis;
	Vector3d b = vp.origin + vp.z_axis + vp.x_axis + vp.y_axis;
	Vector3d c = vp.origin + vp.z_axis + vp.x_axis - vp.y_axis;
	Vector3d d = vp.origin + vp.z_axis - vp.x_axis - vp.y_axis;

	// Compute normals of the clipping plane
	vp.t_normal = ((b-vp.origin)^(a-vp.origin)).Normalize();
	vp.b_normal = ((d-vp.origin)^(c-vp.origin)).Normalize();
	vp.r_normal = ((c-vp.origin)^(b-vp.origin)).Normalize();
	vp.l_normal = ((a-vp.origin)^(d-vp.origin)).Normalize();

	// Normalize axis vectors
	vp.x_axis.Normalize();
	vp.z_axis.Normalize();
	vp.y_axis.Normalize();

	// Generate an object list of the viewpoirt pyramid
	if( glIsList(frustum_list_number) ) glDeleteLists(frustum_list_number, 1);
	frustum_list_number = glGenLists(1);
	glNewList(frustum_list_number, GL_COMPILE_AND_EXECUTE);

		// No face culling
		glDisable(GL_CULL_FACE);

		// Transparent color
		glColor4f(0.8, 0.8, 0.0, 0.2);

		// Sphere at the camera position
		glPushMatrix();
			glTranslated( vp.origin[0], vp.origin[1], vp.origin[2] );
			glCallList( sphere_list_number );
		glPopMatrix();

		// For clipping planes definind the viewport pyramid
		glPushMatrix();

			glBegin(GL_QUADS);

				// Top
				glNormal3dv(vp.t_normal);
				glVertex3dv(vp.origin);
				glVertex3dv(a);
				glVertex3dv(b);
				glVertex3dv(vp.origin);

				// Bottom
				glNormal3dv(vp.b_normal);
				glVertex3dv(vp.origin);
				glVertex3dv(c);
				glVertex3dv(d);
				glVertex3dv(vp.origin);

				// Right
				glNormal3dv(vp.r_normal);
				glVertex3dv(vp.origin);
				glVertex3dv(b);
				glVertex3dv(c);
				glVertex3dv(vp.origin);

				// Left
				glNormal3dv(vp.l_normal);
				glVertex3dv(vp.origin);
				glVertex3dv(d);
				glVertex3dv(a);
				glVertex3dv(vp.origin);

			glEnd();

			glLineWidth(5.0);

			// Transparent color
			glColor4f(0.8, 0.0, 0.8, 0.5);

			glBegin(GL_LINES);

				glVertex3dv(vp.origin);
				glVertex3dv(a);

				glVertex3dv(vp.origin);
				glVertex3dv(b);

				glVertex3dv(vp.origin);
				glVertex3dv(c);

				glVertex3dv(vp.origin);
				glVertex3dv(d);

				glVertex3dv(a);
				glVertex3dv(b);

				glVertex3dv(b);
				glVertex3dv(c);

				glVertex3dv(c);
				glVertex3dv(d);

				glVertex3dv(d);
				glVertex3dv(a);

			glEnd();

			glLineWidth(1.0);

		glPopMatrix();

		// Enable face culling
		glEnable(GL_CULL_FACE);

	glEndList();

	return vp;
}


//--
//
// DrawColorBar
//
//--
void Visualizer::DrawColorBar()
{
	Vector3d color;
	double size = 500;
	// Get viewport
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);

	double start_y = (double)viewport[3]*0.2;
	double size_y = (double)viewport[3]*0.6;

	glDisable(GL_LIGHTING);
	glMatrixMode( GL_PROJECTION );

	glPushMatrix();

		glLoadIdentity();
		gluOrtho2D(0, viewport[2], 0, viewport[3]);

		glMatrixMode( GL_MODELVIEW );

		glPushMatrix();

			glLoadIdentity();
			for( int i=0; i<size; i++ ) {
				color = Double2Color((double)i/(size-1.0));
				glBegin(GL_POLYGON);
					glColor3dv(color);
					glVertex2f(10,start_y+(float)i/size*size_y);
					glVertex2f(30,start_y+(float)i/size*size_y);
					glVertex2f(30,start_y+(float)(i+1)/size*size_y);
					glVertex2f(10,start_y+(float)(i+1)/size*size_y);
				glEnd();
			}

		glPopMatrix();

		glMatrixMode( GL_PROJECTION );

	glPopMatrix();

	glMatrixMode( GL_MODELVIEW );
	glEnable( GL_LIGHTING );
}

void Visualizer::SetTransformation(double rotx, double roty, double rotz)
{
	// Get  external rotations
	external_rotation_x = rotx;
	external_rotation_y = roty;
	external_rotation_z = rotz;
	// Redraw the scene
	redraw();
}
