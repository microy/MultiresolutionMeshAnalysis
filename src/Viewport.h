/***************************************************************************
                                 Viewport.h
                             -------------------
    update               : 2003/11/08
    copyright            : (C) 2003 by Michaël Roy
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

#ifndef _VIEWPORT_
#define _VIEWPORT_

#include "VectorT.h"

//--
//
// Viewport
//
//--
class Viewport
{

	//--
	//
	// Member data
	//
	//--
	public :

		Vector3d origin;
		Vector3d x_axis;
		Vector3d y_axis;
		Vector3d z_axis;
		double fov;
		Vector3d t_normal; // Top normal
		Vector3d b_normal; // Bottom normal
		Vector3d r_normal; // Right normal
		Vector3d l_normal; // Left normal

	//--
	//
	// Member functions
	//
	//--
	public :

		// Default construtor
		inline Viewport()
		: origin(0.0), x_axis(0.0), y_axis(0.0), z_axis(0.0) {
		}

		// Copy construtor
		inline Viewport(const Viewport& vp)
		: origin(vp.origin), x_axis(vp.x_axis), y_axis(vp.y_axis), z_axis(vp.z_axis) {
		}

};

#endif // _VIEWPORT_

