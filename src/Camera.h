/***************************************************************************
                                  Camera.h
                             -------------------
    update               : 2003/11/13
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

#ifndef _CAMERA_
#define _CAMERA_


//--
//
// Camera
//
//--
class Camera
{

	//--
	//
	// Member data
	//
	//--
	protected:

		// The camera's position
		Vector3d eye;

		// The camera's view
		Vector3d center;

	//--
	//
	// Member functions
	//
	//--
	public:

		Camera(double ex=0, double ey=0, double ez=1, double cx=0, double cy=0, double cz=0)
		: eye(ex,ey,ez), center(cx, cy, cz) {}

		Vector3d& Eye() { return eye; }
		const Vector3d& Eye() const { return eye; }
		Vector3d& Center() { return center; }
		const Vector3d& Center() const { return center; }

};

#endif
