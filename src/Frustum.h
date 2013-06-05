/***************************************************************************
                                  Frustum.h
                             -------------------
    update               : 2003/11/09
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

/***************************************************************************
 *                                                                         *
 *   This code was higly inspired by DigiBen                               *
 *   Programming Lesson on Frustum Culling                                 *
 *   digiben@gametutorials.com                                             *
 *   www.gametutorials.com                                                 *
 *                                                                         *
 ***************************************************************************/


#ifndef _FRUSTUM_
#define _FRUSTUM_

#include "VectorT.h"

//--
//
// Frustum
//
//--
// This will allow us to create an object to keep track of our frustum
class Frustum
{

	//--
	//
	// Member data
	//
	//--
	protected:

		// This holds the A B C and D values for each side of our frustum.
		double frustum[6][4];

	//--
	//
	// Member functions
	//
	//--
	public:

		// Call this every time the camera moves to update the frustum
		void CalculateFrustum();

		// This takes a 3D point and returns TRUE if it's inside of the frustum
		bool PointInFrustum(const Vector3d& p);

	protected :

		// This normalizes a plane (A side)
		void NormalizePlane(int side);
/*
		// This takes a 3D point and a radius and returns TRUE if the sphere is inside of the frustum
		bool SphereInFrustum(float x, float y, float z, float radius);

		// This takes the center and half the length of the cube.
		bool CubeInFrustum( float x, float y, float z, float size );
*/
};

#endif // _FRUSTUM_

