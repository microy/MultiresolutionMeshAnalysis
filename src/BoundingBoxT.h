/***************************************************************************
                                BoundingBoxT.h
                             -------------------
    update               : 2003/11/09
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

#ifndef _BOUNDINGBOXT_
#define _BOUNDINGBOXT_

#include "VectorT.h"

#include <cfloat>
#include <iostream>
#include <vector>

#define _USE_OPENGL_
#ifdef _USE_OPENGL_
#include <GL/gl.h>
#endif


//--
//
// BoundingBoxT
//
//--
// Bounding box of points in nD.
// Points must have double precision float values
template<int Size>
class BoundingBoxT
{
	//--
	//
	// Type definition
	//
	//--
	protected :

		typedef VectorT<double,Size> Point;
		
	//--
	//
	// Member data
	//
	//--
	protected :

		Point minpt;
		Point maxpt;
		Point total;
		int point_number;

	//--
	//
	// Member functions
	//
	//--
	public :

		// Default construtor
		inline BoundingBoxT()
		: minpt(DBL_MAX), maxpt(DBL_MIN), total(0.0), point_number(0) {
		}

		// Copy construtor
		inline BoundingBoxT(const BoundingBoxT& bb)
		: minpt(bb.minpt), maxpt(bb.maxpt), total(bb.total), point_number(bb.point_number) {
		}

		//--
		//
		// Data management
		//
		//--
		
		// Reset data
		inline BoundingBoxT<Size>& Reset() {
			minpt = DBL_MAX;
			maxpt = DBL_MIN;
			total = 0.0;
			point_number = 0;
			return *this;
		}

		// Add a point
		inline BoundingBoxT<Size>& AddPoint( const Point& p ) {
			for( int i=0; i<Size; i++ ) {
				if( p[i] < minpt[i] ) minpt[i] = p[i];
				if( p[i] > maxpt[i] ) maxpt[i] = p[i];
			}
			total += p;
			point_number++;
			return *this;
		}

		// Add an array of points
		inline BoundingBoxT<Size>& AddPoints( const std::vector<Point>& v ) {
			for( int i=0; i<(int)v.size(); i++ ) {
				AddPoint( v[i] );
			}
			return *this;
		}
		
		//--
		//
		// Data accessors
		//
		//--

		// Minimum bounds (constant)
		inline const Point& Min() const {
			return minpt;
		}

	     	// Maximum bounds (constant)
		inline const Point& Max() const {
			return maxpt;
		}

		// Size of the bounding box (constant)
		inline Point Length() const {
			return maxpt - minpt;
		}

		// Diagonal of the bounding box (constant)
		inline double Diagonal() const {
			return Length().Length();
		}

		// Barycenter of the bounding box (constant)
		inline Point Barycenter() const {
			return total / (double)point_number;
		}

		// Center of the bounding box (constant)
		inline Point Center() const {
			Point center = minpt;
			for( int i=0; i<Size; i++ ) {
				center[i] += fabs(maxpt[i]-minpt[i])*0.5;
			}
			return center;
		}

		//--
		//
		// Operators
		//
		//--
		inline BoundingBoxT<Size>& operator=(const BoundingBoxT<Size>& bb) {
			minpt = bb.minpt;
			maxpt = bb.maxpt;
			total = bb.total;
			point_number = bb.point_number;
			return *this;
		}

		inline BoundingBoxT<Size>& operator+=(const BoundingBoxT<Size>& bb) {
			for( int i=0; i<Size; i++ ) {
				if( bb.minpt[i] < minpt[i] ) minpt[i] = bb.minpt[i];
				if( bb.maxpt[i] > maxpt[i] ) maxpt[i] = bb.maxpt[i];
			}
			total += bb.total;
			point_number += bb.point_number;
			return *this;
		}
	
		inline BoundingBoxT<Size>& operator+=(const Point& p) {
			AddPoint(p);
			return *this;
		}

		inline BoundingBoxT<Size>& operator+=(const std::vector<Point>& v) {
			AddPoints(v);
			return *this;
		}

		//--
		//
		// OpenGL display
		//
		//--
		#ifdef _USE_OPENGL_
		inline void Draw() const {
			Vector3d length = Length();
			glBegin(GL_QUADS);
				// Top
				glNormal3d(0,1,0);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]+length[2]);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]+length[2]);
				// Bottom
				glNormal3d(0,-1,0);
				glVertex3d(minpt[0], minpt[1], minpt[2]);
				glVertex3d(minpt[0], minpt[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]);
				// Right
				glNormal3d(1,0,0);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]);
				// Left
				glNormal3d(-1,0,0);
				glVertex3d(minpt[0], minpt[1], minpt[2]);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]+length[2]);
				glVertex3d(minpt[0], minpt[1], minpt[2]+length[2]);
				// Front
				glNormal3d(0,0,-1);
				glVertex3d(minpt[0], minpt[1], minpt[2]);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]);
				// Back
				glNormal3d(0,0,1);
				glVertex3d(minpt[0], minpt[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1], minpt[2]+length[2]);
				glVertex3d(minpt[0]+length[0], minpt[1]+length[1], minpt[2]+length[2]);
				glVertex3d(minpt[0], minpt[1]+length[1], minpt[2]+length[2]);
			glEnd();
		}
		#endif
};

//--
//
//Other operators
//
//--
template<int Size>
inline BoundingBoxT<Size> operator+(const BoundingBoxT<Size>& bb1, const BoundingBoxT<Size>& bb2)
{
	BoundingBoxT<Size> result(bb1);
	return result += bb2;
}

template<int Size>
inline std::ostream& operator<<(std::ostream& out, const BoundingBoxT<Size>& bb)
{
	// Output minimum and maximum points of the bounding box
	return out<<"Min: "<<bb.Min()<<" - Max: "<<bb.Max()<<" - Center: "<<bb.Center();
}

//--
//
// Type definitions
//
//--
typedef BoundingBoxT<1> BoundingBox1d;
typedef BoundingBoxT<2> BoundingBox2d;
typedef BoundingBoxT<3> BoundingBox3d;
typedef BoundingBoxT<4> BoundingBox4d;

#endif // _BOUNDINGBOXT_

