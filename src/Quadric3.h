/***************************************************************************
                                  Quadric3.h
                             -------------------
    update               : 2003/11/06
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

#ifndef _QUADRIC3_
#define _QUADRIC3_

#include "VectorT.h"


//--
//
// Quadric3
//
//--
class Quadric3
{

	//--
	//
	// Member variables
	//
	//--
	protected :

		// Ten unique coefficients
		double aa, ab, ac, b1;
		double     bb, bc, b2;
		double         cc, b3;
		double              c;

		double area;


	//--
	//
	// Member functions
	//
	//--
	public :


		// Default Constructor
		Quadric3() :
			aa(0), ab(0), ac(0), b1(0),
			       bb(0), bc(0), b2(0),
			              cc(0), b3(0),
			                      c(0), area(1) {
		}

		// Copy contructor
		Quadric3(const Quadric3& q) :
			aa(q.aa), ab(q.ab), ac(q.ac), b1(q.b1),
			          bb(q.bb), bc(q.bc), b2(q.b2),
			                    cc(q.cc), b3(q.b3),
			                               c(q.c), area(1) {
		}

		// Constructor from a plane
		Quadric3(const Vector3d& n, const double& s) :
			aa(n[0]*n[0]), ab(n[0]*n[1]), ac(n[0]*n[2]), b1(s*n[0]),
			               bb(n[1]*n[1]), bc(n[1]*n[2]), b2(s*n[1]),
			                              cc(n[2]*n[2]), b3(s*n[2]),
			                                              c(s*s), area(1) {
		}

		// Constructor from values
		Quadric3(double aa, double ab, double ac, double b1, double bb, double bc, double b2, double cc, double b3, double c) :
			aa(aa), ab(ab), ac(ac), b1(b1),
			        bb(bb), bc(bc), b2(b2),
			                cc(cc), b3(b3),
			                        c(c), area(1) {
		}

		//
		// Point evaluations
		//
		double Evaluate( const double& x, const double& y, const double& z ) const {
			return x*x*aa + 2.0*x*y*ab + 2.0*x*z*ac + y*y*bb + 2.0*y*z*bc + z*z*cc + 2.0*x*b1 + 2.0*y*b2 + 2.0*z*b3 + c;
		}

		double Evaluate( const Vector3d& v ) const {
			return Evaluate( v[0], v[1], v[2] );
		}

		//
		// Operators
		//
		Quadric3& operator+=(const Quadric3& q) {
			aa += q.aa; ab += q.ab; ac += q.ac; b1 += q.b1;
			            bb += q.bb; bc += q.bc; b2 += q.b2;
			                        cc += q.cc; b3 += q.b3;
			                                    c  += q.c;
			return *this;
		}

		Quadric3& operator*=(const double& x) {
			aa *= x; ab *= x; ac *= x; b1 *= x;
			         bb *= x; bc *= x; b2 *= x;
			                  cc *= x; b3 *= x;
			                           c  *= x;
			return *this;
		}
		
		//
		// Area accessors
		//
		double& Area() {
			return area;
		}
		
		const double& Area() const {
			return area;
		}

		Vector3d Vector() const { return Vector3d(b1,b2,b3); }
	
		double Offset() const { 
			return c;
		}

};

#endif // _QUADRIC3_

