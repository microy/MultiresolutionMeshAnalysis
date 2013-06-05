/***************************************************************************
                                 Curvature.h
                             -------------------
    update               : 2003-06-11
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

#ifndef _CURVATURE_
#define _CURVATURE_

#include "NeighborMesh.h"


//--
//
// Mean Curvature
//
//--
void ComputeCurvatureNormals(NeighborMesh* mesh);
Vector3d MeanCurvatureNormal(const NeighborMesh* mesh, int v);


//--
//
// Gaussian Curvature
//
//--
double GaussianCurvature(const NeighborMesh* mesh, int v);


#endif // _CURVATURE_

