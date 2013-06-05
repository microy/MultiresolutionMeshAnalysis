/***************************************************************************
                             ResolutionLevel.cpp
                             -------------------
    update               : 2003-06-11
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

#include "ResolutionLevel.h"

//
void ResolutionLevel::Resize( int size )
{
	vertex_types.assign( size, UNDEFINED_VERTEX );
	subdivision_weights.resize( size );
}
