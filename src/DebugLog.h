/***************************************************************************
                                 DebugLog.h
                             -------------------
    update               : 2003/10/25
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

#ifndef _DEBUGLOG_
#define _DEBUGLOG_

#include <fstream>
#include <iostream>
#include <iomanip>

// Create a debug log file
static std::ofstream debug("debug.log");

// Use standard io stream
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

// Print out a neighbor list
inline std::ostream& operator<<(std::ostream& os, const NeighborList& l)
{
	for( ConstNeighborIterator it=l.begin(); it!=l.end(); ++it ) {
		os<<*it<<" ";
	}
	return os;
}



#endif // _DEBUGLOG_

