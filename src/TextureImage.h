/***************************************************************************
                               TextureImage.h
                             -------------------
    begin                : 2009/02/19
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

#ifndef _TEXTUREIMAGE_
#define _TEXTUREIMAGE_

#define _HAVE_LIB_JPEG_

#include <cassert>
#include <string>
#include <cstdlib>

//--
//
// Texture
//
//--
// Handle texture image input/output
class TextureImage
{
	public :

		inline TextureImage() : buffer(0) {}
		inline ~TextureImage() { if(buffer) free(buffer); }

		inline void Reset() {
			if( buffer != 0 ) { free(buffer); buffer = 0; }
			height = 0;
			width = 0;
		}

		bool ReadFile(const std::string filename);
		bool WriteFile(const std::string filename);

		inline unsigned char* Buffer() {
			return buffer;
		}

		inline const unsigned char* Buffer() const {
			return buffer;
		}

		inline unsigned char& Buffer(int i) {
			assert( buffer );
			assert( (i>=0) && (i<(width*height*3)) );
			return buffer[i];
		}

		inline const unsigned char& Buffer(int i) const {
			assert( buffer );
			assert( (i>=0) && (i<(width*height*3)) );
			return buffer[i];
		}

		inline int Width() const {
			return width;
		}

		inline int Height() const {
			return height;
		}

	protected :

		unsigned char* buffer; // Points to large array of R,G,B-order data
		int height; //Number of rows in image
		int width; // Number of columns in image
};

#endif // _TEXTUREIMAGE_

