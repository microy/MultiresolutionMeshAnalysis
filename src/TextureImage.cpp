/***************************************************************************
                              TextureImage.cpp
                             -------------------
    begin                : 2003/11/06
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

#include "TextureImage.h"

#include <iostream>

//--
//
// Without JPEG library
//
//--
#ifndef _HAVE_LIB_JPEG_

bool TextureImage::ReadFile(const std::string filename)
{
	std::cerr<<"This software was not compiled with JPEG library support"<<std::endl;
	return false;
}

bool TextureImage::WriteFile(const std::string filename)
{
	std::cerr<<"This software was not compiled with JPEG library support"<<std::endl;
	return false;
}



//--
//
// With JPEG library
//
//--
#else


#include <cstdio>
#include <cstdlib>
#include <iostream>


// Inhibit C++ name-mangling for libjpeg functions but not for system calls.
#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
#include <jpeglib.h>
#ifdef __cplusplus
}
#endif // __cplusplus


bool TextureImage::ReadFile(const std::string filename)
{
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[1];
	int row_stride;
	FILE * file;
	if ((file = fopen(filename.c_str(), "rb")) == NULL)
	{
	//	fprintf(stderr, "can't open %s\n", filename.c_str());
		return false;
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, file);
	jpeg_read_header(&cinfo, TRUE);
	jpeg_start_decompress(&cinfo);
	width = cinfo.output_width;
	height = cinfo.output_height;
	if( cinfo.output_components != 3 )
	{
		jpeg_destroy_decompress(&cinfo);
		fclose(file);
		return false;
	}
	row_stride = cinfo.output_width * cinfo.output_components;
	if( buffer != 0 )
	{
		free(buffer);
	}
	buffer = (unsigned char*)malloc( width * height * 3 );
	while( cinfo.output_scanline < cinfo.output_height )
	{
		row_pointer[0] = &buffer[(cinfo.output_height-cinfo.output_scanline-1) * row_stride];
		jpeg_read_scanlines(&cinfo, row_pointer, 1);
	}
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(file);

	return true;
}

bool TextureImage::WriteFile(const std::string filename)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	int quality = 100;
	FILE* outfile; // target file
	JSAMPROW row_pointer[1]; // pointer to JSAMPLE row[s]
	int row_stride; // physical row width in image buffer
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	if ((outfile = fopen(filename.c_str(), "wb")) == NULL)
	{
	//	fprintf(stderr, "can't open %s\n", filename.c_str());
		return false;
	}
	jpeg_stdio_dest(&cinfo, outfile);
	cinfo.image_width = width; // image width and height, in pixels
	cinfo.image_height = height;
	cinfo.input_components = 3; // # of color components per pixel
	cinfo.in_color_space = JCS_RGB; // colorspace of input image
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);
	jpeg_start_compress(&cinfo, TRUE);
	row_stride = width * 3; // JSAMPLEs per row in image_buffer
	while (cinfo.next_scanline < cinfo.image_height)
	{
		row_pointer[0] = &buffer[(cinfo.image_height-cinfo.next_scanline-1) * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);

	return true;
}

#endif
