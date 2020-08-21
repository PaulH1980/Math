#ifndef MATH_DEFS_H
#define MATH_DEFS_H
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


namespace Enigma
{
#define									INVALID_INDEX				( -1												)
#define									MATH_PI						( 3.1415926535897932384626433832795					)
#define									MATH_EULER					( 2.7182818284590452353602874713527					)		
#define									HALF_PI						( MATH_PI * 0.5										)
#define									QUARTER_PI					( MATH_PI * 0.25									)
#define									TWO_PI						( MATH_PI * 2.0										)
#define									RADIANS_TO_DEGREES			( 180.0 / MATH_PI									)
#define									DEGREES_TO_RADIANS			( MATH_PI / 180.0									)
#define									INT32MIN					(-2147483647 - 1									) 
#define									INT32MAX					( 2147483647										)
#define									UINT32MIN					( 0													)
#define									UINT32MAX					( 4294967295										)
#define									LARGE_FLOAT					( 10000000.0f										)
#define									FLOAT_EPSILON				( 0.00125f											)
#define									MATRIX_INVERT_EPSILON		( 0.0000000000001									)
#define									PLANE_EPSILON				( 0.00125f											)
#define									VEC_COMPARE_EPS				( 0.00125f											)
#define									QUAT_COMPARE_EPS			( 0.00125f											)
#define									INFINITE_FLOAT				( FLT_MAX											)
#define									M_MIN_NEARCLIP				( 0.01f												)
	
}

#endif