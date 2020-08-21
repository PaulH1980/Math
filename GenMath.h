#ifndef GEN_MATH_H
#define GEN_MATH_H
#include <array>
#include <vector>
#include <random>
#include <limits>
#include "MathDecl.h"
#include "CameraRay.h"
#include "MathTypes.h"

namespace Misc
{
	template<typename T>
	static inline T							ConvertTo(const char* dataPtr, u64 offset)
	{
		const T* val = reinterpret_cast<const T*>(dataPtr + offset);
		return *val;
	}
}

namespace Math
{
    
    const auto		FLOAT_MIN		= std::numeric_limits<float>::lowest();
    const auto		FLOAT_MAX		= std::numeric_limits<float>::max();

    const auto		DOUBLE_MIN		= std::numeric_limits<double>::lowest();
    const auto		DOUBLE_MAX		= std::numeric_limits<double>::max();
    
    
    template <class T>
    inline T  MaxVal( const T& a, const T& b )
    {
        if( a > b ) 
            return a;
        return b;
    }

    template <class T>
    inline T  MinVal( const T& a, const T& b )
    {
        if( a < b ) 
            return a;
        return b;
    }

    template <class T>
    inline T  MaxVal( const T& a, const T& b, const T& c )
    {
        return MaxVal( MaxVal(a, b ), c );
    }

    template <class T> 
    inline T GCD( T x, T y )
    {
        static_assert(std::is_integral<T>::value, "Not An Integral");
        if (x < y)
            std::swap(x, y);

        while (y > 0)
        {
            auto f = x % y;
            x = y;
            y = f;
        }
        return x;
    }

    template <class T>
    inline T	 Log( T val, T base = T(2.0))
    {
        return static_cast<T>(log( val ) / log(base));
    }

    template<class T>
    inline constexpr T   Pi()
    {
        return static_cast<T>(3.14159265359);
    }

    template<class T>
    inline constexpr T   HalfPi()
    {
        return static_cast<T>(3.14159265359 / 2.0 );
    }

    template<class T>
    inline constexpr T   TwoPi()
    {
        return static_cast<T>(3.14159265359 * 2.0);
    }



    template <class T>
    inline int Round( T val )
    {
        if( val > T(0) )
            return int( val + T(0.5) );
        return int( val - T(0.5) );
    }	

    template <class T>
    inline T SnapToGrid( const T& inputVal, const T& gridVal )
    {
        int val =  Math::Round( inputVal / gridVal  );
        return val * gridVal;
    }
        
    template <class T>
    inline bool  IsPowerOfTwo( T val )
    {
        static_assert(std::is_integral<T>::value, "Not An Integral");
        return ( ( val & ( val - T(1)) ) == T(0));
    }

    template <class T>
    inline bool  IsOdd( T val )
    {
        static_assert(std::is_integral<T>::value, "Not An Integral");
        return val & T(1);
    }

    template<typename T>
    inline int SolveQuadratic( const T &a, const T &b, const T &c, T &x0, T &x1)
    {
        T discr = b * b - 4 * a * c;
        if (discr < T(0)) 
            return 0; 
        else if (discr == T(0))
        {
            x0 = x1 = T(-0.5) * b / a;
            return 1; //1st solution
        }
        else 
        {
            T q = b > T(0) ? T(-0.5) * (b + sqrt(discr)) :
                    T(-0.5)  * (b - sqrt(discr));
            x0 = q / a;
            x1 = c / q;
        }
        if (x0 > x1) 
            std::swap(x0, x1);
        return 2; //2nd solution
    }


    template<typename T, typename U>
    inline T EvaluateQuadricPatch(const T& v1, const T& v2, const T& v3, const U t )
    {
        auto b    = static_cast<U>(1.0) - t;
		auto b2	  = b * b;
		auto twoB = static_cast<U>(2.0) * b;
		auto t2	  = t * t;
				
        auto result = v1 * (b2) +
					  v2 * (twoB * t) +
					  v3 * (t2);
        return result;
    }
	
	template<typename T, typename U>
	inline T EvaluateCubicPatch(const T& v1, const T& v2, const T& v3, const T& v4, const U t)
	{
		auto b = static_cast<U>(1.0) - t;
		auto b2 = b * b;
		auto b3 = b2 * b;
		
		auto t2 = t * t;
		auto t3 = t2 * t;

		auto result = ( v1 * b ) + 
					  ( v2 * 3 * b2 * t ) +
					  ( v3 * 3 * b * t2 ) + 
					  ( v4 * t3 );
		return result;
	}



    template<typename T>
    inline T Abs( T val )
    {
        if( val < T(0) )
            return -val;
        return val;
    }


    template<typename T>
    inline int Sign( T& val )
    {
        if( val < T(0) )
            return -1;
        return 1;
    }

    inline int NextPowerOfTwo( int input )
    {
        input--;
        input |= input >> 0x01;
        input |= input >> 0x02;
        input |= input >> 0x04;
        input |= input >> 0x08;
        input |= input >> 0x10;
        input++;
        return input;
    }
        
    inline  int PreviousPowerOfTwo( int val )
    {
        return NextPowerOfTwo( val ) >> 1;
    }


    inline int NearestPowerOfTwo( int val )
    {
        auto npt = NextPowerOfTwo( val );
        auto pvt = PreviousPowerOfTwo( val );
        if( Abs( npt - val ) > Abs( pvt - val ) ) 
            return pvt;
        return npt;
    }

    inline unsigned HighestBit( unsigned input )
    {
        input |= (input >>  1);
        input |= (input >>  2);
        input |= (input >>  4);
        input |= (input >>  8);
        input |= (input >> 16);
        return input & ~(input >> 1);
    }

	template<typename T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
	inline T		RandomValue(const T& min, const T& max)
	{
		std::random_device rd;
		std::uniform_int_distribution<T> dist(min, max);
		return static_cast<T>(dist(rd));
	}


	template<typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
	inline T		RandomValue(const T& min, const T& max)
	{
		std::random_device rd;
		std::uniform_real_distribution<T> dist(min, max);
		return static_cast<T>(dist(rd));
	}


    //////////////////////////////////////////////////////////////////////////
    // \Brief: Gets the center of gravity for some points
    //////////////////////////////////////////////////////////////////////////
    template<class Real, template<class> class Dim>
    inline Dim<Real> GetCenter( const  Dim<Real>* points, int numPoints )
    {
         Dim<Real> center;
         for( int i = 0; i < numPoints; i++ )
             center += points[i];
         center /= Real(numPoints);
         return center;
    }
    

    template<typename T>
    T SnapScalarToGrid(const T& val, const T& gridSize, bool floor = true)
    {
        auto roundedVal = val / gridSize;
        roundedVal = floor ? std::floor(roundedVal) : std::ceil(roundedVal);
        return roundedVal * gridSize;
    }
    template<typename T>
    Math::Vector3<T> SnapPointToGrid(const Math::Vector3<T>& point, const T& gridSize, bool floor = true)
    {
        auto newX = SnapScalarToGrid(point[0], gridSize, floor);
        auto newY = SnapScalarToGrid(point[1], gridSize, floor);
        auto newZ = SnapScalarToGrid(point[2], gridSize, floor);
        return Math::Vector3<T>(newX, newY, newZ);
    }
    


    //////////////////////////////////////////////////////////////////////////
    // \Brief convert a rgb to hsv value
    //////////////////////////////////////////////////////////////////////////
    inline Vector3f RgbToHsv( const Vector3f& rgb )
    {
        auto r = rgb[0];
        auto g = rgb[1];
        auto b = rgb[2];
        auto K = 0.f;

        if (g < b)
        {
            std::swap(g, b);
            K = -1.f;
        }

        if (r < g)
        {
            std::swap(r, g);
            K = -2.f / 6.f - K;
        }

        auto chroma = r - MinVal(g, b);
        auto h = fabs(K + (g - b) / (6.f * chroma + 1e-20f));
        auto s = chroma / (r + 1e-20f);
        auto v = r;
        return Vector3f( 1.0f - h, s, v );
    }



    //////////////////////////////////////////////////////////////////////////
    // \Brief: Convert a hsv value to a rgb value
    //////////////////////////////////////////////////////////////////////////
    inline Vector3f HsvToRgb( const Vector3f& hsv )
    {
        int i;
        float f, p, q, t;
        auto v = hsv[2];
        auto h = hsv[0] >= 360.0f ? 0.0f : hsv[0]/60.0f;			// sector 0 to 5
        auto s = hsv[1];
        if( s == 0.0f )
            return Vector3f( v, v, v );

        Vector3f rgba;
        
        i = static_cast<int>( floor( h ) );
        f = h - i;			
        p = v * ( 1 - s );
        q = v * ( 1 - s * f );
        t = v * ( 1 - s * ( 1 - f ) );

        switch( i ) {
        case 0:
            rgba[0] = v;
            rgba[1] = t;
            rgba[2] = p;
            break;
        case 1:
            rgba[0] = q;
            rgba[1] = v;
            rgba[2] = p;
            break;
        case 2:
            rgba[0] = p;
            rgba[1] = v;
            rgba[2] = t;
            break;
        case 3:
            rgba[0] = p;
            rgba[1] = q;
            rgba[2] = v;
            break;
        case 4:
            rgba[0] = t;
            rgba[1] = p;
            rgba[2] = v;
            break;
        default:		
            rgba[0] = v;
            rgba[1] = p;
            rgba[2] = q;
            break;
        }
        return rgba;
    }

    

    //////////////////////////////////////////////////////////////////////////
    // \brief: Returns the number of set bits in this byte
    //////////////////////////////////////////////////////////////////////////
    inline int PopCount8( unsigned char input )
    {
        const static int bitCount[] =
        {
            0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
            3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
            4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
        };
        return bitCount[input];
    }

    //////////////////////////////////////////////////////////////////////////
  // \brief: Returns the number of set bits in this byte
  //////////////////////////////////////////////////////////////////////////
    inline int PopCount32(std::uint32_t input)
    {
        int result = 0;
        result += PopCount8(static_cast<std::uint8_t>(input & 0x000000FF));
        result += PopCount8(static_cast<std::uint8_t>(input & 0x0000FF00));
        result += PopCount8(static_cast<std::uint8_t>(input & 0x00FF0000));
        result += PopCount8(static_cast<std::uint8_t>(input & 0xFF000000));
        return result;
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief: Returns the number of set bits before 'bitIdx'
    //////////////////////////////////////////////////////////////////////////
    inline int BitCountBefore(int bitIdx, unsigned char input )
    {
        std::uint8_t bit = 1 << bitIdx;
        bit += bit - 1; //fill lower bits with '1's
        input &= bit; //cap of most significant bits
        return PopCount8( input );
    }
        

    template<class T>
    inline T Random()
    {
        return (T)( rand() & 0x7fff ) / ( (float) 0x7fff );
    }
    

    template<class T>
    inline T Random( T min, T max )
    {
        T range = T(max) - T(min);
        float tVal	= (float) ( rand() & 0x7fff ) / ( (float) 0x7fff );
        return T( range * tVal + min );
    }

	//////////////////////////////////////////////////////////////////////////
	// \Brief: Calculate average of a range
	//////////////////////////////////////////////////////////////////////////
	template<typename T>
	inline T	CalcAverage(const T* begin, const T* end, size_t stride = 0)
	{
		const std::uint8_t* dataPtr = reinterpret_cast<const std::uint8_t*>(begin);
		T avg = (T) 0.0;
		size_t count = 0;
		size_t offset = stride == 0 ? sizeof(T) : stride;
		do
		{
			avg += *(reinterpret_cast<const T*>(dataPtr));
			dataPtr += offset;
			count++;
		} while (dataPtr != reinterpret_cast<const std::uint8_t*>(end));
		return avg / (T)count;
	}


 
    ///////////////////////////////////////////////////////////////////////////
    // \brief: Make Rectangle 
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline void			MakeRect(const Math::Vector2<T>& minIn, const Math::Vector2<T>& maxIn, 
        const T viewportHeight, std::array<Math::Vector3<T>, 4>& rectangle)
    {
        rectangle[0].setXYZ(static_cast<T>(minIn[0]), static_cast<T>(viewportHeight - minIn[1] ), static_cast<T>(0.0));
        rectangle[1].setXYZ(static_cast<T>(minIn[0]), static_cast<T>(viewportHeight - maxIn[1] ), static_cast<T>(0.0));
        rectangle[2].setXYZ(static_cast<T>(maxIn[0]), static_cast<T>(viewportHeight - maxIn[1] ), static_cast<T>(0.0));
        rectangle[3].setXYZ(static_cast<T>(maxIn[0]), static_cast<T>(viewportHeight - minIn[1] ), static_cast<T>(0.0));
    }

    ///////////////////////////////////////////////////////////////////////////
    // \brief: Make Rectangle 
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline void								QuadToTriangles(const std::array<Math::Vector3<T>, 4>& rectangle, std::array<Math::Vector3<T>, 6>& triangles)
    {
        static const int QUAD_TO_TRIANGLE[6] = { 0, 1, 2, 0, 2, 3 };
        for (int i = 0; i < 6; ++i) {
            auto idx = QUAD_TO_TRIANGLE[i];
            triangles[i] = rectangle[idx];
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    // \brief: Make Rectangle 
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline void								MakeRect(int left, int right, int top, int bottom, std::array<Math::Vector3<T>, 4>& rectangle )
    {
      
		const T one = static_cast<T>(0.0);
		
		rectangle[0].setXYZ((T)left,  (T)top,		one );	//top left
        rectangle[1].setXYZ((T)left,  (T)bottom,	one );	//bottom left
        rectangle[2].setXYZ((T)right, (T)bottom,	one ); //bottom right
        rectangle[3].setXYZ((T)right, (T)top,		one );  //top right
    }


    template<class T>
    inline bool IsNaN( T val )
    {
        return val != val;
    }

    template<class T>
    inline T Sin( T angleRadians )
    {
        return  std::sin( angleRadians );
    }

    template<class T>
    inline T Cos( T angleRadians )
    {
        return  Sin( angleRadians + (T)HALF_PI );
    }


    template <class T> 
    inline T Atan2ToRads(T y, T x)
    {
        auto val = std::atan2(y, x);
        if (val < static_cast<T>(0.0))
            val += static_cast<T>(TWO_PI);
        return val;        
    }


    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    Vector2<T> SignNotZero(const Vector2<T> &v) {
        return Vector2<T>( (v.getX() >= (T)0.0) ? (T)1.0 : (T)-1.0,
                           (v.getY() >= (T)0.0) ? (T)1.0 : (T)-1.0);
    }
    
    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////

    template<class T>
    Vector2<T> float32x3_to_oct(const Vector3<T>& v) 	// Assume normalized input. Output is on [-1, 1] for each component.
    {
        // Project the sphere onto the octahedron, and then onto the xy plane
        Vector2<T> p = v.toVector2();
        p *= (T)1.0 / (Abs(v[0]) + Abs(v[1]) + Abs(v[2]) );
        // Reflect the folds of the lower hemisphere over the diagonals
        if (v.getZ() <= (T)0.0)
            p = ( Vector2<T>((T)1.0) - Vector2<T>( Abs(p[1]), Abs(p[0]) ) )* SignNotZero(p);	
        
        return p;
    }
    
    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    Vector3<T> oct_to_float32x3(const Vector2<T>& e) 
    {
        Vector3<T> v( e[0], e[1], (T)1.0 - Abs(e[0]) - Abs(e[1]) );
        if (v[2] < (T)0.0) 
        {
            Vector2<T> xy = (Vector2<T>((T)1.0) - Vector2<T>(Abs(v[1]), Abs(v[0]))) * SignNotZero(Vector2<T>(e[0], e[1]));
            v[0] = xy[0];
            v[1] = xy[1];			
        }
        return v.getNormalized();
    }

    template<class T>
    inline T Clamp(T minVal, T maxVal, T inputVal)
    {
        if (inputVal < minVal)
            return minVal;
        else if (inputVal > maxVal)
            return maxVal;
        return inputVal;
    }

   

    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    inline void quantizeOct(const Math::Vector2f& octData, int numBits, T& uOut, T& vOut )
    {
        auto scale = float((1 << numBits) - 1);
        auto norm = ((octData * 0.5f) + 0.5f) * scale;
        norm[0] = Clamp(0.0f, scale, norm[0]);
        norm[1] = Clamp(0.0f, scale, norm[1]);
        uOut = (T)norm[0];
        vOut = (T)norm[1];
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    inline void quantizeOct(const Math::Vector2f& octData, T& uOut, T& vOut)
    {
        quantizeOct( octData, static_cast<int>(sizeof(T)), uOut, vOut );
    }


    inline Vector4ui_2_10_10_10 NormalToOctahedron(const Math::Vector3f& normal)
    {
        Vector4ui_2_10_10_10 result;        
        const auto octData = Math::float32x3_to_oct(normal);
        std::uint32_t u, v;
        quantizeOct(octData, 10, u, v);
        result.m_x = u;
        result.m_y = v;
        return result;
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    Vector2f normalizeOct( T uIn, T vIn, int numBits )
    {
        const auto scale = float((1 << numBits) - 1);
        return (Vector2f(uIn / scale, vIn / scale) * 2.0f) - 1.0f;
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief: See paper octahedron normals
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    Vector2f normalizeOct( T uIn, T vIn )
    {
        int numBits = sizeof(T) * 8;
        float scale = float((1 << numBits) - 1);
        return (Vector2f(uIn / scale, vIn / scale) * 2.0f) - 1.0f;
    }

    //////////////////////////////////////////////////////////////////////////
  // \Brief: See paper octahedron normals
  //////////////////////////////////////////////////////////////////////////

    inline Vector3f         OctahedronToNormal(const Vector4ui_2_10_10_10& normalEnc)
    {
        const auto norm = normalizeOct(normalEnc.m_x, normalEnc.m_y, 10);
        return oct_to_float32x3(norm);
    }

    

   

    template<class T>
    inline void SinCos( T angleRadians, T& sinOut, T& cosOut )
    {
        sinOut = sin( angleRadians );
        cosOut = cos( angleRadians );
    }
	
	template<typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    inline T ToRadians( T degrees )
    {
        return degrees * static_cast<T>( DEGREES_TO_RADIANS );
    }

	template<typename T, typename std::enable_if<std::is_floating_point<T>::value>::type* = nullptr>
    inline T ToDegrees( T radians )
    {
        return radians * static_cast<T>(RADIANS_TO_DEGREES);
    }
    
    template<class T>
    inline T NormalizeRange( const T& value, const T& minValue, const T& maxValue)
    {
        if (value < minValue)
            return static_cast<T>(0.0);
        else if (value > maxValue)
            return static_cast<T>(1.0);
        auto range = maxValue - minValue;
        auto offset = value - minValue;
        return offset / range;
    }


    template <typename T, typename U>
    inline U NormalizeRange(T value)
    {
        static const U oneOverMin = static_cast<U>(1.0) / std::numeric_limits<T>::min();
        static const U oneOverMax = static_cast<U>(1.0) / std::numeric_limits<T>::max();

        return value < 0
            ? -static_cast<U>(value) * oneOverMin
            :  static_cast<U>(value) * oneOverMax;
    }

    template <typename T, typename U>
    inline U ConvertRange(T value)
    {
        return value < 0
            ? -static_cast<U> (value *  std::numeric_limits<U>::min())
            :  static_cast<U> (value *  std::numeric_limits<U>::max());
    }


        
    template<class T>
    inline T Lerp( T start, T end, float time )
    {
        if( time <= 0.0f )
            return start;
        else if( time >= 1.0f )
            return end;
        return  start + (end - start) * (T)time ;
    }


	template <class T>
	inline T Lerp(const std::vector<T>& input, const float t)
	{
		int idxMin = 0;
		int idxMax = static_cast<int>(input.size()) - 1;

		if (t <= 0.0f)
			return input[idxMin];
		else if (t >= 1.0f)
			return input[idxMax];

		float idxVal = t * idxMax;
		idxMin = static_cast<int>(std::floor(idxVal));
		idxMax = static_cast<int>(std::ceil(idxVal));
		float fracPart = idxVal - static_cast<float>(idxMin);
		return Lerp(input[idxMin], input[idxMax], fracPart);
	}

    ///*
    //    
    //*/
    //template <class T>
    //inline T NonLerp( const std::vector<T>& input, const std::vector<float>& values, const float t) {
    //    assert(!input.empty() && input.size() == values.size());

    //    const auto idxMin = 0;
    //    const auto idxMax = static_cast<int>(input.size()) - 1;

    //    if (t <= values[idxMin])
    //        return input[idxMin];
    //    else if ( t >= values[idxMax] )
    //        return input[idxMax];
    //    for (int i = 0; i < idxMax; ++i) {
    //        const auto nxt = i + 1;
    //        const auto firstVal = values[i];
    //        const auto nextVal = values[nxt];
    //        if( t>= firstVal && t <=nextVal)

    //        
    //    }


    //}



    template<class T>
    inline Vector3<T> CalculateNormal( const Math::Vector3<T>& v1,  const Math::Vector3<T>& v2, const Math::Vector3<T>& v3)
    {
        auto v1v3 = v3 - v1;
        auto v1v2 = v2 - v1;
        auto normal = v1v2.cross( v1v3 );
        normal.normalize();
        return normal;
    }


    

    //////////////////////////////////////////////////////////////////////////
    // \Brief: Compresses a normal to an std::uint32_t
    //////////////////////////////////////////////////////////////////////////
    inline void		CompressNormal( const Math::Vector3f& normal, std::uint32_t& resultOut )
    {
        std::uint32_t x = (int) ( 1023.0f * ( 0.5f * ( normal[0] + 1.0f) ) );
        std::uint32_t y = (int) ( 1023.0f * ( 0.5f * ( normal[1] + 1.0f) ) );
        std::uint32_t z = (int) ( 1023.0f * ( 0.5f * ( normal[2] + 1.0f) ) );
        resultOut = x | (y << 10) | (z << 20);		
    }


    template <class T>
    inline bool		Equals( const T& left, const T& right, const T& epsilon )
    {
        return Math::Abs( left - right ) <= epsilon;
    }


    template <class T>
    inline bool	PointInPoly(const Vector3<T>& point, const Vector3<T>* triangleVerts, int numVerts, size_t stride = 0 )
    {
        
        if (stride == 0)
            stride = sizeof(Vector3<T>);

        const char* dataPtr = reinterpret_cast<const char*>(triangleVerts);
        Math::Plane<T> polyPlane(
            *reinterpret_cast<const Vector3<T>*>(dataPtr + stride * 0),
            *reinterpret_cast<const Vector3<T>*>(dataPtr + stride * 1),
            *reinterpret_cast<const Vector3<T>*>(dataPtr + stride * 2));

        if (std::fabs(polyPlane.distanceToPlane(point)) > PLANE_EPSILON)
            return false;
		
        //loop through edges point should be in the back of all edges
        for( int i = 0; i < numVerts; i++ )
        {
            //int first, second,  next;
            int first = i;
            int second = (i+1) % numVerts;

            size_t offsetFirst  = first * stride;
            size_t offsetSecond = second * stride;
            

            Math::Vector3<T> p1 = *(reinterpret_cast<const Vector3<T>*>( dataPtr + offsetFirst  ));
            Math::Vector3<T> p2 = *(reinterpret_cast<const Vector3<T>*>( dataPtr + offsetSecond ));
        
            Math::Vector3<T> polyEdge	 = p2 - p1;
            Math::Vector3<T> pointToEdge = point -  p1;
            Math::Vector3<T> perp = pointToEdge.cross( polyEdge );

            T side = polyPlane.getNormal().dot( perp );
            if (side < (T) 0.0) 
                return false; //  outside
        }
        return true;
    }

    //////////////////////////////////////////////////////////////////////////
    //\Brief: Return if the input vertices are sorted in a CCW manner
    //////////////////////////////////////////////////////////////////////////
    template <class T>
    bool VerticesCounterClockWise( const std::vector<Vector3<T>>& points, const Vector3<T>&  normal )
    {
#if DEBUG
        if (!points.size())
            throw std::exception("No data");
#endif
        
        Math::Vector3<T> center;
        for( int i = 0; i < points.size(); ++i )
            center += points[i];

        center /= static_cast<T>( points.size() );
            
        
        bool cw = false;
        for( int i = 0; i < points.size() ; ++i )
        {
            int next = (i + 1) % points.size();

            auto p1 = points[i];
            auto p2 = points[next];

            auto cProd = (p1 - center).cross( p2 - center );
            auto dProd = normal.dot( cProd );
            if( dProd < (T)0.0 )
            {
                cw = true;
                break;
            }
        }	
        return !cw;
    }


    template<class T> 
    Vector3<T> GetPolyCenter( const Vector3<T>* points, int numPoints, size_t stride = 0 )
    {
#if DEBUG
        if (!numPoints)
            throw std::exception("No data");
#endif

        Vector3<T> res;
        stride = stride != 0 ? stride : sizeof(Vector3<T>);
        T oneOverPoints = static_cast<T>(1.0) / numPoints;
        auto dataPtr = reinterpret_cast<const char*>(points);
        for (int i = 0; i < numPoints; ++i) 
        {
            res += *(reinterpret_cast< const Vector3<T>*>(dataPtr));
            dataPtr += stride;
        }
        res *= oneOverPoints;
        return res;
    }

    template<class T>
    inline bool   AABBInsideFrustum(const Math::BBox3<T>& bounds, const Math::Plane<T>* planeList, int numPlanes)
    {
        auto center = bounds.getCenter();
        auto edge = center - bounds.getMin();

        for (int i = 0; i < numPlanes; ++i)
        {
            const auto& plane = planeList[i];
            T dist = plane.getNormal().dot(center) + plane.getDistance();
            auto planeAbs = plane.getNormal();

            planeAbs.absolute();
            T absDist = planeAbs.dot(edge);
            //all behind this plane, so are the others
            if (dist < -absDist)
                return false;
        }
        return true;
    }



    template<class T>
    std::vector<int> SortVerticesCCW( const Vector3<T>* points, int numPoints, const Vector3<T>&  normal )
    {
        struct AngleIndices
        {
            AngleIndices(){}
            int   idx;
            T     angle;
            bool  operator <= ( const AngleIndices& rhs ) const {
                return angle <= rhs.angle;
            }

            bool  operator < ( const AngleIndices& rhs ) const {
                return angle < rhs.angle;
            }
        };

        Vector3<T> center;
        for( int i = 0; i < numPoints; ++i )
            center += points[i];

        center /= static_cast<T>( numPoints );

        
        std::vector<AngleIndices> angles;
        //calculate angle for each vertex				
        for( int i = 0; i < numPoints; i += 1 )
        {
            auto p1 = points[i] - center;
            auto p2 = points[0] - center;

            p1.normalize();
            p2.normalize();

            auto angle = p1.angleBetween( p2 );

            auto cProd = (p1 ).cross( p2 );
            auto dProd = normal.dot( cProd );
            if( dProd > (T)0.0 )
                angle = -angle;
            if( angle < (T)0.0 )
                angle = (T)TWO_PI + angle;


            AngleIndices angleI;
            angleI.angle = Math::ToDegrees( angle  );
            angleI.idx = i;
            angles.push_back( angleI ) ;			
        }
        std::sort( angles.begin(), angles.end() );
        
        std::vector<int> newIndices;   	
        for( unsigned i = 0; i < angles.size(); i++ )
            newIndices.push_back( angles[i].idx );			

        return newIndices;
    }




    //////////////////////////////////////////////////////////////////////////
    //\Brief: Returns the barycentric coordinate for 'Point'
    //////////////////////////////////////////////////////////////////////////
    template <class T>
    inline Vector3<T> GetBaryCentric(const Vector3<T>* triangle, const Vector3<T>& point, size_t stride = 0 )
    {
        stride = stride != 0 ? stride : 0;       

        auto p1 = Misc::ConvertTo <Vector3<T>>(reinterpret_cast<const char*>(triangle), stride * 0);
        auto p2 = Misc::ConvertTo <Vector3<T>>(reinterpret_cast<const char*>(triangle), stride * 1);
        auto p3 = Misc::ConvertTo <Vector3<T>>(reinterpret_cast<const char*>(triangle), stride * 2);
        
        
        Vector3<T> v0 = p2		- p1, 
                   v1 = p3		- p1, 
                   v2 = point	- p1;

        T d00 = v0.dot( v0 );
        T d01 = v0.dot( v1 );
        T d11 = v1.dot( v1 );
        T d20 = v2.dot( v0 );
        T d21 = v2.dot( v1 );

        T denom = d00 * d11 - d01 * d01;
        T v = (d11 * d20 - d01 * d21) / denom;
        T w = (d00 * d21 - d01 * d20) / denom;
        T u = T(1.0) - v - w;
        return Vector3<T>( u, v, w );

    }
    

    template <class T>
    inline T		TriangleArea(const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c )
    {
        auto ab = b - a;
        auto ac = c - a;
        auto cp = ab.cross(ac);
        return cp.length() * (T) 0.5;
    }

    template <class T>
    inline T		TriangleArea(const Vector3<T>* triangle)
    {
        return TriangleArea(triangle[0], triangle[1], triangle[2]);

    }
    

    //////////////////////////////////////////////////////////////////////////
    // \Brief: Ray/AABB intersection
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    bool IntersectsRay( const BBox3<T>& box, const Vector3<T>& orig, const Vector3<T>& dir ) 
    {
        Vector3<T> T_1, T_2; // vectors to hold the T-values for every direction
        double t_near = -FLOAT_MAX; // maximums defined in float.h
        double t_far   = FLOAT_MAX;

        for (int i = 0; i < 3; i++) { //we test slabs in every direction
            if (dir[i] == 0) { // ray parallel to planes in this direction
                if ((orig[i] < box.m_min[i]) || (orig[i] > box.m_max[i])) {
                    return false; // parallel AND outside box : no intersection possible
                }
            }
            else { // ray not parallel to planes in this direction
                T_1[i] = (box.m_min[i] - orig[i]) / dir[i];
                T_2[i] = (box.m_max[i] - orig[i]) / dir[i];

                if (T_1[i] > T_2[i]) { // we want T_1 to hold values for intersection with near plane
                    std::swap(T_1, T_2);
                }
                if (T_1[i] > t_near) {
                    t_near = T_1[i];
                }
                if (T_2[i] < t_far) {
                    t_far = T_2[i];
                }
                if ((t_near > t_far) || (t_far < 0)) {
                    return false;
                }
            }
        }        
        return true; // if we made it here, there was an intersection - YAY
    }



    //////////////////////////////////////////////////////////////////////////
    // \Brief: Ray/Triangle intersection
    //////////////////////////////////////////////////////////////////////////
    template<class T>
    inline bool	IntersectsRay(const Vector3<T>& rayStart, const Vector3<T>& rayDir,
        const Vector3<T>* triangleVerts, int numVerts,
        Vector3<T>& intSectOut, size_t stride = 0 )
    {
        
        stride = stride == 0 ? sizeof(Vector3<T>) : stride;
        auto GetPoint = [triangleVerts, stride ](int index) ->Vector3<T>
        {
            auto offset = index * stride;
            const char* data = reinterpret_cast <const char*> (triangleVerts);
            auto result = *reinterpret_cast<const Vector3<T>*>( data + offset );
            return result;
        };
        
        
        Plane<T> polyPlane(GetPoint(0), GetPoint(1), GetPoint(2));
        Vector3<T> pointOnPlane;
        if( !polyPlane.intersect( rayStart, rayDir, pointOnPlane ) )
            return false;
        //loop through edges point should be in the back of all edges
        for( int i = 0; i < numVerts; i++ )
        {
            //int first, second,  next;
            int first = i;
            int second = (i+1) % numVerts;

            Vector3<T> polyEdge	 = GetPoint(second) - GetPoint(first);
            Vector3<T> pointToEdge = pointOnPlane - GetPoint(first);
            Vector3<T> perp = pointToEdge.cross( polyEdge );

            T side = polyPlane.getNormal().dot( perp );
            if (side < (T) 0.0) 
                return false; //  outside
        }
        intSectOut = pointOnPlane;
        return true;		
    }

    

    template<class T>
    inline EulerAngle<T> QuaternionToEulerAngle( const Quat<T>& quat)
    {
        return  
            Matrix3ToEulerAngle( QuaternionToMatrix3( quat ) );
    }

    template<class T>
    inline Matrix3<T> QuaternionToMatrix3( const Quat<T>& quat )
    {
        Matrix3<T> result;
        const T* q = quat.toConstPointer();

		const T one = static_cast<T>(1.0);
		const T two = static_cast<T>(2.0);

        const T xx  = q[0] * q[0];
        const T xy  = q[0] * q[1];
        const T xz  = q[0] * q[2];
        const T xw  = q[0] * q[3];

        const T yy  = q[1] * q[1];
        const T yz  = q[1] * q[2];
        const T yw  = q[1] * q[3];

        const T zz  = q[2] * q[2];
        const T zw  = q[2] * q[3];

        result[0][0] = one - two * ( yy + zz );
        result[0][1] =		 two * ( xy - zw );
        result[0][2] =		 two * ( xz + yw );

        result[1][0] =       two * ( xy + zw );
        result[1][1] = one - two * ( xx + zz );
        result[1][2] =		 two * ( yz - xw );

        result[2][0] =       two * ( xz - yw );
        result[2][1] =		 two * ( yz + xw );
        result[2][2] = one - two * ( xx + yy );
        
		return result;
    }

    template<class T>
    inline Matrix4<T> QuaternionToMatrix4( const Quat<T>& quat )
    {
        Matrix4<T> result;
        const T* q = quat.toConstPointer();
		
		const T zero = static_cast<T>(0.0);
		const T one  = static_cast<T>(1.0);
		const T two  = static_cast<T>(2.0);

		const T xx = q[0] * q[0];
		const T xy = q[0] * q[1];
		const T xz = q[0] * q[2];
		const T xw = q[0] * q[3];

		const T yy = q[1] * q[1];
		const T yz = q[1] * q[2];
		const T yw = q[1] * q[3];

		const T zz = q[2] * q[2];
		const T zw = q[2] * q[3];

        result[0][0] = one  - two * ( yy + zz );
        result[0][1] =		  two * ( xy - zw );
        result[0][2] =		  two * ( xz + yw );
        result[0][3] = zero;

        result[1][0] =        two * ( xy + zw );
        result[1][1] = one  - two * ( xx + zz );
        result[1][2] =		  two * ( yz - xw );
        result[1][3] = zero;

        result[2][0] =        two * ( xz - yw );
        result[2][1] =		  two * ( yz + xw );
        result[2][2] = one  - two * ( xx + yy );
        result[2][3] = zero;

        result[3][0] = zero;
        result[3][1] = zero;
        result[3][2] = zero;
        result[3][3] = one;
		
		return result;
    }

    template<class T>
    inline  AxisAngle<T> QuaternionToAxisAngle( const Quat<T>& quat )
    {
		const T two  = static_cast<T>(2.0);
		const T one  = static_cast<T>(1.0);
		const T eps  = static_cast<T>(0.000125);

		const T* q = quat.toConstPointer();

        AxisAngle<T> result;	
        T angle      = acos( q[3] ) * two;
        T sin_angle  = sqrt( one - q[3] * q[3] );

        if( sin_angle < eps )
        {
            result.m_axis[0] = q[0];
            result.m_axis[1] = q[1];
            result.m_axis[2] = q[2];
        }
        else
        {
            T invScale =  one / sin_angle;

            result.m_axis[0] = q[0] * invScale;
            result.m_axis[1] = q[1] * invScale;
            result.m_axis[2] = q[2] * invScale;
        }
        result.m_angle = angle;
        return result;
    }	

    template<class T>
    inline  Matrix3<T> AxisAngleToMatrix3( const AxisAngle<T>& aAngle )
    {
        return QuaternionToMatrix3( AxisAngleToQuaternion( aAngle ) );
    }

    template<class T>
    inline Matrix4<T> AxisAngleToMatrix4( const AxisAngle<T>& aAngle )
    {
        return QuaternionToMatrix4( AxisAngleToQuaternion( aAngle ) );
    }

    template<class T>
    inline EulerAngle<T> AxisAngleToEulerAngle( const AxisAngle<T>& aAngle )
    {
        return Matrix3ToEulerAngle( QuaternionToMatrix3( AxisAngleToQuaternion( aAngle ) ) );
    }

    template<class T>
    inline Quat<T> AxisAngleToQuaternion( const AxisAngle<T>& aAngle )
    {
        T s, c;
        sinCos( aAngle.getAngle() *  T(0.5), s, c );
        Quat<T> result;
        T* resPtr = result.toPointer();
        resPtr[0] = aAngle.m_axis[0] * s;
        resPtr[1] = aAngle.m_axis[1] * s;
        resPtr[2] = aAngle.m_axis[2] * s;
        resPtr[3] = c;	
        return result;
    }

    template<class T>
    inline EulerAngle<T> Matrix3ToEulerAngle( const Matrix3<T>& mat3 )
    {
        EulerAngle<T> result;
        T epsilon = T(0.0001), ky, kz, kx;
        if (mat3[0][2] < 1 - epsilon &&  mat3[0][2] > -1 + epsilon) 
        {
            ky = -asin( mat3[0][2] );
            T c = (T)1.0/cos(ky);
            kx = atan2( mat3[1][2] * c, mat3[2][2] * c);
            kz = atan2( mat3[0][1] * c, mat3[0][0] * c);
        } 
        else 
        {
            kz = -atan2(mat3.getElement(2), T(0));
            ky = 0;
            kx = atan2(-mat3.getElement(7), mat3.getElement(4));
        }

        result[0] = kx;
        result[1] = ky;
        result[2] = kz;	

        return result;
    }

    template<class T>
    inline Matrix4<T> Matrix3ToMatrix4( const Matrix3<T>& mat3 )
    {
         Matrix4<T> result;
        
        for( int i = 0; i < 3; ++i )
        {
            result[i][0] = mat3[i][0];
            result[i][1] = mat3[i][1];
            result[i][2] = mat3[i][2];
            result[i][3] = T(0);
        }
        result[3].setXYZW( T(0), T(0), T(0), T(1));

        return result;
    }

    template<class T>
    inline  AxisAngle<T> Matrix3ToAxisAngle( const Matrix3<T>& mat3  )
    {
        return QuaternionToAxisAngle( Matrix3ToQuaternion( mat3) );
    }

    template<class T>
    inline Quat<T> Matrix3ToQuaternion( const Matrix3<T>& mat3  )
    {
        T trace = mat3.trace();
        Quat<T> quat;
        T* result  = quat.toPointer();
        if( trace > 0 )
        {
            T s = T(0.5) / sqrt( trace + T(1.0) );
            result[3] = T(0.25) / s;
            result[0] = ( mat3[2][1] - mat3[1][2] ) * s;
            result[1] = ( mat3[0][2] - mat3[2][0] ) * s;
            result[2] = ( mat3[1][0] - mat3[0][1] ) * s;			
        }
        else
        {
            //find highest axis in diagonal
            if( mat3[0][0] > mat3[1][1] && mat3[0][0] > mat3[2][2] ) //x
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat3[0][0] - mat3[1][1] - mat3[2][2]);		
                result[0] = (T)0.25 * s;
                
                s = (T)1.0 / s;
                                
                result[1] = ( mat3[0][1] + mat3[1][0] ) * s;
                result[2] = ( mat3[0][2] + mat3[2][0] ) * s;
                result[3] = ( mat3[2][1] - mat3[1][2] ) * s;
                
                
            }
            else if( mat3[1][1] > mat3[2][2] ) //y
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat3[1][1] - mat3[0][0] - mat3[2][2]);
                result[1] = (T)0.25 * s;
                
                s = (T)1.0 / s;
                            
                result[0] = ( mat3[0][1] + mat3[1][0] ) * s;				
                result[2] = ( mat3[1][2] + mat3[2][1] ) * s;	
                result[3] = ( mat3[0][2] - mat3[2][0] ) * s;
                
            }
            else //z
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat3[2][2] - mat3[0][0] - mat3[1][1]);
                result[2] = (T)0.25 * s;
                s = (T)1.0 / s;

                result[0] = ( mat3[0][2] + mat3[2][0] ) * s;
                result[1] = ( mat3[1][2] + mat3[2][1] ) * s;
                result[3] = ( mat3[1][0] - mat3[0][1] ) * s;									
            }
        }
        return result;
    }	

    template<class T>
    inline EulerAngle<T> Matrix4ToEulerAngle( const Matrix4<T>& mat4 )
    {
        EulerAngle<T> result;
        T epsilon = (T)0.0001, ky, kz, kx;
        
        if (mat4[0][2] < 1 - epsilon &&  mat4[0][2] > -1 + epsilon) 
        {
            ky = -asin( mat4[0][2] );
            float c = 1.0f/cos(ky);
            kx = atan2( mat4[1][2] * c, mat4[2][2] * c);
            kz = atan2( mat4[0][1] * c, mat4[0][0] * c);
        } 
        else 
        {
            //TODO recalculate this better!!!
            kz = 0;
            ky = -atan2(mat4.getElement(2), (T)0.0);
            kx = atan2(-mat4.getElement(9), mat4.getElement(5));
        }

        result[0] = kx;
        result[1] = ky;
        result[2] = kz;		

        return result;
    }

    template<class T>
    inline  Matrix3<T> Matrix4ToMatrix3( const Matrix4<T>& mat4 )
    {
        Matrix3<T> result;
        for( int i = 0; i < 3; ++i )
        {
            result[i][0] = mat4[i][0];
            result[i][1] = mat4[i][1];
            result[i][2] = mat4[i][2];
        }
        return result;
    }




    template<class T>
    inline AxisAngle<T> Matrix4ToAxisAngle( const Matrix4<T>& mat4  )
    {
        return QuaternionToAxisAngle( Matrix4ToQuaternion( mat4 ) );
    }

    template<class T>
    inline Quat<T> Matrix4ToQuaternion( const Matrix4<T>& mat4 )
    {
        Quat<T> result;
        T trace = mat4.trace();
        T* qPtr  = result.toPointer();
        if( trace > 0 )
        {
            T s = (T)0.5 / sqrt( trace );
            qPtr[3] = (T) 0.25 / s;
            qPtr[0] = ( mat4[2][1] - mat4[1][2] ) * s;
            qPtr[1] = ( mat4[0][2] - mat4[2][0] ) * s;
            qPtr[2] = ( mat4[1][0] - mat4[0][1] ) * s;
        }
        else
        {
            //find highest axis in diagonal
            if( mat4[0][0] > mat4[1][1] && mat4[0][0] > mat4[2][2] ) //x
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat4[0][0] - mat4[1][1] - mat4[2][2]);					
                qPtr[3] = ( mat4[2][1] - mat4[1][2] ) / s;
                qPtr[0] = (T)0.25 * s;
                qPtr[1] = ( mat4[0][1] + mat4[1][0] ) / s;
                qPtr[2] = ( mat4[0][2] + mat4[2][0] ) / s;


            }
            else if( mat4[1][1] > mat4[2][2] ) //y
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat4[1][1] - mat4[0][0] - mat4[2][2]);
                qPtr[3] = ( mat4[0][2] - mat4[2][0] ) / s;				
                qPtr[0] = ( mat4[0][1] + mat4[1][0] ) / s;
                qPtr[1] = (T)0.25 * s;
                qPtr[2] = ( mat4[1][2] + mat4[2][1] ) / s;		

            }
            else //z
            {
                T s = (T)2.0 * sqrt( (T)1.0 + mat4[2][2] - mat4[0][0] - mat4[1][1]);

                qPtr[3] = ( mat4[1][0] - mat4[0][1] ) / s;
                qPtr[0] = ( mat4[0][2] + mat4[2][0] ) / s;
                qPtr[1] = ( mat4[1][2] + mat4[2][1] ) / s;
                qPtr[2] = (T)0.25 * s;		
            }
        }

        return result;
    }



    template<class T>
    inline Matrix3<T> EulerAngleToMatrix3( const EulerAngle<T>& eAngle )
    {
         Matrix3<T> result;
        
        T  sr, sp, sy, cr, cp, cy;		
        
        sinCos( eAngle.getPitch(), sp, cp);
        sinCos( eAngle.getYaw(),   sy, cy);
        sinCos( eAngle.getRoll(),  sr, cr);

        result[0].setXYZ( cr * cy,					cr * sy,						-sr );
        result[1].setXYZ( sp * sr * cy + cp * -sy,	sp * sr * sy + cp * cy,		sp * cr );
        result[2].setXYZ( cp * sr * cy + -sp * -sy,	cp * sr * sy + -sp * cy,	cp * cr );

        return result;
    
        ////intermediate result;
        //Matrix3<T> xyRot(	   cr,		 0,			-sr, 
        //				  sp * sr,		cp,		sp * cr,
        //				  cp * sr,		-sp,		cp *cr );
        //total result
    /*	Matrix3<T> xyzRot ( cr * cy,					cr * sy,						-sr, 
                            sp * sr * cy + cp * -sy,	sp * sr * sy + cp * cy,		sp * cr,
                            cp * sr * cy + -sp * -sy,	cp * sr * sy + -sp * cy,	cp * cr 
                          );*/
                            
        //result = xRot * yRot * zRot;
                         
        //result = xRot * yRot;
        
        //bool equals = xyzRot.equals( result, (T)0.001 );
    //	assert( equals );
    }

    template<class T>
    inline Matrix4<T> EulerAngleToMatrix4( const EulerAngle<T>& eAngle )
    {
        return Matrix3ToMatrix4( EulerAngleToMatrix3( eAngle ) );
    }

    template<class T>
    inline AxisAngle<T> EulerAngleToAxisAngle( const EulerAngle<T>& eAngle )
    {
        return QuaternionToAxisAngle( EulerAngleToQuaternion( eAngle ) );
    }

    template<class T>
    inline Quat<T> EulerAngleToQuaternion( const EulerAngle<T>& eAngle )
    {
 		const T half = static_cast<T>(0.5);
		Quat<T> result;
        T sp, cp, sy, cy, sr, cr;
        T* qPtr = result.toPointer();				

        sinCos( eAngle.getPitch() *  half, sp, cp );  
        sinCos( eAngle.getYaw()   *  half, sy, cy );		
        sinCos( eAngle.getRoll()  *  half, sr, cr );
        
        T spcr = sp * cr;
        T crcp = cr * cp;
        T srsp = sr * sp;
        T cpsr = cp * sr;
            
        qPtr[0] =  cpsr * sy - spcr * cy;
        qPtr[1] = -cpsr * cy - spcr * sy;
        qPtr[2] =  srsp * cy - crcp * sy;
        qPtr[3] =  crcp * cy + srsp * sy;

        return result;
    }		



    //////////////////////////////////////////////////////////////////////////
    // \Brief: Build a transform matrix
    //////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline Math::Matrix4<T>  BuildTransform(
		const Math::Quat<T>& rot, 
        const Math::Vector3<T>& translation,
        const Math::Vector3<T>& scale)
    {
        Math::Matrix4<T> transMat, rotationMat, scaleMat;
        rotationMat = Math::QuaternionToMatrix4(rot);
        transMat.setTranslation(translation);
        scaleMat.setScale(scale);
        auto complete = scaleMat * rotationMat * transMat;
        return complete;
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief: Builds a transformation matrix, with a pivot rotation point
    //////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline Math::Matrix4f  BuildTransform(const Math::Quat<T>& rot, const Math::Vector3<T>& translation,
        const Math::Vector3<T>& scale, const Math::Vector3<T>& pivotPoint)
    {
		const T negOne = static_cast<T>(-1.0);
		
		Math::Matrix4<T> transMat, pivotMat, rotationMat, scaleMat, invPivotMat;
        rotationMat = Math::QuaternionToMatrix4(rot);

        const auto pivotTrans = pivotPoint - translation;
        pivotMat.setTranslation(pivotTrans * negOne);

        invPivotMat.setTranslation(pivotTrans);
        transMat.setTranslation(translation);
        scaleMat.setScale(scale);

        const auto complete = pivotMat * scaleMat * rotationMat * invPivotMat * transMat;
        return complete;
    }

    //////////////////////////////////////////////////////////////////////////
    //\Brief: Returns a rotation given 2 input directions
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline Quat<T> QuatFromRotation( const Math::Vector3<T>& start, const Math::Vector3<T>& end )
    {
		const T eps  = static_cast<T>	(  0.0001);
		const T half = static_cast<T>	(  0.5 );
		const T one  = static_cast<T>	(  1.0 );
		const T two  = static_cast<T>	(  2.0 );
		const T negOne = static_cast<T>	( -1.0 );
		
		auto normStart = start;
        auto normEnd   =   end;
        T d = normStart.dot( normEnd );
        Quat<T> result;
        if ( d > negOne + eps )
        {
            auto c = normStart.cross(normEnd);
            T  s	   = (T)sqrtf( (one + d ) * two );
            T invS = one / s;
            
			result.m_xyzw[0] = c[0] * invS;
            result.m_xyzw[1] = c[1] * invS;
            result.m_xyzw[2] = c[2] * invS;
            result.m_xyzw[3] = half * s;			
        }
        else
        {
            Vector3<T> axis = Vector3<T>::RIGHT.cross(normStart);
            if (axis.length() < eps )
                axis = Vector3<T>::UP.cross(normStart);
            AxisAngle<T> aa( axis, ToRadians( (T)180.0 ) );

            result = AxisAngleToQuaternion( aa );
        }
        return result;
    }


    //////////////////////////////////////////////////////////////////////////
    //\Brief: Returns a quaternion from 3 (orthogonal) axes
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline Quat<T> QuatFromAxes(const Vector3<T>& xAxis, const Vector3<T>& yAxis, const Vector3<T>& zAxis)
    {
        return Matrix3ToQuaternion( Matrix3<T>( xAxis, yAxis, zAxis ) );
    }

	/*
	@brief: Generates a look at matrix, equivalent with gluLookAt
	*/
	template<typename T>
    inline Matrix4<T> LookAtMatrix( const Vector3<T>& eye, const Vector3<T>& target, 
		const Vector3<T>& up = Vector3<T>::UP )
    {
		const T negOne = static_cast<T>(-1.0);

		auto dir     = ( target - eye ).getNormalized();
        auto side    = dir.cross( up ).getNormalized();
        auto newUp	 = side.cross( dir ).getNormalized();
        
		auto eyeOrig = eye * negOne;
		auto negDir  = dir * negOne;
		auto eyeDot  = Vector3<T>(eyeOrig.dot(side), eyeOrig.dot(newUp), eyeOrig.dot(negDir));

		Matrix4<T> result;
        result.setDirection( negDir );
        result.setRight( side );
        result.setUp( newUp );
        result.setTranslation( eyeDot );
        return result.transposed();
    }


    //////////////////////////////////////////////////////////////////////////
    // \Brief : Generate 6 look at matrices, useful for omni lights & cubemaps
    //////////////////////////////////////////////////////////////////////////
	template<typename T>
	inline std::array<Matrix4<T>, 6>       GenerateLookAtMatrices( const Math::Vector3<T>& pos )
    {
		const T zero   = static_cast<T>(0.0);
		const T one	   = static_cast<T>(1.0);
		const T negOne = static_cast<T>(-1.0);

		std::array<Matrix4<T>, 6> result;
        //TODO: optimize, only the translation part varies here!
        result[0] = LookAtMatrix( pos, Vector3<T>( one, zero, zero		),	Vector3<T>( zero, negOne, zero	) );
        result[1] = LookAtMatrix( pos, Vector3<T>( negOne, zero, zero	),  Vector3<T>( zero, negOne, zero	) );
        result[2] = LookAtMatrix( pos, Vector3<T>( zero, one, zero		),	Vector3<T>( zero,  zero,  one	) );
        result[3] = LookAtMatrix( pos, Vector3<T>( zero, negOne, zero	),  Vector3<T>( zero,  zero, negOne ) );
        result[4] = LookAtMatrix( pos, Vector3<T>( zero, zero,  one		),  Vector3<T>( zero, negOne, zero	) );
        result[5] = LookAtMatrix( pos, Vector3<T>( zero, zero, negOne	),  Vector3<T>( zero, negOne, zero	) );
        return result;
    }


    //////////////////////////////////////////////////////////////////////////
    // \Brief : Generate a quaternion from a view direction & up vector
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline Quat<T> QuatFromLookRotation(const Vector3<T>& direction, const Vector3<T>& up = Vector3<T>::UP)
    {
        const T negOne  = static_cast<T>(-1.0);
        auto dir		= direction.getNormalized();
        auto side		= dir.cross(up).getNormalized();
        auto newUp		= side.cross(dir).getNormalized();
        return QuatFromAxes(side, newUp, dir * negOne);
    }

    //////////////////////////////////////////////////////////////////////////
    // \Brief : Generate a reflection quaternion
    //////////////////////////////////////////////////////////////////////////
    //https://math.stackexchange.com/questions/2015960/finding-the-matrix-of-a-reflection-in-a-plane
    template<class T>
    inline Quat<T>		Reflect(const Quat<T>& orig, const Plane<T>& reflPlane,
        const Vector3<T>& upVec = Vector3<T>::UP, const Vector3<T>& FORWARD = Vector3<T>( (T)0, (T)0, (T)-1 ) )
    {
        const auto& norm = reflPlane.getNormal();
        const T& a = norm[0];
        const T& b = norm[1];
        const T& c = norm[2];

        const T one = static_cast<T>(1.0);
        const T two = static_cast<T>(2.0);
        const T aa = a * a;
        const T bb = b * b;
        const T cc = c * c;
        const T ab2 = -two * a * b;
        const T ac2 = -two * a * c;
        const T bc2 = -two * b * c;

        Matrix3<T> refMat =
        {
            { one - (two * aa),		ab2,				ac2				 },
            { ab2,					one - (two * bb),	bc2				 },
            { ac2,					bc2,				one - (two * cc) }
        };

        auto view = orig * FORWARD;   //calculate world direction   
        return QuatFromLookRotation(refMat * view, upVec).conjugated(); //conjugated because of column/major conversion       
    }
    


    //////////////////////////////////////////////////////////////////////////
    // \Brief: Generate a random rotation
    //////////////////////////////////////////////////////////////////////////
    template<typename T>
    inline Quat<T> RandomRotation()
    {
        AxisAngle<T> rotX( (T)1.0, (T)0.0, (T)0.0f, Math::Random( (T) 0.0, (T) TWO_PI ) );
        AxisAngle<T> rotY( (T)0.0, (T)1.0, (T)0.0f, Math::Random( (T) 0.0, (T) TWO_PI ) );
        AxisAngle<T> rotZ( (T)0.0, (T)0.0, (T)1.0f, Math::Random( (T) 0.0, (T) TWO_PI ) );

        Quat<T> qX, qY, qZ;

        qZ = AxisAngleToQuaternion( rotZ );
        qY = AxisAngleToQuaternion( rotY );
        qX = AxisAngleToQuaternion( rotX );

        return qX * qY * qZ;
    }

    //////////////////////////////////////////////////////////////////////////
    //\ Brief: Generates a new orthonormal basis, see paper of the same name
    //////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline void BranchlessONB( const Vector3<T> &n,  Vector3<T> &b1,  Vector3<T> &b2 )
    {
		const T one		= static_cast<T>(  1.0 );
		const T negOne	= static_cast<T>( -1.0 );
		
		const T x = n.getX();
        const T y = n.getY();
        const T z = n.getZ();
        
        const T sign = copysign( one, z );

        const T a = negOne / ( sign + z );
        const T b = x * y * a;
        b1 = Vector3<T>( one + sign * x * x * a, sign * b, -sign * x);
        b2 = Vector3<T>( b, sign + y * y * a, -y );
    }

   

    //////////////////////////////////////////////////////////////////////////
    // \brief: Project point 'w' on line 'v'
    //////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline  Vector3<T> ProjectPointOnLine( const Vector3<T>& w, const Vector3<T>& v )
    {
        return( v * w.dot( v ) ); //assumes v is normalized
    }
    

    template <typename T>
    inline Matrix4<T> GetRotatedMatrix(const Vector3<T>& dir1, const  Vector3<T>& planeNormal, const Vector3<T>& pivot)
    {
        Matrix4<T> rotMat;
        rotMat.setRight(dir1);
        rotMat.setUp(dir1.cross(planeNormal));
        rotMat.setDirection(planeNormal);
        rotMat.setTranslation(pivot);
        return rotMat;
    }


    //////////////////////////////////////////////////////////////////////////
    //	Brief: Returns a rotated point rotated around pivot		
    //////////////////////////////////////////////////////////////////////////
    template <typename T>
    inline Vector3<T> GetRotatedPoint( const Vector3<T>& dir1, const  Vector3<T>& planeNormal, 
        const Vector3<T>& pivot, T angle, T length )
    {
        T angleRad = ToRadians( angle );
        
        //calculate point @ the origin without any transform
        T x = cos( angleRad ) * length;
        T y = sin( angleRad ) * length;
        T z = static_cast<T>(0.0);
        
        //calculate rotation matrix
        Matrix4<T> rotMat = GetRotatedMatrix( dir1, planeNormal, pivot );
        return rotMat.multiply( Vector3<T>( x, y, z) );
    }

	template<class T>
	static inline Math::Vector3<T>		AngleToDirection(T angleDegrees)
	{
		const auto angleRad = Math::ToRadians(angleDegrees);
		return Math::Vector3<T>(std::cos(angleRad), std::sin(angleRad), static_cast<T>(0.0));
	}


   template<class T>
   inline  Quat<T> GetRotationTo( const  Vector3<T>& src, const Vector3<T>& dest,
        const Vector3<T>& fallbackAxis = Vector3<T>::ZERO )
    {
	    const T eps		= static_cast<T>(  0.000125 );
	    const T zero    = static_cast<T>(  0.0 );
		const T half	= static_cast<T>(  0.5 );
	    const T one     = static_cast<T>(  1.0 );
		const T two		= static_cast<T>(  2.0 );
		const T negOne  = static_cast<T>( -1.0 );
		const T quatEps = static_cast<T>( 1e-6 );
		const T PI		= static_cast<T>( MATH_PI );
			    
        // Based on Stan Melax's article in Game Programming Gems
        Quat<T> result;
        Vector3<T> v0 = src;
        Vector3<T> v1 = dest;
        v0.normalize();
        v1.normalize();

        T d = v0.dot( v1 );
        // If dot == 1, vectors are the same
        if( d >= one )
            return  Quat<T>::IDENTITY;
        if(d < ( quatEps - one ) )
        {
            if( fallbackAxis != Vector3<T>::ZERO )
                result = AxisAngleToQuaternion( Math::AxisAngle<T>( fallbackAxis, PI) );
            else
            {
               Vector3<T> axis = Vector3<T>::RIGHT.cross( src ); // vec3( 1, 0, 0 )
                if(axis.length() <= eps ) // pick another if collinear
                    axis = Vector3<T>::FORWARD.cross( src ); // vec3(0, 1, 0 )
                axis.normalize();
                result = AxisAngleToQuaternion( Math::AxisAngle<T>( axis, PI) );
            }
        }
        else
        {
            T s = sqrtf( ( one + d ) * two );
            T invs = one / s;
            Vector3<T> c = v0.cross( v1 ) * invs;
            result.m_xyzw = Math::Vector4<T>( c, s * half ).getNormalized();            
        }
        return result;
   }

}; //namespace
#endif