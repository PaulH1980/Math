#ifndef COLOR_H
#define COLOR_H
#include "MathDecl.h"
#include "Vector.h"

// Defined by Windows headers
#undef TRANSPARENT

namespace Math
{
	static const std::uint32_t RGBA_OFFSETS[] =
	{
		0,  8, 16, 24
	};
	

	static inline std::uint16_t	RGB888ToRGB565(std::uint8_t r, std::uint8_t g, std::uint8_t b)
	{
		std::uint16_t red = r >> 3;
		std::uint16_t green = g >> 2;
		std::uint16_t blue = b >> 3;

		return ( ( red << 11 ) | ( green << 5 ) | ( blue ) );
	}

	
	class Color
	{
	public:

		Color(){
			m_rgba[0] = m_rgba[1] = m_rgba[2] = m_rgba[3] = 1.0f;
		}

		Color( float rgba )
		{
			m_rgba[0] = m_rgba[1] = m_rgba[2] = m_rgba[3] = rgba;
		}


		Color( float r, float g, float b, float a )
		{
			m_rgba[0] = r; m_rgba[1] = g; m_rgba[2] = b; m_rgba[3] = a;
		}
		Color( float r, float g, float b )
		{
			m_rgba[0] = r; m_rgba[1] = g; m_rgba[2] = b; m_rgba[3] = 1.0f;
		}
		Color( const Vector4<float>& rhs )
			: m_rgba( rhs ){}
		
		Color( const Color& rhs )
			: m_rgba( rhs.m_rgba ){}
		
		Color( const Color& rhs, float a )
		{ 
			m_rgba[0] = rhs.m_rgba[0]; m_rgba[1] = rhs.m_rgba[1];  m_rgba[2] = rhs.m_rgba[2]; m_rgba[3] = a;
		}
		
		Color( unsigned rgba )
		{
			for( int i = 0; i < 4; ++i )
				m_rgba[i] = (float)((rgba >> RGBA_OFFSETS[i]) & 0xFF )/255.0f;
		}

		Color( std::uint16_t rgb565 )
		{
			m_rgba[0] =  (float)(((rgb565 >> 11 ) &  ( 31  )) << 3) / 255.0f;
			m_rgba[1] =  (float)(((rgb565 >>  5 ) &  ( 63  )) << 2) / 255.0f;
			m_rgba[2] =  (float)(((rgb565 >>  0 ) &  ( 31  )) << 3) / 255.0f;
			m_rgba[3] = 0.0f;	
		}

		Color( Vector4<std::uint8_t> val )
		{
			for( int i = 0; i < 4; ++i )
				m_rgba[i] = val[i] / 255.0f;
		}


		
		bool operator == ( const Color& rhs )
		{
			return m_rgba == rhs.m_rgba;
		}

		bool operator != ( const Color& rhs )
		{
			return !(*this == rhs );
		}

		Color  operator * ( const Color& rhs ) const{ return Color( m_rgba * rhs.m_rgba );}
		Color  operator - ( const Color& rhs ) const{ return Color( m_rgba - rhs.m_rgba );}
		Color  operator + ( const Color& rhs ) const{ return Color( m_rgba + rhs.m_rgba );}
		Color  operator * ( float val		 ) const{ return Color( m_rgba * val );}

		bool  Equals( const Color& rhs){ return m_rgba.equals( rhs.m_rgba, 0.0f ); }

		float sumRGB()  const { return m_rgba[0] + m_rgba[1] + m_rgba[2]; }
		float average() const { return sumRGB() * 0.333333f; }
		float luma()    const { return (m_rgba[0] * 0.299f +  m_rgba[1] * 0.587f + m_rgba[2] * 0.114f); }
		float chroma()  const
		{
			auto range = minmax();
			return range[1] - range[0];
		}

		float hue() const;

		std::string toString() const{ return m_rgba.toString(); }
		bool fromString( const std::string& rgba ) { return m_rgba.fromString( rgba ); }
		
		std::uint32_t ToUint()  const {
			std::uint32_t val = 0;
			for( int i = 0; i < 4; ++i )
				val |= ( (unsigned)( m_rgba[i] * 255.0f ) ) << RGBA_OFFSETS[i];
			return val;
		}

		Color   absolute() const {
			auto rgba = m_rgba;
			rgba.absolute();
			return Color( rgba );
		}

		float Red() const{
			return m_rgba[0];
		}

		float Green() const{
			return m_rgba[1];
		}

		float Blue() const{
			return m_rgba[2];
		}

		float Alpha() const{
			return m_rgba[3];
		}


		std::uint16_t ToRGBB565() const;

		Vector4<std::uint8_t>  ToRGBA8() const;




		void clip( bool clipAlpha = false );
		void invert( bool invertAlpha = false )
		{
			int range = invertAlpha ? 4 : 3;
			for( int i = 0; i < range; ++i )
				m_rgba[i] = 1.0f - m_rgba[i];
		}

		void abs(){  m_rgba.absolute(); }

		float range() const 
		{
			auto range = minmax( false );
			return range[1] - range[0];
		}

		float lightness() const
		{
			auto range = minmax( false );
			return ( range[1] + range[0] ) * 0.5f;
		}

		Vector2<float> minmax  ( bool clip = true ) const;
		
		Vector4<float>		m_rgba;


		/// Opaque white color.
		static const Color WHITE;
		static const Color GRAY;
		static const Color BLACK;
		static const Color RED;
		static const Color GREEN;
		static const Color BLUE;
		static const Color CYAN;
		static const Color MAGENTA;
		static const Color YELLOW;
		static const Color TRANSPARENT;
	};



}

#endif