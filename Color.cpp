#include <Math/GenMath.h>
#include <Math/Color.h>


namespace Math
{
	const Color Color::WHITE;
	const Color Color::GRAY			(0.5f, 0.5f, 0.5f);
	const Color Color::BLACK		(0.0f, 0.0f, 0.0f);
	const Color Color::RED			(1.0f, 0.0f, 0.0f);
	const Color Color::GREEN		(0.0f, 1.0f, 0.0f);
	const Color Color::BLUE			(0.0f, 0.0f, 1.0f);
	const Color Color::CYAN			(0.0f, 1.0f, 1.0f);
	const Color Color::MAGENTA		(1.0f, 0.0f, 1.0f);
	const Color Color::YELLOW		(1.0f, 1.0f, 0.0f);
	const Color Color::TRANSPARENT	(0.0f, 0.0f, 0.0f, 0.0f);

	

	float Color::hue() const
	{
		auto val = minmax( true );
		auto chroma = val[1] - val[0];
		if( chroma < (float) FLOAT_EPSILON )
			return 0.0f;
		// Calculate and return hue
		if ( Math::Equals(m_rgba[1], val[1], (float) FLOAT_EPSILON ) )
			return ( m_rgba[2] + 2.0f * chroma - m_rgba[0] ) / ( 6.0f * chroma );
		else if ( Math::Equals( m_rgba[2], val[1], (float) FLOAT_EPSILON ) )
			return ( 4.0f * chroma - m_rgba[1] + m_rgba[0] ) / (6.0f * chroma );
		else
		{
			float r = (m_rgba[1] - m_rgba[2] ) / (6.0f * chroma );
			return (r < 0.0f) ? 1.0f + r : ((r >= 1.0f) ? r - 1.0f : r);
		}
	}

	void Color::clip(bool clipAlpha /*= false */)
	{
		int range = clipAlpha ? 4 : 3;
		for( int i = 0; i < range; ++i )
			m_rgba[i] = Math::Clamp( 0.0f, 1.0f, m_rgba[i] );
	}

	Math::Vector2f Color::minmax(bool clip /*= true */) const
	{
		float maxVal = -1000.0f;
		float minVal =  1000.0f;
		for( int i = 0; i < 3; ++i )
		{
			if( m_rgba[i] > maxVal)
				maxVal = m_rgba[i];
			if( m_rgba[i] < minVal )
				minVal = m_rgba[i];
		}

		if( clip )
		{
			minVal = Math::Clamp( 0.0f, 1.0f, minVal );
			maxVal = Math::Clamp( 0.0f, 1.0f, maxVal );
		}

		return Vector2f( minVal, maxVal );
	}

	Math::Vector4ub Color::ToRGBA8() const
	{
		Color cpy = *this;
		cpy.clip( true );

		Math::Vector4ub ret;

		for( int i = 0; i < 4; ++i )
			ret[i] = (std::uint8_t)(cpy.m_rgba[i] * 255.0f);
		return ret;
	}

	std::uint16_t Color::ToRGBB565() const
	{
		auto rgba8 = ToRGBA8();
		//auto red = ((int)255) >> 3;
		std::uint16_t r = (((int)rgba8[0] >> 3) << 11 );
		std::uint16_t g = (((int)rgba8[1] >> 2) << 5 );
		std::uint16_t b = (((int)rgba8[2] >> 3) << 0 );

		//std::uint16_t ret = ((rgba8[0] >> 3) << 11 ) | ((rgba8[1] >> 2) << 5) | (rgba8[2] >> 3);
		return  r | g | b;
	}

}