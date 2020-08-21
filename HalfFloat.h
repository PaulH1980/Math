#ifndef HALF_FLOAT_H
#define HALF_FLOAT_H
#include "MathDecl.h"
#include "Vector.h"

namespace Math
{
	struct HalfFloatLUT
	{
		//////////////////////////////////////////////////////////////////////////
		// \brief: 
		//////////////////////////////////////////////////////////////////////////
		HalfFloatLUT												(												);
		
		std::uint32_t					m_mantissaTable[2048];
		std::uint32_t					m_offsetTable[64];
		std::uint32_t					m_exponentTable[64];
		std::uint32_t					m_baseTable[1024];
		std::uint32_t					m_shiftTable[512];
	};
	
	//////////////////////////////////////////////////////////////////////////
	// \brief: 16 bits precision floating point, 
	//		   see details here: ftp://www.fox-toolkit.org/pub/fasthalffloatconversion.pdf
	//////////////////////////////////////////////////////////////////////////
	class HalfFloat
	{
	public:
		HalfFloat													(												);
		HalfFloat													( float val										);												
		HalfFloat													( const HalfFloat& other						);
		
		inline float					toFloat						(												) const;

		inline void						fromFloat					( float val										);

		inline HalfFloat				operator +					(  const HalfFloat& other						) const;
		inline HalfFloat				operator -					(  const HalfFloat& other						) const;
		inline HalfFloat				operator *					(  const HalfFloat& other						) const;
		inline HalfFloat				operator /					(  const HalfFloat& other						) const;

		inline HalfFloat&				operator +=					(  const HalfFloat& other						);
		inline HalfFloat&				operator -=					(  const HalfFloat& other						);
		inline HalfFloat&				operator *=					(  const HalfFloat& other						);
		inline HalfFloat&				operator /=					(  const HalfFloat& other						);

		inline HalfFloat&				operator =					(  const HalfFloat& other						);

		inline bool						operator ==					(  const HalfFloat& other						) const;
		inline bool						operator !=					(  const HalfFloat& other						) const;

		inline std::uint16_t*			toPointer					(												);
		inline const std::uint16_t*		toConstPointer				(												) const;

		inline std::uint16_t			getData						(												) const;
		inline void						setData						( std::uint16_t data									);

	protected:

	private:
		std::uint16_t					m_floatData;
		static	HalfFloatLUT			g_LUTs;		
	};	


	

	float HalfFloat::toFloat() const
	{
		int floatBits = g_LUTs.m_mantissaTable[ g_LUTs.m_offsetTable[ m_floatData >> 10 ] 
						+ (m_floatData & 0x3FF ) ] 
						+ g_LUTs.m_exponentTable[ m_floatData >> 10 ];  
		return *( ( float* )&floatBits );
	}

	void HalfFloat::fromFloat( float val )
	{
		std::uint32_t floatBits = *((std::uint32_t*)&val);
		m_floatData = (std::uint16_t)g_LUTs.m_baseTable[( floatBits >> 23 ) & 0x1FF];
		m_floatData += (std::uint16_t)( ( floatBits & 0x007FFFFF ) >> g_LUTs.m_shiftTable[ ( floatBits >> 23 ) & 0x1FF ] );
	}

	std::uint16_t* HalfFloat::toPointer()
	{
		return &m_floatData;
	}

	const std::uint16_t* HalfFloat::toConstPointer() const
	{
		return &m_floatData;
	}

	std::uint16_t HalfFloat::getData() const
	{
		return m_floatData;
	}

	void HalfFloat::setData( std::uint16_t data )
	{
		m_floatData = data;
	}
	
	HalfFloat HalfFloat::operator+( const HalfFloat& other ) const
	{
		HalfFloat result( this->toFloat() + other.toFloat() );
		return result;	
	}

	HalfFloat HalfFloat::operator-( const HalfFloat& other ) const
	{
		HalfFloat result( this->toFloat() - other.toFloat() );
		return result;	
	}

	HalfFloat HalfFloat::operator*( const HalfFloat& other ) const
	{
		HalfFloat result( this->toFloat() * other.toFloat() );
		return result;	
	}

	HalfFloat HalfFloat::operator/( const HalfFloat& other ) const
	{
		HalfFloat result( this->toFloat() / other.toFloat() );
		return result;	
	}

	HalfFloat& HalfFloat::operator+=( const HalfFloat& other )
	{
		fromFloat( this->toFloat() + other.toFloat() );
		return *this;
	}

	HalfFloat& HalfFloat::operator-=( const HalfFloat& other )
	{
		fromFloat( this->toFloat() - other.toFloat() );
		return *this;
	}

	HalfFloat& HalfFloat::operator*=( const HalfFloat& other )
	{
		fromFloat( this->toFloat() * other.toFloat() );
		return *this;
	}

	HalfFloat& HalfFloat::operator/=( const HalfFloat& other )
	{
		fromFloat( this->toFloat() / other.toFloat() );
		return *this;
	}

	HalfFloat& HalfFloat::operator=( const HalfFloat& other )
	{
		m_floatData = other.m_floatData;
	}

	bool HalfFloat::operator==( const HalfFloat& other ) const
	{
		return (m_floatData == other.m_floatData);
	}

	bool HalfFloat::operator!=( const HalfFloat& other ) const
	{
		return (m_floatData != other.m_floatData);
	}
}; //namespace


#endif