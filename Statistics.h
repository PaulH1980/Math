#ifndef STATISTICS_H
#define STATISTICS_H
#include <vector>
#include "MathDecl.h"
#include "GenMath.h"

namespace Math
{
	
	template<class T>
	class StandarDeviation
	{
	public:
		StandarDeviation( const std::vector<T>& values );

		T				  m_mean;
		T				  m_variance;
		T                 m_stDev;
		std::vector<T>	  m_values;


		T				  pdf( const T& val );
		
	private:
		void			  calcMean();
		void			  calcVariance();
		void			  calcSTDev();

	};

	

	template<class T>
	Math::StandarDeviation<T>::StandarDeviation( const std::vector<T>& values)
		: m_values( values )
	{
		if ( m_values.size())
		{
			calcMean();
			calcVariance();
			calcSTDev();
		}
		else
		{
			m_mean		= (T) 0.0;
			m_variance	= (T) 1.0;
			m_stDev		= (T) 0.0;
		}
	}




	template<class T>
	void Math::StandarDeviation<T>::calcSTDev()
	{
		m_stDev = sqrt( m_variance );
	}

	template<class T>
	void Math::StandarDeviation<T>::calcVariance()
	{
		m_variance = (T) 0.0;		
		for( int i = 0; i < m_values.size(); ++i )
		{
			T tmp = m_values[i] - m_mean;
			m_variance += tmp * tmp;
		}
        m_variance /= ( size_t ) m_values.size();
	}

	template<class T>
	void Math::StandarDeviation<T>::calcMean()
	{
		m_mean = (T) 0.0;
		for( int i = 0; i < m_values.size(); ++i )
			m_mean += m_values[i];
		m_mean /= (size_t)m_values.size();
	}


	template<class T>
	T Math::StandarDeviation<T>::pdf(const T& val)
	{
		T a = (T) 1.0 / ( sqrt( (T) 2.0 * (T) MATH_PI * m_variance ) );
		T numer = (val - m_mean) * (val - m_mean);
		T denom = (T) 2.0 * m_variance;
#if DEBUG
		assert( denom != (T) 0.0 );
#endif
		T b = -numer / denom;
		return a * std::pow( (T) MATH_EULER, b );		
	}
	
}

#endif