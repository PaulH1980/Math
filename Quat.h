#ifndef QUATERNION_MATH_H
#define QUATERNION_MATH_H

#include "MathDecl.h"
#include "Vector.h"
#include "Matrix.h"


namespace Math
{
	template<typename T>
	class Quat
	{
	public:
		Quat														(												);
		Quat														( T x, T y, T z, T w							);
		Quat														( const Quat& other								);
		Quat														( const T* vals									);
		
		inline void							identity				(												);
		inline void							normalize				(												);
		inline T							length					(												) const;

		Quat<T>								getNormalized			(												) const;

		inline T							getX					(												) const;
		inline T							getY					(												) const;
		inline T							getZ					(												) const;
		inline T							getW					(												) const;

		inline void							setX					( T val											);
		inline void							setY					( T val											);
		inline void							setZ					( T val											);
		inline void							setW					( T val											);

		inline void							setXYZW					( T x, T y, T z, T w							);
		
		inline void							slerp					( const Quat<T>& q1, const Quat<T>& q2, T time	);

		inline T*							toPointer				(												);
		inline const T*						toConstPointer			(												) const;
		
		inline T							dot						( const Quat<T>& other							) const;

		inline void							conjugate				(												);

		inline  Quat<T>						conjugated				(												) const;

		inline void							invert					(												);
        inline  Quat<T>						inverted                (                                               ) const;

		inline T&							operator[]				( int index										);
		inline const T&						operator[]				( int index										) const;
		inline  Quat<T>						operator *				( const Quat<T>& other							) const; 
		inline  Quat<T>						operator *				( T scale										) const; 
		inline  Quat<T>						operator +				( const Quat<T>& other							) const;
		inline  Quat<T>						operator -				( const Quat<T>& other							) const;
		inline  Quat<T>						operator /				( T scale										) const;
		inline  Vector3<T>					operator *				( const Vector3<T>& vec							) const;	


		inline bool							operator ==				( const Quat<T>& other							) const;
        inline bool                         operator !=             ( const Quat<T>& other							) const;

		inline bool							equals					( const Quat<T>& other, const T& eps, bool compSign = false) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transform a point
		//////////////////////////////////////////////////////////////////////////
		inline Vector3<T>					transformPoint			( const Vector3<T>& inputVert					) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transform a normal
		//////////////////////////////////////////////////////////////////////////
		inline Vector3<T>					transformNormal			( const Vector3<T>& normal						) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply to quaternions
		////////////////////////////////////////////////////////////////////////// 
		void								multiply				( const Quat<T>& q1, const Quat<T> q2			);


        //////////////////////////////////////////////////////////////////////////
        //\brief: Multiply this quaternion by rhs, return a copy
        //////////////////////////////////////////////////////////////////////////
        Quat<T>                             multiplied              ( const Quat<T>& rhs                            ) const;


		/////////////////////////////////////////////////////////////////////////
		// \brief: Returns string represenation
		//////////////////////////////////////////////////////////////////////////
		inline std::string					toString				(												) const;
		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool								fromString				( const std::string& val				        );			

        bool                                fromJsonObject(const IO::JSonObject& obj)
        {
            return m_xyzw.fromJsonObject(obj);               
        }


        IO::JSonObject                      toJsonObject() const
        {
            return m_xyzw.toJsonObject();
        }


	public:
		static	Quat<T>						IDENTITY;
		Vector4<T>							m_xyzw;
	};


	using Quatf = Quat<float>;
	using Quatd = Quat<double>;
	
	
	template<class T>
	Quat<T>::Quat()
	{
		m_xyzw[0] = m_xyzw[1] = m_xyzw[2] = 0.0; 
		m_xyzw[3] = 1.0;
	}

	template<class T>
	Quat<T>::Quat( const Quat& other )
	{
		m_xyzw = other.m_xyzw;
	}

	template<class T>
	Quat<T>::Quat( T x, T y, T z, T w )
	{
		m_xyzw[0] = x; 
		m_xyzw[1] = y;
		m_xyzw[2] = z; 
		m_xyzw[3] = w;
	}

	template<class T>
	Quat<T>::Quat( const T* vals )
	{
		m_xyzw[0] = vals[0]; 
		m_xyzw[1] = vals[1]; 
		m_xyzw[2] = vals[2]; 
		m_xyzw[3] = vals[3];
	}


	template<class T>
	void Quat<T>::normalize()
	{
		m_xyzw.normalize();
	}


	template<typename T>
	Quat<T> Quat<T>::getNormalized() const
	{
		return Quat<T>(m_xyzw.getNormalized().toConstPointer());
	}


	template<class T>
	T Quat<T>::length() const
	{
		return m_xyzw.length();
	}


	template<class T>
	T Quat<T>::getX() const
	{
		return m_xyzw[0];
	}

	template<class T>
	T Quat<T>::getY() const
	{
		return m_xyzw[1];
	}

	template<class T>
	T Quat<T>::getZ() const
	{
		return m_xyzw[2];
	}

	template<class T>
	T Quat<T>::getW() const
	{
		return m_xyzw[3];
	}


	template<class T>
	void Quat<T>::setX( T val )
	{
		 m_xyzw[0] = val;
	}

	template<class T>
	void Quat<T>::setY( T val )
	{
		m_xyzw[1] = val;
	}

	template<class T>
	void Quat<T>::setZ( T val )
	{
		m_xyzw[2] = val;
	}

	template<class T>
	void Quat<T>::setW( T val )
	{
		m_xyzw[3] = val;
	}


	template<typename T>
	void Quat<T>::setXYZW( T x, T y, T z, T w )
	{
		m_xyzw.setXYZW( x, y, z, w );
	}


	template<class T>
	const T* Quat<T>::toConstPointer() const
	{
		return m_xyzw.toConstPointer();
	}

	template<class T>
	T* Quat<T>::toPointer()
	{
		return m_xyzw.toPointer();
	}

	template<typename T>
	std::string Quat<T>::toString() const
	{
		return m_xyzw.toString();
	}


	template<typename T>
	bool Quat<T>::fromString(const std::string& val)
	{
		return m_xyzw.fromString( val );
	}


	





	template<class T>
	void Quat<T>::slerp( const Quat<T>& q1, const Quat<T>& q2, T time )
	{
		//clamp
		if( time <= (T)0.0 )
		{
			*this = q1;
		}
		else if( time >= (T)1.0 ){
			*this = q2;
		}
		//okay do an interpolation
		else
		{
			T dotProd = q1.dot( q2 );
			Quat<T> tmp = q2;
			//rotate over smallest angle
			if( dotProd < 0.0 )
			{
				tmp.invert();
				dotProd = -dotProd;
			}
			//no need to do spherical interpolation..
			if( dotProd > 0.99 )
			{
				m_xyzw.lerp( q1.m_xyzw, tmp.m_xyzw, time );
			}
			else
			{
				T angle = acos( ( T )dotProd );			

				*this =  ( q1 * sin( angle * ( (T)1.0 - time ) ) + tmp * sin( angle * time ) ) / sin( angle );
			}	
		}		
	}

	template<typename T>
	Quat<T> Quat<T>::operator*(const Quat<T>& other) const
	{
		Quat<T> ret;
		ret.multiply( *this, other );
		return ret;
	}


	template<typename T>
	Vector3<T> Quat<T>::operator*(const Vector3<T>& vec) const
	{
		/*Vector3<T> qVec(m_xyzw[0],m_xyzw[1],m_xyzw[2]);
		Vector3<T> cross1(qVec.cross(vec));
		Vector3<T> cross2(qVec.cross(cross1));
		auto temp = ( ( cross1 * m_xyzw[3] ) + cross2 ) * 2.0f;
		return vec + temp;*/
		return transformPoint(vec);
	}

	template<typename T>
	Quat<T> Quat<T>::inverted() const
	{
		Quat<T> ret = *this;
		ret.invert();
		return ret;
	}



	template<typename T>
	bool Quat<T>::equals(const Quat<T>& other, const T& eps, bool compSign ) const
	{
		if (compSign)
			return m_xyzw.equals(other.m_xyzw, eps);		
		//compare rotation
		for (auto i = 0; i < 4; ++i)
		{
			auto dif = std::abs(m_xyzw[i]) - std::abs(other.m_xyzw[i]);
			if (std::abs(dif) > eps)
				return false;
		}
		return true;
	}



	template<class T>
	T Quat<T>::dot( const Quat<T>& other ) const
	{
		return m_xyzw.dot( other.m_xyzw );
	}



	template<class T>
	void Quat<T>::conjugate()
	{
		m_xyzw[0] = -m_xyzw[0];
		m_xyzw[1] = -m_xyzw[1];
		m_xyzw[2] = -m_xyzw[2];			
	}

	template<class T>
	Quat<T> Quat<T>::conjugated() const
	{
		Quat<T> result(-m_xyzw[0], -m_xyzw[1], -m_xyzw[2], m_xyzw[3] );	
		return result;
	}


	template<typename T>
	void Quat<T>::invert()
	{
		for( int i =0; i < 4; ++i )
			m_xyzw[i] = -m_xyzw[i];
	}

	template<typename T>
	const T& Quat<T>::operator[]( int index ) const
	{
		return m_xyzw[index];
	}

	template<typename T>
	T& Quat<T>::operator[]( int index)
	{
		return m_xyzw[index];
	}


	template<typename T>
	 Quat<T> Quat<T>::operator+( const Quat<T>& other ) const
	{
		return Quat<T>( (m_xyzw + other.m_xyzw).toConstPointer() );
	}

	template<typename T>
	 Quat<T> Quat<T>::operator- ( const Quat<T>& other ) const
	{
		return Quat<T>( (m_xyzw - other.m_xyzw).toConstPointer() );
	}


	template<typename T>
	 Quat<T> Quat<T>::operator*( T scale ) const
	{
		return Quat<T>( m_xyzw[0] * scale, m_xyzw[1] * scale, m_xyzw[2] * scale, m_xyzw[3] * scale );
	}	

	template<typename T>
	Quat<T> Quat<T>::operator / ( T scale ) const
	{
		return Quat<T>( m_xyzw[0] / scale, m_xyzw[1] / scale, m_xyzw[2] / scale, m_xyzw[3] / scale );
	}	


	template<typename T>
	bool Quat<T>::operator==( const Quat<T>& other ) const
	{
		return this->equals(other, (T)0.0, false);
	}


    template<typename T>
    bool Quat<T>::operator!=(const Quat<T>& other) const
    {
        return !(*this == other);
    }



	template<typename T>
	void Quat<T>::identity()
	{
		m_xyzw[0] = m_xyzw[1] = m_xyzw[2] = (T)0.0;
		m_xyzw[3] = (T) 1.0;
	}

	template<typename T>
	void Quat<T>::multiply( const Quat<T>& q1, const Quat<T> q2 )
	{
		m_xyzw[0] = q1[3] * q2[0] + q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1];	
		m_xyzw[1] = q1[3] * q2[1] + q1[1] * q2[3] + q1[2] * q2[0] - q1[0] * q2[2];	
		m_xyzw[2] = q1[3] * q2[2] + q1[2] * q2[3] + q1[0] * q2[1] - q1[1] * q2[0];  
		m_xyzw[3] = q1[3] * q2[3] - q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2];	
	}

    template<typename T>
    Quat<T> Quat<T>::multiplied(const Quat<T>& rhs) const
    {
        Quat<T> tmp;
        tmp.multiply(*this, rhs);
        return tmp;
    }



	template<typename T>
	Vector3<T> Quat<T>::transformPoint( const Vector3<T>& inputVert ) const
	{
		T two = (T)2.0;
		Vector3<T> uv, uuv;
		Vector3<T> qVec( m_xyzw[0], m_xyzw[1], m_xyzw[2] );

		uv  = qVec.cross( inputVert );
		uuv = qVec.cross( uv );

		uv  *= ( two * m_xyzw[3] );
		uuv *= two;

		return inputVert + uv + uuv;
	}

	template<typename T>
	Vector3<T> Quat<T>::transformNormal( const Vector3<T>& normal ) const
	{
		Quat<T> tmp ( *this );
		tmp.conjugate();
		return tmp.transformPoint( normal );
	}



    void to_json(IO::JSonObject& obj, const Quat<float>& v);
    void from_json(const IO::JSonObject& obj, Quat<float>& v);
	
	 template<typename T>	Quat<T> Math::Quat<T>::IDENTITY = Quat<T>((T) 0.0, (T) 0.0, (T) 0.0, (T) 1.0 );
};
#endif