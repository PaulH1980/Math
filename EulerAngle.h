#ifndef EULER_ANGLE_H
#define EULER_ANGLE_H
#include "MathDecl.h"
#include "Vector.h"

namespace Math
{


	/*
		Euler angle based rotation, pitch is around local x axis, roll y axis and yaw is the z axis
	*/
	template<class T>
	class EulerAngle
	{
	public:
	
		EulerAngle													(													);
		EulerAngle													( T pitch, T yaw, T roll							);
		EulerAngle													( const Vector3<T>& pry								);
		EulerAngle													( const EulerAngle& other							);
		
	
		inline void							setPitch				( T val												);
		inline void							setYaw					( T val												);
		inline void							setRoll					( T val												);

		inline void							setPitchYawRoll			( const Vector3<T>& pry								);
		inline void							setPitchYawRoll			( T pitch, T yaw, T roll							);

		inline const Vector3<T>&			getPitchYawRoll			(													) const;
	
		inline T							getPitch				(													) const;
		inline T							getYaw					(													) const;
		inline T							getRoll					(													) const;


		inline T&							operator	[]			( unsigned int index								);
		inline const T&						operator	[]			( unsigned int index								) const;

		inline bool							equals					( const EulerAngle<T>& other, T epsilon = (T)0.0	) const;

		/////////////////////////////////////////////////////////////////////////
		// \brief: Returns string representation
		//////////////////////////////////////////////////////////////////////////
		inline std::string					toString				(													) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool								fromString				( const std::string& val							);


        bool                                fromJsonObject(const IO::JSonObject& obj);
        IO::JSonObject                      toJsonObject() const;



	public:

		Vector3<T>							m_pry;			//in radians
	};

 
	
	typedef EulerAngle<float>	EulerAnglef;
	typedef EulerAngle<double>	EulerAngled;


	template<class T>
	Math::EulerAngle<T>::EulerAngle()
	{
		m_pry[0] = (T) 0.0;
		m_pry[1] = (T) 0.0;
		m_pry[2] = (T) 0.0;
	}

	template<class T>
	Math::EulerAngle<T>::EulerAngle( const Vector3<T>& pry )
	{
		m_pry = pry;
	}

	template<class T>
	Math::EulerAngle<T>::EulerAngle( T pitch, T yaw, T roll )
	{
		m_pry[0] = pitch;
		m_pry[1] = roll;
		m_pry[2] = yaw;
		
	}

	template<class T>
	Math::EulerAngle<T>::EulerAngle( const EulerAngle& other )
	{
		m_pry = other.m_pry;
	}

	
	template<class T>
	void Math::EulerAngle<T>::setPitch( T val )
	{
		m_pry[0] = val;
	}


	template<class T>
	void Math::EulerAngle<T>::setYaw( T val )
	{
		m_pry[2] = val;
	}

	template<class T>
	void Math::EulerAngle<T>::setRoll( T val )
	{
		m_pry[1] = val;
	}



	template<class T>
	void Math::EulerAngle<T>::setPitchYawRoll( T pitch, T yaw, T roll )
	{
		m_pry[0] = pitch;
		m_pry[1] = roll;
		m_pry[2] = yaw;
	}

	template<class T>
	void Math::EulerAngle<T>::setPitchYawRoll( const Vector3<T>& pry )
	{
		m_pry = pry;
	}


	template<class T>
	T Math::EulerAngle<T>::getRoll() const
	{
		return m_pry[1];
	}

	template<class T>
	T Math::EulerAngle<T>::getYaw() const
	{
		return m_pry[2];
	}

	template<class T>
	T Math::EulerAngle<T>::getPitch() const
	{
		return m_pry[0];
	}

	template<class T>
	const Vector3<T>& Math::EulerAngle<T>::getPitchYawRoll() const
	{
		return m_pry;
	}

	template<class T>
	bool Math::EulerAngle<T>::equals( const EulerAngle<T>& other, T epsilon /*= (T)0.0 */ ) const
	{
		for( int i = 0; i < 3; ++i )
			if( Abs( m_pry[i] - other.m_pry[i] ) > epsilon )
				return false;
		return true;
	}


	template<class T>
	const T& Math::EulerAngle<T>::operator[]( unsigned int index ) const
	{
		return m_pry[index];
	}

	template<class T>
	T& Math::EulerAngle<T>::operator[]( unsigned int index )
	{
		return m_pry[index];
	}

	template<class T>
	std::string Math::EulerAngle<T>::toString() const
	{
		return m_pry.toString();
	}

	template<class T>
	bool Math::EulerAngle<T>::fromString(const std::string& val)
	{
		return m_pry.fromString( val );
	}


    template<class T>
    IO::JSonObject Math::EulerAngle<T>::toJsonObject() const
    {
        IO::JSonObject result;
        result["EulerAngle"] = m_pry.toJsonObject();       
        return result;
    }

    template<class T>
    bool Math::EulerAngle<T>::fromJsonObject(const IO::JSonObject& obj)
    {
        bool succeed = true;
        succeed &= m_pry.fromJsonObject(obj["EulerAngle"]);      
        return succeed;
    }

}




#endif