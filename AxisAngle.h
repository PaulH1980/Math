#ifndef AXIS_ANGLE_MATH_H
#define AXIS_ANGLE_MATH_H

#include "MathDecl.h"
#include "Vector.h"

namespace Math
{
	template<class T>
	class AxisAngle
	{
	public:
		AxisAngle													(														);
		AxisAngle													( const Vector3<T>& axis, T angle						);
		AxisAngle													( T xAxis, T yAxis, T zAxis, T angle					);
		AxisAngle													( const AxisAngle<T>& other								);

		inline void							setAxis					( const Vector3<T>& axis								);
		inline void							setAngle				( T angle												);

		inline const Vector3<T>&			getAxis					(														) const;
		inline T							getAngle				(														) const;

		std::string					        toString                (                                                       ) const;
		bool								fromString              ( const std::string& val                                );
    
        bool                                fromJsonObject          ( const IO::JSonObject& obj                             );
        IO::JSonObject                      toJsonObject            (                                                       ) const;

	
	public:
		Vector3<T>							m_axis;
		T									m_angle;
	};


	using AxisAnglef = AxisAngle<float>;
	using AxisAngled = AxisAngle<double>;

	template<class T>
	Math::AxisAngle<T>::AxisAngle()
	{
		m_axis.setXYZ( (T)0.0, (T)0.0, (T)1.0 );
		m_angle = (T)0.0;
	}

	template<class T>
	Math::AxisAngle<T>::AxisAngle( T xAxis, T yAxis, T zAxis, T angle )
	{
		m_axis.setXYZ( xAxis, yAxis, zAxis );
		m_angle = angle;
	}

	template<class T>
	Math::AxisAngle<T>::AxisAngle( const Vector3<T>& axis, T angle )
	{
		m_axis = axis;
		m_angle = angle;
	}

	template<class T>
	Math::AxisAngle<T>::AxisAngle( const AxisAngle<T>& other )
	{
		m_axis = other.m_axis;
		m_angle = other.m_angle;
	}


	template<class T>
	T Math::AxisAngle<T>::getAngle() const
	{
		return m_angle;
	}

	template<class T>
	const Vector3<T>& Math::AxisAngle<T>::getAxis() const
	{
		return m_axis;
	}

	template<class T>
	void Math::AxisAngle<T>::setAngle( T angle )
	{
		m_angle = angle;
	}

	template<class T>
	void Math::AxisAngle<T>::setAxis( const Vector3<T>& axis )
	{
		m_axis = axis;
	}


    template<class T>
    std::string AxisAngle<T>::toString() const
    {
        return Math::Vector4<T>(m_axis, m_angle).toString();
    }

    template<class T>
    bool AxisAngle<T>::fromString(const std::string& val)
    {
        Math::Vector4<T> result;
        auto succeed = result.fromString(val);
        if (!succeed)
            return false;
        m_axis = result.toVector3();
        m_angle = result.getW();
        return true;
    }

    template<class T>
    IO::JSonObject AxisAngle<T>::toJsonObject() const
    {
        IO::JSonObject result;
        result["Axis"] = m_axis.toJsonObject();
        result["Angle"] = m_angle;
        return result;
    }

    template<class T>
    bool AxisAngle<T>::fromJsonObject(const IO::JSonObject& obj)
    {
        bool succeed = true;
        succeed &= m_axis.fromJsonObject(obj["Axis"]);
        m_angle = obj.value("Angle", (T) 0.0);
        return succeed;
    }
};



#endif