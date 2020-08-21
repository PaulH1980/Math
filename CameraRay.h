 #pragma once

#include "MathDecl.h"
#include "Vector.h"
#include "AABB.h"


namespace Math
{
	
	template<class T>
	class CameraRay 
	{
		public:
			CameraRay(const Math::Vector3<T>& orig, const Math::Vector3<T>& dir);
			CameraRay(const CameraRay<T>& rhs);



			CameraRay<T>& operator = (const CameraRay<T>& rhs);
			CameraRay<T>   transformed( const Math::Matrix4<T>& val ) const;
			CameraRay<T>&  transform(const Math::Matrix4<T>& val );


			
		
			Vector4<T>		m_rayDir;
			Vector4<T>		m_invRayDir;
			Vector4<T>		m_orig;
			int				m_sign[4];
	};


	class CameraRay2x2
	{
		Vector4<float>			m_oX, m_oY, m_oZ;					//origin
		Vector4<float>			m_dirX, m_dirY, m_dirZ;				//direction
		Vector4<float>			m_invDirX, m_invDirY, m_invDirZ;	//inv direction
		Vector4<float>			m_rayFlags;
	};
		
	
	template<class T>
	Math::CameraRay<T>::CameraRay(const Math::Vector3<T>& orig, const Math::Vector3<T>& dir)
		: m_rayDir( dir )
		, m_orig( orig )
	{
		for (int i = 0; i < 3; ++i)
		{
			m_invRayDir[i]	= (T)1.0 / m_rayDir[i];
			m_sign[i] = m_invRayDir[i] < (T) 0.0 ? 1 : 0;
		}		
		
	}


	template<class T>
	Math::CameraRay<T>& Math::CameraRay<T>::operator=(const Math::CameraRay<T>& rhs)
	{
		m_orig = rhs.m_orig;
		m_invRayDir = rhs.m_invRayDir;
		m_rayDir = rhs.m_rayDir;

		
		return *this;
	}

	template<class T>
	Math::CameraRay<T>::CameraRay(const CameraRay<T>& rhs)
		: m_orig(rhs.m_orig)
		, m_rayDir(rhs.m_rayDir)
		, m_invRayDir(rhs.m_invRayDir)
	{

	}


	template<class T>
	CameraRay<T> Math::CameraRay<T>::transformed(const Math::Matrix4<T>& val ) const
	{
		CameraRay<T> result = this->transform(val);
		return result;
	}

	template<class T>
	CameraRay<T>& Math::CameraRay<T>::transform(const Math::Matrix4<T>& val )
	{
		CameraRay<T> result(val * m_orig, val.transformNormal(m_rayDir));
		*this = result;
		return *this;
	}


}
