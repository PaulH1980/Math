#ifndef PLANE_MATH_H
#define PLANE_MATH_H
#include "MathDecl.h"
#include "Vector.h"
#include "Matrix.h"

namespace Math
{
	enum ePlaneClassification
	{
        CLASSIFCATION_INVALID = -1,
        BEHIND_PLANE        = 0,
        ON_PLANE            = 1,
        FRONT_PLANE         = 2,
        INTERSECT_PLANE     = 3 // only for edges meshes etc.
	};
	
	template <class T>
	class Plane
	{
	public:
		Plane														(																		);
		Plane														( T x, T y, T z, T distance												);
		Plane														( const Vector3<T>& normal, T distance									);
		Plane														( const Plane<T>& other													);
		Plane														( const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& p3		);
		Plane														( const Vector3<T>& normal, const Vector3<T>& pointOnPlane				);
        Plane														( const Vector4<T>& normal_distance                     				);

        IO::JSonObject                      toJsonObject            (                                                                       ) const;
        bool                                fromJsonObject          ( const IO::JSonObject& obj                                             );
        
		/////////////////////////////////////////////////////////////////////////
		// \brief: Returns string representation
		//////////////////////////////////////////////////////////////////////////
		inline std::string			        toString				(																		) const;
		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool								fromString				( const std::string& val						);
		inline void							set						( const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& p3		);
		
		inline void							normalize				(																		);
		inline void							setDistance				( T distance															);
		inline void							setNormal				( const Vector3<T>& normal												);
		inline T							getDistance				(																		) const;
		inline const Vector3<T>&			getNormal				(																		) const;

		inline bool							intersect				( const Vector3<T>& rayStart, const Vector3<T>& rayDir, Vector3<T>& out	) const;
		inline T							distanceToPlane			( const Vector3<T>& point												) const;
		inline Vector3<T>					projectPointOnPlane		( const Vector3<T>& point												) const;

        inline ePlaneClassification         classifySphere          ( const Vector3<T>& pos, T radius ) const;
		inline ePlaneClassification 		classifyPoint			( const Vector3<T>& point, const T& epsilon =(T) 0.0					) const;
        inline void                         invert                  (                                                                       );
        inline Plane<T>                     inverted                (                                                                       ) const;
		inline bool							clip					( const Plane<T>& p1, const Plane<T>& p2, Vector3<T>& result			) const;
		inline Vector4<T>					toVector4				(																		) const;
        void                                fromVector4             ( const Vector4<T>& xyzw                                                );

        inline bool                         equals                  ( const Plane<T>& rhs, const T& dEps =(T) 0.00125, const T& nEps=(T) 0.00125 ) const;
        inline bool                         operator    ==          ( const Plane<T>& rhs                                                   ) const;
        inline bool                         operator    !=          ( const Plane<T>& rhs                                                   ) const;
		Matrix4<T>							reflectionMatrix		(																		) const;

	public:
		Vector3<T>					m_normal;
		T							m_distance;		
		static  Plane<T>			UP;			
	};

    template <class T>
    Math::Plane<T>::Plane(const Vector4<T>& normal_distance)
    {
        fromVector4(normal_distance);
    }
    

	using Planef = Plane<float>;
	using Planed = Plane<double>;
   

	template <class T>
	Plane<T>::Plane()
	{
		m_normal.setXYZ( (T)0.0, (T)0.0, (T)1.0 );
		m_distance = 0.0;
	}

	template <class T>
	Plane<T>::Plane( T x, T y, T z, T distance )
	{
		m_normal.setXYZ(x, y, z);
		m_distance =  distance;
	}

	template <class T>
	Plane<T>::Plane( const Vector3<T>& normal, T distance )
	{
		m_normal	= normal;
		m_distance	= distance;
	}

	template <class T>
	Plane<T>::Plane( const Plane<T>& other )
	{
		m_normal	= other.m_normal;
		m_distance	= other.m_distance;
	}
	template <class T>
	Plane<T>::Plane ( const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& p3 )
	{
		set( p1, p2, p3 );
	}

	template <class T>
	Plane<T>::Plane( const Vector3<T>& normal, const Vector3<T>& pointOnPlane )
	{
		m_normal	= normal;
		m_distance	= -m_normal.dot( pointOnPlane );
	}

	template <class T>
	void Plane<T>::normalize()
	{
		T invLength = ( (T)1.0 ) / m_normal.length();

		m_normal	*= invLength;
		m_distance	*= invLength;
	}

	template <class T>
	void Plane<T>::invert()
	{
		m_normal = m_normal.inverted();
		m_distance = -m_distance;
	}


    template <class T>
    void Math::Plane<T>::fromVector4(const Vector4<T>& xyzw)
    {
        m_normal = xyzw.toVector3();
        m_distance = xyzw[3];
    }

	template <class T>
	bool Plane<T>::operator!=(const Plane<T>& rhs) const
	{
		return !(*this == rhs);
	}

	template <class T>
	bool Plane<T>::operator==(const Plane<T>& rhs) const
	{
		return this->equals(rhs, (T) 0.0, (T) 0.0);
	}

	template <class T>
	bool Plane<T>::equals(const Plane<T>& rhs, const T& dEps, const T& nEps) const
	{
		return m_normal.equals(rhs.m_normal, nEps) && (fabs(m_distance - rhs.m_distance) <= dEps);
	}


	template <class T>
	void Plane<T>::setDistance( T distance )
	{
		m_distance = distance;
	}

	template <class T>
	void Plane<T>::setNormal( const Vector3<T>& normal )
	{
		m_normal = normal;
	}

	template <class T>
	T Plane<T>::getDistance() const
	{
		return m_distance;
	}
	
	template <class T>
	const Vector3<T>&  Plane<T>::getNormal() const
	{
		return m_normal;
	}

	template <class T>
	T Plane<T>::distanceToPlane( const Vector3<T>& point ) const
	{
		return (  m_normal.dot( point ) + m_distance );
	}


	template <class T>
	Vector3<T> Plane<T>::projectPointOnPlane( const Vector3<T>& point ) const
	{
		auto distance = distanceToPlane( point ) * (T) -1.0;
		return point + (  m_normal * distance );
	}


	template <class T>
	void Plane<T>::set(const Vector3<T>& p1, const Vector3<T>& p2, const Vector3<T>& p3)
	{
		m_normal		= ( p3 - p1 ) ^ ( p2 - p1 );
		m_normal.normalize();
		m_distance		= -m_normal.dot( p1 );	
	}



	template <class T>
	std::string Plane<T>::toString() const
	{
		Math::Vector4<T> xyzw( m_normal[0], m_normal[1], m_normal[2], m_distance );
		return xyzw.toString();
	}


	template <class T>
	bool Plane<T>::fromString(const std::string& val)
	{
		Vector4<T> res;
		if( res.fromString( val ) )
		{
			m_normal[0] = (T)res[0];
			m_normal[1] = (T)res[1];
			m_normal[2] = (T)res[2];
			m_distance  = (T)res[3];
			return true;
		}
		return false;
	}

    template <class T>
    Plane<T> Plane<T>::inverted() const
    {
        Plane<T> result = *this;
        result.invert();
        return result;
    }



	template <class T>
	Vector4<T> Plane<T>::toVector4() const
	{
		return Vector4<T>( m_normal[0], m_normal[1], m_normal[2], m_distance );
	}

	template <class T>
	bool Plane<T>::intersect(const Vector3<T>& rayStart, const Vector3<T>& rayDir, Vector3<T>& out) const
	{
		//Math::Plane<T> polyPlane( triangleVerts[0], triangleVerts[1], triangleVerts[2] );
		T NdotDir = rayDir.dot( m_normal );
		if( NdotDir == (T) 0.0 ) // parallel
			return false;
		
		T t =-( m_normal.dot( rayStart ) + m_distance )  / NdotDir;
		if ( t < (T) 0.0 )
			return false; // the triangle is behind 

		out =  rayStart + ( rayDir * t );
		return true;
	}

    template <class T>
    bool Plane<T>::fromJsonObject(const IO::JSonObject& obj)
    {
        bool succeed = true;
        succeed &= m_normal.fromJsonObject(obj["Normal"]);
        m_distance = obj.value("Distance", (T) 0.0);
        return succeed;
    }

    template <class T>
    IO::JSonObject Plane<T>::toJsonObject() const
    {
        IO::JSonObject result;
        result["Normal"] = m_normal.toJsonObject();
        result["Distance"] = m_distance;
        return result;
    }



    template <class T>
    ePlaneClassification Plane<T>::classifySphere( const Vector3<T>& pos, T radius ) const
    {
         T dist = distanceToPlane( pos );

        if(dist > radius)
            return FRONT_PLANE; // The sphere is in front of the plane
        else if(dist < -radius)
            return BEHIND_PLANE; // The sphere is behind the plane

        return INTERSECT_PLANE; // The sphere collides/straddles with the plane
    }


	
	template <class T>
	ePlaneClassification Plane<T>::classifyPoint( const Vector3<T>& point, const T& epsilon ) const
	{
		T dist = distanceToPlane( point );
		if( dist < (T) -epsilon )
			return BEHIND_PLANE;
		else if (dist > (T)epsilon )
			return FRONT_PLANE;
		return ON_PLANE;
	}
	
	template <class T>
	bool  Plane<T>::clip( const Plane<T>& p1, const Plane<T>& p2, Vector3<T>& result ) const
	{
		Vector3<T> c1 = p1.getNormal() ^ p2.getNormal();
		T den = this->m_normal.dot(c1);
		if( Abs( den ) < (T)0.000125 )
			return false;
		c1 *= -m_distance;
		
		Vector3<T> c2 = ( p2.getNormal() ^ m_normal ) * p1.getDistance();
		Vector3<T> c3 = ( m_normal ^ p1.getNormal() ) * p2.getDistance();;

		c1 -= c2;
		c1 -= c3;
		c1 /= den;
		result = c1;	
		return true;
	}


	template <class T>
	Matrix4<T> Plane<T>::reflectionMatrix() const
	{
		T a = m_normal[0];
		T b = m_normal[1];
		T c = m_normal[2];
		T d = m_distance;
		//http://ami.ektf.hu/uploads/papers/finalpdf/AMI_40_from175to186.pdf
		return Matrix4<T>
		(
			(T)-2.0 * a * a + (T)1.0,	(T)-2.0 * b * a,				(T)-2.0 * c * a,			(T)0.0,
			(T)-2.0 * a * b,			(T)-2.0 * b * b + (T)1.0,		(T)-2.0 * c * b,			(T)0.0,
			(T)-2.0 * a * c,			(T)-2.0 * b * c,				(T)-2.0 * c * c + (T)1.0,	(T)0.0,
			(T)-2.0 * a * d,			(T)-2.0 * b * d,				(T)-2.0 * c * d,			(T)1.0 
		);
	}


	
	template <class T>	Plane<T> Plane<T>::UP = Math::Plane<T>( (T)0.0, (T)0.0, (T)1.0, (T)0.0 );


    void to_json(IO::JSonObject& obj, const Plane<float>& v);
    void from_json(const IO::JSonObject& obj, Plane<float>& v);

}; //namespace








#endif