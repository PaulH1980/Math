#ifndef AXIS_ALIGNED_BBOX_MATH_H
#define AXIS_ALIGNED_BBOX_MATH_H
#include <array>
#include <limits>
#include <string>

#include <IO/SerializableBase.h>
#include <IO/JSonSerializer.h>
#include "MathDecl.h"
#include "Matrix.h"
#include "Vector.h"

namespace Math
{
	static const int AABB_POINTS = 8;

	//////////////////////////////////////////////////////////////////////////
	// \brief: 2D Axis Aligned Bounding Box
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class BBox2 
	{
	public:

		BBox2														(												);
		BBox2														( const Vector2<T>& min, const Vector2<T>& max	);
		BBox2														( T minX, T minY, T maxX, T maxY				);
		BBox2														( const BBox2<T>& other							);
        ~BBox2                                                      (                                               ) = default;

		inline Vector2<T>					getSize					(												) const;
		inline Vector2<T>					getHalfSize				(												) const;
		inline void							clearBounds				(												);
		inline void							updateBounds			( const Vector2<T>& point						);
		inline void							updateBounds			( const BBox2<T>& other							);

		inline BBox2<T>&					operator +=				( const BBox2<T>& other							);
		inline BBox2<T>&					operator +=				( const Vector2<T>& translate					);

        inline bool                         operator ==             ( const BBox2<T>& other                         ) const;
        inline bool                         operator !=             ( const BBox2<T>& other                         ) const;

        
		inline bool							intersects				( const BBox2<T>& other							) const;
		inline bool							insideBounds			( const Vector2<T>& point						) const;
		
		inline bool							isValid					(												) const;

		inline const Vector2<T>&			getMin					(												) const;
		inline const Vector2<T>&			getMax					(												) const;

		inline void							expand					( const T& amount								);
		inline BBox2<T>						expanded				( const T& amount								) const;
		inline void							setMin					( const Vector2<T>& min							);
		inline void							setMax					( const Vector2<T>& max							);
        inline bool                         equals                  ( const BBox2<T>& other, const T eps = (T) 0.00125 ) const;


		inline Vector2<T>					getCenter				(												) const;
		
		inline std::string					toString				(												) const;
		bool								fromString				( const std::string& val						);
        
        IO::JSonObject                      toJsonObject() const;
        bool                                fromJsonObject( const IO::JSonObject& obj );

	public:

		Vector2<T>							m_min,
											m_max;
	};

   
    

	//////////////////////////////////////////////////////////////////////////
	// \brief: 3D Axis Aligned Bounding Box
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class BBox3 
	{
        
	public:
		BBox3														(												);
		BBox3														( const Vector3<T>& min, const Vector3<T>& max	);
		BBox3														( T minX, T minY, T minZ, T maxX, T maxY, T maxZ);
		BBox3														( const BBox3<T>& other							);
		BBox3														( const Vector3<T>* points, unsigned numPoints	);
		BBox3														( const BBox2<T>& bounds						);
		BBox3														( const Vector3<T>& center, T radius            );
        ~BBox3                                                      (                                               ) = default;

		inline Vector3<T>					getSize					(												) const;
		inline Vector3<T>					getHalfSize				(												) const;

		inline void							clearBounds				(												);
		inline void							updateBounds			( const Vector3<T>& point						);
		inline void							updateBounds			( const BBox3<T>& other							);

		inline bool							intersects				( const BBox3<T>& other							) const;
		inline bool							insideBounds			( const Vector3<T>& point						) const;

        inline void                         scale                   ( const T& factor                               );
        inline BBox3<T>	                    scaled                  ( const T& factor                               ) const;


        inline bool                         equals                  ( const Math::BBox3<T>& other, const T eps =(T) 0.00125 ) const;

		inline bool							isValid					(												) const;
		inline void							expand					( const T& amount								);
		inline BBox3<T>						expanded				( const T& amount								) const;	

		inline Vector3<T>					getMin					(												) const;
		inline Vector3<T>					getMax					(												) const;

		inline BBox3<T>&					operator +=				( const BBox3<T>& other							);	
		inline BBox3<T>&					operator +=				( const Vector3<T>& other                       );

        inline bool                         operator ==             ( const BBox3<T>& other                         ) const;
        inline bool                         operator !=             ( const BBox3<T>& other                         ) const;



		inline BBox2<T>						toBBox2					(												) const;


        inline bool                         isBigger                 ( const T value                                ) const;





		inline void							setMin					( const Vector3<T>& min							);
		inline void							setMax					( const Vector3<T>& max							);
		inline void							set						( const Vector3<T>& min, const Vector3<T>& max	);
		inline void							set						( const Vector3<T>* points, unsigned numPoints	);

		inline		Vector3<T>				getCenter				(												)const;

		inline std::string					toString				(												) const;
		bool								fromString				( const std::string& val						);

		inline BBox3<T>						transform				( const Matrix4<T>& trans						) const;

        
		
        inline std::array<Vector3<T>, AABB_POINTS>    
											toPoints                (                                               ) const;


        IO::JSonObject                      toJsonObject()  const;

        bool                                fromJsonObject(const IO::JSonObject& obj);

	    public:

		Vector4<T>							m_min,
											m_max;
	};

   
	
	
	
	using BBox2i = BBox2<int>;
	using BBox2f = BBox2<float>;
	using BBox3i = BBox3<int>;
	using BBox3f = BBox3<float>;	


	template <class T>
	BBox2<T>::BBox2()
	{
		clearBounds();
	}

	template <class T>
	BBox2<T>::BBox2( const BBox2<T>& other )
	{
		m_min = other.m_min;
		m_max = other.m_max;
	}

	template <class T>
	BBox2<T>::BBox2( T minX, T minY, T maxX, T maxY )
	{
		m_min.setXY( minX, minY );
		m_max.setXY( maxX, maxY );
	}

	template <class T>
	BBox2<T>::BBox2( const Vector2<T>& min, const Vector2<T>& max )
	{
		m_min = min;
		m_max = max;
	}


	template <typename T>
	Vector2<T> BBox2<T>::getHalfSize() const
	{
		return getSize() * T(0.5f);
	}


	template <class T>
	void BBox2<T>::clearBounds()
	{
		m_min[0] = m_min[1] = std::numeric_limits<T>::max();
		m_max[0] = m_max[1] = std::numeric_limits<T>::lowest();
	}

	template <typename T>
	void BBox2<T>::expand(const T& amount)
	{
		for( int i = 0; i < 2; ++i )
		{
			m_min[i] -= amount;
			m_max[i] += amount;
		}
	}

	template <typename T>
	BBox2<T> BBox2<T>::expanded(const T& amount) const
	{
		BBox2<T> result = *this;
		result.expand(amount);
		return result;
	}

	

	template <typename T>
	Vector2<T> BBox2<T>::getSize() const
	{
		return m_max - m_min;
	}



	template <class T>
	void BBox2<T>::updateBounds( const BBox2<T>& other )
	{
		updateBounds( other.m_min );
		updateBounds( other.m_max );
	}


	template <class T>
	void BBox2<T>::updateBounds( const Vector2<T>& point )
	{
		for( int i = 0; i < 2; ++i )
		{
			if( point[i] < m_min[i] )
				m_min[i] = point[i];
			if( point[i] > m_max[i] )
				m_max[i] = point[i];
		}
	}

	template <typename T>
	void BBox3<T>::set(const Math::Vector3<T>& min, const Math::Vector3<T>& max)
	{
		setMin( min);
		setMax( max );
	}

	template <class T>
	void BBox2<T>::setMax( const Vector2<T>& max )
	{
		m_max = max;
	}

	template <class T>
	void BBox2<T>::setMin( const Vector2<T>& min )
	{
		m_min = min;
	}

	template <class T>
	const Vector2<T>& BBox2<T>::getMax() const
	{
		return m_max;
	}

	template <class T>
	const Vector2<T>& BBox2<T>::getMin() const
	{
		return m_min;
	}

    template <typename T>
    bool BBox2<T>::operator!=(const BBox2<T>& other) const
    {
        return !(*this == other);
    }

    template <typename T>
    bool BBox2<T>::operator==(const BBox2<T>& other) const
    {
        return m_min == other.m_min &&
            m_max == other.m_max;
    }



	template <class T>
	bool BBox2<T>::isValid() const
	{
		for( int i = 0; i < 2; ++i )
			if( m_min[i] >= m_max[i] )
				return false;
		return true;
	}

	template <class T>
	bool BBox2<T>::insideBounds( const Vector2<T>& point ) const
	{
		for( int i = 0; i < 2; ++i )
		{
			if( point[i] < m_min[i] || point[i] > m_max[i] )
				return false;
		}
		return true;
	}


    template <typename T>
    bool BBox2<T>::equals(const BBox2<T>& other, const T eps /*= (T) 0.00125 */) const
    {
        return getMin().equals(other.getMin(), eps) && 
               getMax().equals(other.getMax(), eps);
    }
	


	template <typename T>
	BBox2<T>& BBox2<T>::operator+=(const BBox2<T>& other)
	{
		this->updateBounds( other );
		return *this;
	}

	template <typename T>
	BBox2<T>& BBox2<T>::operator+=(const Vector2<T>& translate )
	{
		m_min += translate;
		m_min += translate;
		return *this;
	}



	template <typename T>
	Vector2<T> BBox2<T>::getCenter() const
	{
		return( m_min + ( m_max - m_min ) * (T)0.5 );
	}


    template <typename T>
    bool BBox2<T>::fromJsonObject(const IO::JSonObject& obj)
    {
        bool succeed = true;

        succeed &= m_min.fromJsonObject(obj["Min"]);
        succeed &= m_max.fromJsonObject(obj["Max"]);

        return succeed;
    }

    template <typename T>
    IO::JSonObject BBox2<T>::toJsonObject() const
    {
        IO::JSonObject result;
        result["Min"] = m_min.toJsonObject();
        result["Max"] = m_max.toJsonObject();
        return result;
    }





/*
	template <typename T>
	eBoundsIntersection BBox2<T>::classifyBounds(const BBox2<T>& other) const
	{
		if( other.m_min[0] >= m_min[0] && other.m_max[0] <= m_max[0] &&
			other.m_max[1] >= m_min[1] && other.m_max[1] <= m_max[1]  )
			return INSIDE;

        eBoundsIntersection code = intersects( other ) ? INTERSECT : OUTSIDE;
		return code;
	}*/

	template <class T>
	bool BBox2<T>::intersects( const BBox2<T>& other ) const
	{
		for( int i = 0; i < 2; ++i )
		{
			if( other.m_min[i] >= m_max[i] || other.m_max[i] < m_min[i] )
				return false;

		}
		return true;
	}


	template <typename T>
	std::string BBox2<T>::toString() const
	{
        std::string buf;
        buf += m_min.toString() + " " + m_max.toString();
        return buf;
	}



	template <typename T>
	bool BBox2<T>::fromString(const std::string& val)
	{
		if( !m_min.fromString( val ) )
			return false;
		return m_max.fromString( val );
	}



	template <class T>
	BBox3<T>::BBox3()
	{
		clearBounds();
	}


	template <class T>
	BBox3<T>::BBox3( const BBox3<T>& other )
		: m_min( other.m_min )
		, m_max( other.m_max )
	{
		
	}


	template <typename T>
	BBox3<T>::BBox3(const Vector3<T>& center, T radius)
	{

		setMin(center - Math::Vector3<T>(radius));
		setMax(center + Math::Vector3<T>(radius));
	}


	template <typename T>
	BBox3<T>::BBox3(const BBox2<T>& bounds)
		: m_min( bounds.m_min.toVector3() )
		, m_max( bounds.m_max.toVector3() )
	{

	}

	template <typename T>
	BBox3<T>::BBox3(const Vector3<T>* points, unsigned numPoints)
	{
		set( points, numPoints );
	}
	
	
	template <typename T>
	BBox2<T> BBox3<T>::toBBox2() const
	{
		return BBox2<T>(getMin().toVector2(), getMax().toVector2());
	}

	template <typename T>
	BBox3<T> BBox3<T>::transform(const Matrix4<T>& trans) const
	{
		auto points = toPoints();
        BBox3<T> bounds;
        for(int i = 0; i < 8; ++i)
        {
            points[i] = trans.multiply( points[i] );         
            bounds.updateBounds( points[i] );
        }
		return bounds;
	
		
	}

	template <class T>
	BBox3<T>::BBox3( T minX, T minY, T minZ, T maxX, T maxY, T maxZ )
	{
		m_min.setXYZW(minX, minY, minZ, static_cast<T>(0.0));
		m_max.setXYZW(maxX, maxY, maxZ, static_cast<T>(0.0));
	}

	template <typename T>
	void BBox3<T>::set(const Vector3<T>* points, unsigned numPoints)
	{
		clearBounds();
		while( numPoints-- )
			updateBounds( points[numPoints] );
	}


    template <typename T>
    bool BBox3<T>::equals(const Math::BBox3<T>& other, const T eps ) const
    {
        return getMin().equals(other.getMin(), eps) && 
               getMax().equals(other.getMax(), eps) ;
    }

	

	template <class T>
	BBox3<T>::BBox3( const Vector3<T>& min, const Vector3<T>& max )
	{
		setMin(min);
		setMax(max);
	}

	template <typename T>
	Vector3<T> BBox3<T>::getHalfSize() const
	{
		return getSize() * T(0.5);
	}


	template <class T>
	void BBox3<T>::setMax( const Vector3<T>& max )
	{
		for (int i = 0; i < 3; ++i)		
			m_max[i] = max[i];
	}

	template <class T>
	void BBox3<T>::setMin( const Vector3<T>& min )
	{
		//m_min = min;
		for (int i = 0; i < 3; ++i)
			m_min[i] = min[i];
	}

	template <class T>
	Vector3<T> BBox3<T>::getMax() const
	{
		return m_max.toVector3();
	}

	template <class T>
	Vector3<T> BBox3<T>::getMin() const
	{
		return m_min.toVector3();
	}

	template <class T>
	bool BBox3<T>::isValid() const
	{
		for( int i = 0; i < 3; ++i )
			if( m_min[i] >= m_max[i] )
				return false;
		return true;
	}
	template <typename T>
	Vector3<T> BBox3<T>::getSize() const
	{
		return (m_max - m_min).toVector3();
	}

	template <typename T>
    std::array<Vector3<T>, AABB_POINTS> BBox3<T>::toPoints() const
	{
        std::array<Vector3<T>, AABB_POINTS> ret;

        ret[0] = Math::Vector3<T>( m_max[0], m_max[1], m_min[2] );
        ret[1] = Math::Vector3<T>( m_max[0], m_min[1], m_min[2] );
        ret[2] = Math::Vector3<T>( m_min[0], m_min[1], m_min[2] );
        ret[3] = Math::Vector3<T>( m_min[0], m_max[1], m_min[2] );

        ret[4] = Math::Vector3<T>( m_max[0], m_max[1], m_max[2] );
        ret[5] = Math::Vector3<T>( m_max[0], m_min[1], m_max[2] );
        ret[6] = Math::Vector3<T>( m_min[0], m_min[1], m_max[2] );
        ret[7] = Math::Vector3<T>( m_min[0], m_max[1], m_max[2] );          

		return ret;
	}

	template <typename T>
	BBox3<T> BBox3<T>::expanded(const T& amount) const
	{
		BBox3<T> result = *this;
		result.expand(amount);
		return result;
	}


	template <typename T>
	BBox3<T>& BBox3<T>::operator+=( const BBox3<T>& other )
	{
		this->updateBounds( other );
		return *this;
	}


	template <typename T>
	BBox3<T>& BBox3<T>::operator+=(const Vector3<T>& translate )
	{
		m_min += translate.toVector4();
		m_max += translate.toVector4();
		return *this;
	}



    template <typename T>
    BBox3<T> BBox3<T>::scaled(const T& factor) const
    {
        auto copy = *this;
        copy.scale(factor);
        return copy;
    }

    template <typename T>
    void BBox3<T>::scale(const T& factor)
    {
        const auto center = getCenter();
        const auto offset = (getMax() - center) * factor;;
        
        for (int i = 0; i < 3; ++i) {
            m_min[i] = center[i] - offset[i];
            m_max[i] = center[i] + offset[i];
        }
    }


	template <class T>
	bool BBox3<T>::insideBounds( const Vector3<T>& point ) const
	{
		for( int i = 0; i < 3; ++i )
		{
			if( point[i] < m_min[i] || point[i] > m_max[i] )
				return false;
		}
		return true;
	}

	template <class T>
	bool BBox3<T>::intersects( const BBox3<T>& other ) const
	{
		for( int i = 0; i < 3; ++i )
		{
			if( other.m_min[i] > m_max[i] || other.m_max[i] < m_min[i] )
				return false;
		}
		return true;
	}


    template <typename T>
    bool BBox3<T>::isBigger(const T value) const
    {
        for (int i = 0; i < 3; ++i) {
            if (m_max[i] - m_min[i] < value)
                return false;
        }
        return true;
    }

    template <typename T>
    bool BBox3<T>::operator!=(const BBox3<T>& other) const
    {
        return !(*this == other);
    }

    template <typename T>
    bool BBox3<T>::operator==(const BBox3<T>& other) const
    {
        return m_min == other.m_min &&
            m_max == other.m_max;
    }

  
	

	template <class T>
	void BBox3<T>::expand(const T& amount)
	{
		for( int i = 0; i < 3; ++i )
		{
			m_min[i] -= amount;
			m_max[i] += amount;
		}
	}


	template <class T>
	void BBox3<T>::updateBounds( const BBox3<T>& other )
	{
		updateBounds( other.m_min.toVector3() );
		updateBounds( other.m_max.toVector3() );
	}

	template <class T>
	void BBox3<T>::updateBounds( const Vector3<T>& point )
	{
		for( int i = 0; i < 3; ++i )
		{
			if( point[i] < m_min[i] )
				m_min[i] = point[i];
			if( point[i] > m_max[i] )
				m_max[i] = point[i];
		}
	}

	template <class T>
	void BBox3<T>::clearBounds()
	{
		m_min[0] = m_min[1] = m_min[2] = std::numeric_limits<T>::max();
		m_max[0] = m_max[1] = m_max[2] = std::numeric_limits<T>::lowest();
	}

	template <typename T>
	Vector3<T> BBox3<T>::getCenter() const
	{
		return( m_min + ( m_max - m_min ) * (T)0.5 ).toVector3();
	}

	template <typename T>
	std::string BBox3<T>::toString() const
	{
        std::string buf;
        buf += m_min.toString() + " " + m_max.toString();
		return buf;
	}



	template <typename T>
	bool BBox3<T>::fromString(const std::string& val)
	{
		if( !m_min.fromString( val ) )
			return false;
		return m_max.fromString( val );
	}


    template <typename T>
    bool BBox3<T>::fromJsonObject(const IO::JSonObject& obj)
    {
        bool succeed = true;

        succeed &= m_min.fromJsonObject(obj["Min"]);
        succeed &= m_max.fromJsonObject(obj["Max"]);

        return succeed;
    }

    template <typename T>
    IO::JSonObject BBox3<T>::toJsonObject() const
    {
        IO::JSonObject result;
        result["Min"] = m_min.toJsonObject();
        result["Max"] = m_max.toJsonObject();
        return result;
    }

    /*namespace ns*/
   // {
        template <typename T>
        void to_json(IO::JSonObject& obj, const BBox2<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, BBox2<T>& v) {
            v.fromJsonObject(obj);
        }

        template <typename T>
        void to_json(IO::JSonObject& obj, const BBox3<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, BBox3<T>& v) {
            v.fromJsonObject(obj);
        }
   // }
   

};
#endif