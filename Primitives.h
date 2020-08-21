#ifndef PRIMITIVES_H
#define PRIMITIVES_H
#include "MathDecl.h"
#include "Vector.h"
#include "GenMath.h"  //quadratic solver
#include "Sphere.h"
namespace Math
{
	
	//////////////////////////////////////////////////////////////////////////
	//\brief: 2D Line class
	//////////////////////////////////////////////////////////////////////////
	template<class T>
	class Line2D
	{
	public:
		Line2D												(													);
		Line2D												( const Vector2<T>& p1, const Vector2<T>& p2		);
		~Line2D												(													);

		T						distanceToPoint				( const Vector2<T>& point							);

		void					set							( const Vector2<T>& p1, const Vector2<T>& p2		);
		const Vector2<T>&		getP1						(													) const;
		const Vector2<T>&		getP2						(													) const;
		Vector2<T>				getDirection				( bool normalize = true								) const;

	
		
	private:

		Vector2<T>				m_p1,
								m_p2;
	};
	
	//////////////////////////////////////////////////////////////////////////
	//\brief: Generic Line class
	//////////////////////////////////////////////////////////////////////////
	template<class Real, template<class> class Dim>
	class LineGeneric
	{
	public:
		LineGeneric											(														);	
		LineGeneric											( const Dim<Real>& p1, const Dim<Real>& p2				);
		~LineGeneric										(														);

		Real					distanceToPoint				( const Dim<Real>& point								);
		bool					closePointToRay				( const Dim<Real>& rayStart, const Dim<Real>& rayDir,	
															  Dim<Real>& pointOnLineOut, Dim<Real>& pointOnRayOut	) const;

		bool					segmentIntersection			( const LineGeneric<Real,Dim>& other, 
															  Dim<Real>& pointOnLineOut, Dim<Real>& pointOnRayOut	);

		void					set							( const Dim<Real>& p1, const Dim<Real>& p2				);
		const Dim<Real>&		getP1						(														) const;
		const Dim<Real>&		getP2						(														) const;
		Dim<Real>				getDirection				( bool normalize = true									) const;

		void					setSegment					( bool	val												);
		
	private:
		Dim<Real>				m_p1,
								m_p2;
		bool					m_isSegment;
	
	};

	
	typedef LineGeneric<float, Vector3> Line3D;
	typedef LineGeneric<float, Vector2> LineGeneric2D;		
	

	//////////////////////////////////////////////////////////////////////////
	// \brief: Basis cube class
	//////////////////////////////////////////////////////////////////////////
	template<class T>
	class Cuboid
	{
	public:
		Cuboid												( const Vector3<T>& center, T sideLength			);
		Cuboid												(													);
		~Cuboid												(													);

		int						intersect					( const Vector3<T>& rayPos, const Vector3<T>& rayDir,
															  Vector3<T>& p1Out, Vector3<T>& p2Out				);
		T                       getSideLength				(													) const;
		const Vector3<T>&		getCenter					(													) const;


		void					setSideLength				( const T& sideLength									);
		void					setCenter					( const Vector3<T>& center							);

	private:
		Vector3<T>				m_center;
		T						m_sideLength;

	};

	


	//////////////////////////////////////////////////////////////////////////
	// \brief: 2D Line implemnation
	//////////////////////////////////////////////////////////////////////////

	template<class T>
	Math::Line2D<T>::Line2D()
	{

	}
	
	
	template<class T>
	Math::Line2D<T>::Line2D(const Vector2<T>& p1, const Vector2<T>& p2) 
		: m_p1( p1 )
		, m_p2( p2 )
	{

	}
	
	
	template<class T>
	Math::Line2D<T>::~Line2D()
	{

	}

	template<class T>
	T Math::Line2D<T>::distanceToPoint(const Vector2<T>& point)
	{
		
		Vector2<T> pointToP1 = point - m_p1;
		Vector2<T> p2Top1    = m_p2 - m_p1;
		T dot = p2Top1.dot( pointToP1 );
		T lengthSqr = p2Top1.lengthSquared();
		T param = dot/lengthSqr;
		//T newX, newY;
		Vector2<T> newPoint;
		if( param < (T)0.0)
			newPoint = m_p1;
		else if( param > (T) 1.0 )
			newPoint =  m_p2;
		else 
			newPoint = m_p1 + ( p2Top1 * param );

		return newPoint.distance( point );	
	}

	template<class T>
	Vector2<T> Math::Line2D<T>::getDirection(bool normalize ) const
	{
		Vector2<T> dir = m_p2 - m_p1;
		T length =  dir.length();
		if( normalize && ( length > (T) 0.001 ) )
			dir /= length;
		return dir;
	}


	template<class T>
	const Vector2<T>& Math::Line2D<T>::getP2() const
	{
		return m_p2;
	}

	template<class T>
	const Vector2<T>& Math::Line2D<T>::getP1() const
	{
		return m_p1;
	}

	template<class T>
	void Math::Line2D<T>::set(const Vector2<T>& p1, const Vector2<T>& p2)
	{
		m_p1 = p1; 
		m_p2 = p2;;
	}


	

	//////////////////////////////////////////////////////////////////////////
	// \Brief: Cuboid Implementation
	//////////////////////////////////////////////////////////////////////////
	template<class T>
	Math::Cuboid<T>::Cuboid()
		: m_center( Vector3<T>( (T)0.0, (T)0.0, (T)0.0) )
		, m_sideLength((T) 1.0 )
	{

	}

	template<class T>
	Math::Cuboid<T>::Cuboid(const Vector3<T>& center, T sideLength)
		: m_center( center )
		, m_sideLength( sideLength )
	{

	}


	template<class T>
	Math::Cuboid<T>::~Cuboid()
	{

	}


	template<class T>
	void Math::Cuboid<T>::setCenter(const Vector3<T>& center)
	{
		m_center = center;
	}

	template<class T>
	void Math::Cuboid<T>::setSideLength(const T& sideLength)
	{
		m_sideLength = sideLength;
	}

	template<class T>
	const Vector3<T>& Math::Cuboid<T>::getCenter() const
	{
		return m_center;
	}

	template<class T>
	int Math::Cuboid<T>::intersect(const Vector3<T>& rayPos, const Vector3<T>& rayDir, Vector3<T>& p1Out, Vector3<T>& p2Out)
	{
		return 0;
	}



	//////////////////////////////////////////////////////////////////////////
	// \brief: Generic Line class
	//////////////////////////////////////////////////////////////////////////
	template<class Real, template<class> class Dim>
	Math::LineGeneric<Real, Dim>::LineGeneric()
		: m_isSegment( false )
	{

	}


	template<class Real, template<class> class Dim>
	Math::LineGeneric<Real, Dim>::LineGeneric(const Dim<Real>& p1, const Dim<Real>& p2)
		: m_p1 ( p1 )
		, m_p2 ( p2 )
		, m_isSegment( false )
	{

	}


	template<class Real, template<class> class Dim>
	Math::LineGeneric<Real, Dim>::~LineGeneric()
	{

	}


	template<class Real, template<class> class Dim>
	Dim<Real> Math::LineGeneric<Real, Dim>::getDirection(bool normalize /*= true */) const
	{
		Dim<Real> dir = m_p2 - m_p1;
		Real length =  dir.length();
		if( normalize && ( length > (Real) 0.0001 ) )
			dir /= length;
		return dir;
	}

	template<class Real, template<class> class Dim>
	const Dim<Real>& Math::LineGeneric<Real, Dim>::getP1() const
	{
		return m_p1;
	}

	template<class Real, template<class> class Dim>
	const Dim<Real>& Math::LineGeneric<Real, Dim>::getP2() const
	{
		return m_p2;
	}



	template<class Real, template<class> class Dim>
	void Math::LineGeneric<Real, Dim>::set(const Dim<Real>& p1, const Dim<Real>& p2)
	{
		m_p1  = p1;
		m_p2  = p2;
	}


	template<class Real, template<class> class Dim>
	Real Math::LineGeneric<Real, Dim>::distanceToPoint(const Dim<Real>& point)
	{
		Dim<Real> pointToP1 = point - m_p1;
		Dim<Real> p2Top1    = m_p2 - m_p1;
		Real dot = p2Top1.dot( pointToP1 );
		Real lengthSqr = p2Top1.lengthSquared();
		Real param = dot/lengthSqr;
		//T newX, newY;
		Dim<Real> newPoint;
		if( param < (Real)0.0)
			newPoint = m_p1;
		else if( param > (Real) 1.0 )
			newPoint =  m_p2;
		else 
			newPoint = m_p1 + ( p2Top1 * param );

		return newPoint.distance( point );	
	}


	template<class Real, template<class> class Dim>
	void Math::LineGeneric<Real, Dim>::setSegment(bool val)
	{
		m_isSegment = val;
	}


	template<class Real, template<class> class Dim>
	bool Math::LineGeneric<Real, Dim>::closePointToRay( const Dim<Real>& rayStart, const Dim<Real>& rayDir ,
		Dim<Real>& pointOnLineOut, Dim<Real>& pointOnRayOut) const
	{
		Dim<Real> d1	= getDirection( true );
		Dim<Real> d2    = rayDir;
		Real a, b, c, d, e, f;

		Dim<Real> r = m_p1 - rayStart;
		a = d1.dot( d1  );
		b = d1.dot( d2  );
		c = d1.dot( r   );
		e = d2.dot( d2  );
		f = d2.dot( r   );
		d = a*e -( b * b );

		if( d <= (Real) 0.0001)
			return false;

		Real s = (b*f - c*e)/d;
		Real t = (a*f - b*c)/d;

		if( m_isSegment )
		{
			Real lngth = m_p1.distance( m_p2 );
			s = Math::Clamp( (Real)0.0, (Real)lngth, s );
		}
		pointOnLineOut = m_p1 + ( d1 * s );
		pointOnRayOut =  rayStart + ( rayDir * t );



		return true;
	}

	template<class Real, template<class> class Dim>
	bool Math::LineGeneric<Real, Dim>::segmentIntersection(const LineGeneric<Real,Dim>& other, Dim<Real>& pointOnLineOut, Dim<Real>& thisLineOut)
	{
		Dim<Real> d1	= getDirection( true );
		Dim<Real> d2    = other.getDirection( true );
		Real a, b, c, d, e, f;

		Dim<Real> r = m_p1 - other.m_p1;
		a = d1.dot( d1  );
		b = d1.dot( d2  );
		c = d1.dot( r   );
		e = d2.dot( d2  );
		f = d2.dot( r   );
		d = a * e -( b * b ); //determinant

		if( d <= (Real) 0.0001 )
			return false;
		
		
		
		Real s = (b*f - c*e) / d;
		Real t = (a*f - b*c) / d;
		
		//s = Math::clamp((Real) 0.0, (Real) 1.0, s );
		//t = Math::clamp((Real) 0.0, (Real) 1.0, t );
		
		thisLineOut = m_p1 + ( d1 * s );
		pointOnLineOut = other.m_p1 + ( d2 * t );
		

		return true;

	}








	

}

#endif