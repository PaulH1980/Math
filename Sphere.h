#ifndef SPHERE_H
#define SPHERE_H
#include "MathDecl.h"
#include "AABB.h"
#include "Vector.h"

namespace Math
{
	template <typename T>
	class Sphere
	{
	public:
		Sphere();
		Sphere												( const Vector3<T>* points, std::uint32_t numPoint		);
		Sphere												( const Vector3<T>& pos, T radius					);

		std::uint32_t			intersect					( const Vector3<T>& rayPos, const Vector3<T>& rayDir,
															  Vector3<T>& p1Out, Vector3<T>& p2Out				) const;
		T                       getRadius					(													) const;
		const Vector3<T>&		getCenter					(													) const;

		void					setRadius					( const T& radius									);
		void					setCenter					( const Vector3<T>& center							);

		void					set							( const Vector3<T>* points, std::uint32_t numPoints, 
															 bool clear_ = false							    );

		BBox3<T>				toBounds					(													) const;
		void					clear						(													);
		T						distance					( const Vector3<T>& pos								) const;
		
		Vector3<T>				m_position;
		T						m_radius;
	};

	template <typename T>
	Sphere<T>::Sphere() : m_radius( (T) 0.0 )
	{

	}


	template <typename T>
	Sphere<T>::Sphere(const Vector3<T>* points, std::uint32_t numPoint)
	{
		set( points, numPoint );
	}


	template <typename T>
	Sphere<T>::Sphere(const Vector3<T>& pos, T radius) : m_position( pos )
		, m_radius( radius )
	{

	}


	template <typename T>
	void Math::Sphere<T>::setCenter(const Vector3<T>& center)
	{
		m_position = center;
	}

	template <typename T>
	void Math::Sphere<T>::setRadius(const T& radius)
	{
		m_radius = radius;
	}

	template <typename T>
	const Vector3<T>& Math::Sphere<T>::getCenter() const
	{
		return m_position;
	}

	template <typename T>
	T Math::Sphere<T>::getRadius() const
	{
		return m_radius;
	}
       

	template <typename T>
	void Sphere<T>::set(const Vector3<T>* points, std::uint32_t numPoints, bool clear_ /*= false */)
	{
		if( clear_ )
			clear();

		while( numPoints-- )
		{
			const auto& point = *points++;
            const auto offset = point - m_position;
			const auto dist = offset.length();

			if (dist > m_radius)
			{
				const auto half = (dist - m_radius) * static_cast<T>(0.5);
				m_radius += half;
				m_position += (half / dist) * offset;
			}
		}
	}


	template <typename T>
	BBox3<T> Sphere<T>::toBounds() const
	{
		auto min_ = m_position + Vector3<T>(-m_radius, -m_radius, -m_radius);
		auto max_ = m_position + Vector3<T>( m_radius,  m_radius,  m_radius);

		return BBox3<T>( min_, max_ );
	}


	template <typename T>
	void Sphere<T>::clear()
	{
		m_position = Vector3<T>::ZERO;
		m_radius   = static_cast<T>(0.0);
	}


	template <typename T>
	T Sphere<T>::distance(const Vector3<T>& pos) const
	{
		T dist = pos.distance( m_position );
		return dist < m_radius ? static_cast<T>(0.0) : dist - m_radius;
	}



	

	template <typename T>
	std::uint32_t Sphere<T>::intersect(const Vector3<T>& rayPos, const Vector3<T>& rayDir, Vector3<T>& p1Out, Vector3<T>& p2Out) const
	{
		Vector3<T> L = rayPos - m_position;
		const auto a = static_cast<T>(1.0);//ray.dir.dot(ray.dir);
		const auto b = rayDir.dot( L ) * static_cast<T>(2.0);
		const auto c = L.dot(L) - (m_radius * m_radius);
		T t1, t2; 
		const auto solutions = (std::uint32_t) SolveQuadratic( a, b, c, t1, t2 );	
		if( !solutions )
			return 0;
		p1Out = rayPos + rayDir * t1;
		p2Out = rayPos + rayDir * t2;
		return solutions;
	}

}

#endif