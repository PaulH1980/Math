#pragma once
#include "MathDecl.h"
#include "GenMath.h"

namespace Math
{

	template<class Real, template<class> class Dim>
	static inline Dim<Real> EvalBezier(const Dim<Real>& p1, const Dim<Real>& p2, const Dim<Real>& p3, Real t)
	{
		Real u = (Real)1.0 - t;
		Real uu = u * u;
		Real tt = t *t;
		Real u2 = (Real) 2.0 * u;
		Real tu2 = t * u2;
		Dim<Real> res = (p1 * uu) + (p2 * tu2) + (p3 * tt);
		return res;
	}

	template<class Real, template<class> class Dim>
	static inline Dim<Real> EvalBezierTangent(const Dim<Real>& p1, const Dim<Real>& p2, const Dim<Real>& p3, Real t)
	{
		Real u = (Real) 1.0 - t;
		Real u2 = (Real) 2.0 * u;   //(1.0 - t)^2 --> 2*(1.0 - t) //first coeff
		Real t2 = (Real) 2.0 * t;
		Dim<Real> res = (p2 - p1) * u2 + (p3 - p2) * t2;
		return res;
	}

	template<class Real, template<class> class Dim>
	static inline Dim<Real> EvalBezier(const Dim<Real>& p1, const Dim<Real>& p2, const Dim<Real>& p3, const Dim<Real>& p4, Real t)
	{

		Real u = (Real)1.0 - t;
		Real tt = t*t;
		Real uu = u*u;
		Real uuu = uu * u;
		Real ttt = tt * t;

		Dim<Real> p = p1 * uuu; //first term
		p += p2 * ((Real) 3.0 * uu * t); //second term
		p += p3 * ((Real) 3.0 * u * tt); //third term
		p += p4 * ttt; //fourth term

		return p;
	}

}