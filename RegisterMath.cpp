#include <Common/ReflectionRegister.h>
#include <ScriptEngine/ScriptEngine.h>
#include <ScriptEngine/chaiscript.hpp>

#include "Vector.h"
#include "Quat.h"
#include "AABB.h"
#include "Range.h"
#include "EulerAngle.h"
#include "GenMath.h"
#include "RegisterMath.h"

using namespace Script;
using namespace chaiscript;


//non templated functions
void RegisterGenMath() {
    using namespace rttr;
    using namespace Math;
    //constants
    {
        registration::property_readonly("FLOAT_MIN", &FLOAT_MIN);
        registration::property_readonly("FLOAT_MAX", &FLOAT_MAX);

        registration::property_readonly("DOUBLE_MIN", &DOUBLE_MIN);
        registration::property_readonly("DOUBLE_MAX", &DOUBLE_MAX);
    }
    //methods
    {
        registration::method("NextPowerOfTwo", &NextPowerOfTwo);
        registration::method("PreviousPowerOfTwo", &PreviousPowerOfTwo);
        registration::method("NearestPowerOfTwo", &NearestPowerOfTwo);
        registration::method("BitCountBefore", &BitCountBefore);
        registration::method("HighestBit", &HighestBit);
        registration::method("PopCount8", &PopCount8);
        registration::method("PopCount32", &PopCount32);

        registration::method("NormalToOctahedron", &NormalToOctahedron);
        registration::method("OctahedronToNormal", &OctahedronToNormal);
        registration::method("CompressNormal", &CompressNormal);
        registration::method("RgbToHsv", &RgbToHsv);
        registration::method("HsvToRgb", &HsvToRgb);       
    }
}

template<typename T>
void RegisterIntGenMath(const std::string& suffix) {
    static_assert(std::is_integral<T>::value, "Not An Integral");
    using namespace rttr;
    using namespace Math;
    //int specializations
    registration::method("Gcd" + suffix, &GCD<T>);
    registration::method("IsOdd" + suffix, &IsOdd<T>);
}


template<typename T>
void RegisterFloatGenMath(const std::string& suffix) {
    static_assert(std::is_floating_point<T>::value, "Not A Floating Point");
    using namespace rttr;
    using namespace Math;
    registration::method("MinVal" + suffix, static_cast<T(*)(const T&, const T&)>(&MinVal<T>));
    registration::method("MaxVal" + suffix, static_cast<T(*)(const T&, const T&)>(&MaxVal<T>));
    registration::method("MaxVal" + suffix, static_cast<T(*)(const T&, const T&, const T&)>(&MaxVal<T>));
  
    registration::method("Log" + suffix, static_cast<T(*)( T, T )>(&Math::Log<T>));
    registration::method("Pi" + suffix, &Pi<T>);
    registration::method("HalfPi" + suffix, &HalfPi<T>);
    registration::method("TwoPi" + suffix, &TwoPi<T>);
    registration::method("Round" + suffix, &Round<T>);
    registration::method("Abs" + suffix, &Abs<T>);
    registration::method("Sign" + suffix, &Sign<T>);
   
    registration::method("IsNan" + suffix, &IsNaN<T>);
    registration::method("Sin" + suffix, &Sin<T>);
    registration::method("Cos" + suffix, &Cos<T>);
    registration::method("SinCos" + suffix, &SinCos<T>);
    registration::method("Clamp" + suffix, &Clamp<T>);
    registration::method("ToRadians" + suffix, &ToRadians<T>);
    registration::method("ToDegrees" + suffix, &ToDegrees<T>);
    registration::method("NormalizeRange" + suffix, &NormalizeRange<T>);
    registration::method("Lerp" + suffix, static_cast<T(*)(const std::vector<T>&, float)>(&Lerp<T>));
    registration::method("Lerp" + suffix, static_cast<T(*)(T, T, float)>(&Lerp<T>));
    registration::method("Atan2ToRads" + suffix, &Atan2ToRads<T>);

    registration::method("SnapToGrid" + suffix, &SnapToGrid<T>);
    registration::method("SnapPointToGrid" + suffix, &SnapPointToGrid<T>);
    registration::method("SnapScalarToGrid" + suffix, &SnapScalarToGrid<T>);    
    registration::method("SolveQuadratic" + suffix, &SolveQuadratic<T>);   
    registration::method("RandomValue" + suffix, &RandomValue<T>);
    registration::method("Random" + suffix, static_cast<T(*)(T, T)>(&Random<T>));
    registration::method("Random" + suffix, static_cast<T(*)()>(&Random<T>));
    registration::method("CalcAverage" + suffix, &CalcAverage<T>);
    registration::method("MakeRect" + suffix, static_cast<void(*)(int, int, int, int, std::array<Vector3<T>,4>& )>(&Math::MakeRect<T>));
    registration::method("MakeRect" + suffix, static_cast<void(*)( const Vector2<T>&, const Vector2<T>&, const T, std::array<Vector3<T>, 4>&)>(&Math::MakeRect<T>));
    registration::method("QuadToTriangles" + suffix, &QuadToTriangles<T>);
    registration::method("SignNotZero" + suffix, &SignNotZero<T>);
   
    registration::method("CalculateNormal" + suffix, &CalculateNormal<T>);
    registration::method("Equals" + suffix, &Equals<T>);
    registration::method("PointInPoly" + suffix, &PointInPoly<T>);
    registration::method("VerticesCounterClockWise" + suffix, &VerticesCounterClockWise<T>);
    registration::method("GetPolyCenter" + suffix, &GetPolyCenter<T>);
    registration::method("AABBInsideFrustum" + suffix, &AABBInsideFrustum<T>);
    registration::method("SortVerticesCCW" + suffix, &SortVerticesCCW<T>);
    registration::method("GetBaryCentric" + suffix, &GetBaryCentric<T>);
    registration::method("TriangleArea" + suffix,  static_cast<T(*)(const Vector3<T>&, const Vector3<T>&, const Vector3<T>&)>(&TriangleArea<T>));
    registration::method("IntersectsRay" + suffix, static_cast<bool(*)(const BBox3<T>&, const Vector3<T>&, const Vector3<T>&)>( &IntersectsRay<T>));
    registration::method("IntersectsRay" + suffix, static_cast<bool(*)(const Vector3<T>&, const Vector3<T>&, const Vector3<T>*, int, Vector3<T>&, size_t)>(&IntersectsRay<T>));

    registration::method("QuaternionToEulerAngle" + suffix, &QuaternionToEulerAngle<T>);
    registration::method("QuaternionToMatrix3" + suffix, &QuaternionToMatrix3<T>);
    registration::method("QuaternionToMatrix4" + suffix, &QuaternionToMatrix4<T>);
    registration::method("QuaternionToAxisAngle" + suffix, &QuaternionToAxisAngle<T>);

    registration::method("AxisAngleToMatrix3" + suffix, &AxisAngleToMatrix3<T>);
    registration::method("AxisAngleToMatrix4" + suffix, &AxisAngleToMatrix4<T>);
    registration::method("AxisAngleToEulerAngle" + suffix, &AxisAngleToEulerAngle<T>);
    registration::method("AxisAngleToQuaternion" + suffix, &AxisAngleToQuaternion<T>);

    registration::method("Matrix3ToEulerAngle" + suffix, &Matrix3ToEulerAngle<T>);;
    registration::method("Matrix3ToMatrix4" + suffix, &Matrix3ToMatrix4<T>);
    registration::method("Matrix3ToQuaternion" + suffix, &Matrix3ToQuaternion<T>);
    registration::method("Matrix3ToAxisAngle" + suffix, &Matrix3ToAxisAngle<T>);

    registration::method("Matrix4ToEulerAngle" + suffix,    &Matrix4ToEulerAngle<T>);
    registration::method("Matrix4ToMatrix3" + suffix,       &Matrix4ToMatrix3<T>);
    registration::method("Matrix4ToQuaternion" + suffix,    &Matrix4ToQuaternion<T>);
    registration::method("Matrix4ToAxisAngle" + suffix,     &Matrix4ToAxisAngle<T>);

    registration::method("LookAtMatrix" + suffix, &LookAtMatrix<T>);
    registration::method("GenerateLookAtMatrices" + suffix, &GenerateLookAtMatrices<T>);
    registration::method("Reflect" + suffix, &Reflect<T>);
    registration::method("BranchlessONB" + suffix, &BranchlessONB<T>);
    registration::method("ProjectPointOnLine" + suffix, &ProjectPointOnLine<T>);

}

template<typename T>
void RegisterAxisAngle( const std::string& className ) {
    static_assert(std::is_floating_point<T>::value, "Not A Floating Point");
    using namespace rttr;
    using namespace Math;
    using Type = AxisAngle<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<const Vector3<T>&, T>();
    val.constructor<T, T, T, T>();
    val.constructor<const Type&>();
    val.method("setAxis", &Type::setAxis);
    val.method("setAngle", &Type::setAngle);
    val.method("getAxis", &Type::getAxis);
    val.method("getAngle", &Type::getAngle);   

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    val.property("Axis", &Type::m_axis);
    val.property("Angle", &Type::m_angle);
};

template<typename T>
void RegisterEulerAngle(const std::string& className) {
    static_assert(std::is_floating_point<T>::value, "Not A Floating Point");
    using namespace rttr;
    using namespace Math;

    using Type = EulerAngle<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<const Vector3<T>&>();
    val.constructor<T, T, T>();
    val.constructor<const Type&>();

    val.method("setPitch", &Type::setPitch);
    val.method("setYaw", &Type::setYaw);
    val.method("setRoll", &Type::setRoll);
    val.method("setPitchYawRoll", static_cast<void(Type::*)(const Vector3<T>&)>(&Type::setPitchYawRoll));
    val.method("setPitchYawRoll", static_cast<void(Type::*)(T, T, T)>(&Type::setPitchYawRoll));

    val.method("getPitch", &Type::getPitch);
    val.method("getYaw", &Type::getYaw);
    val.method("getRoll", &Type::getRoll);
    val.method("getPitchYawRoll", &Type::getPitchYawRoll);  

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    val.property("PitchYawRoll", &Type::m_pry);
};

template<typename T>
void RegisterQuaternion(const std::string& className) {
    static_assert(std::is_floating_point<T>::value, "Not A Floating Point");
    using namespace rttr;
    using namespace Math;
    using Type = Quat<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();    
    val.constructor<T, T, T, T>();
    val.constructor<const Type&>();
    val.method("identity", &Type::identity);
    val.method("normalized", &Type::normalize);
    val.method("length", &Type::length);
    val.method("getNormalized", &Type::getNormalized);
    val.method("setXYZW", &Type::setXYZW);   
    val.method("slerp", &Type::slerp);
    val.method("dot", &Type::dot);
    val.method("conjugate", &Type::conjugate);
    val.method("conjugated", &Type::conjugated);
    val.method("invert", &Type::invert);
    val.method("inverted", &Type::inverted);
    val.method("transformPoint", &Type::transformPoint);
    val.method("transformNormal", &Type::transformNormal);
    val.method("multiply", &Type::multiply);
    val.method("multiplied", &Type::multiplied);  
  
    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);
};

template<typename T>
void RegisterClampedRange(const std::string& className) {
    using namespace rttr;
    using namespace Math;

    using Type = ClampedValue<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T, T>();
    val.constructor<T, T, T>();
    val.constructor<const Type&>();

    val.method("getMin", &Type::getMin);
    val.method("getMax", &Type::getMax);
    val.method("setMin", &Type::setMin);
    val.method("setMax", &Type::setMax);
    val.method("setValue", &Type::setValue);
    val.method("getValue", &Type::getValue);
    val.method("setFromPercentage", &Type::setFromPercentage);
};

template<typename T>
void RegisterRange(const std::string& className) {
    using namespace rttr;
    using namespace Math;
    using Type = Range<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T, T>();
    val.constructor<const Type&>();

    val.method("set", &Type::set);
    val.method("setMin", &Type::setMin);
    val.method("setMax", &Type::setMax);

    val.method("getMin", &Type::getMin);
    val.method("getMax", &Type::getMax);

    val.method("update", &Type::update);
    val.method("valid", &Type::valid);
    val.method("reset", &Type::reset);
    val.method("intersects", &Type::intersects);
};

template<typename T>
void RegisterPlane(const std::string& className) {
    static_assert(std::is_floating_point<T>::value, "Not A Floating Point");
    
    using namespace rttr;
    using namespace Math;
    using Type = Plane<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T, T, T, T>();
    val.constructor<const Vector3<T>&, T>();
    val.constructor<const Vector3<T>&, const Vector3<T>&>();
    val.constructor<const Vector3<T>&, const Vector3<T>&, const Vector3<T>>();
    val.constructor<const Vector4<T>&>();
    val.constructor<const Type&>();

    val.method("set", &Type::set);
    val.method("normalize", &Type::normalize);
    val.method("setDistance", &Type::setDistance);
    val.method("setNormal", &Type::setNormal);
    val.method("getDistance", &Type::getDistance);
    val.method("getNormal", &Type::getNormal);


    val.method("intersect", &Type::intersect);
    val.method("distanceToPlane", &Type::distanceToPlane);
    val.method("projectPointOnPlane", &Type::projectPointOnPlane);

    val.method("classifySphere", &Type::classifySphere);
    val.method("classifyPoint", &Type::classifyPoint);
    val.method("invert", &Type::invert);
    val.method("inverted", &Type::inverted);

    val.method("toVector4", &Type::toVector4);
    val.method("fromVector4", &Type::fromVector4);
    val.method("reflectionMatrix", &Type::reflectionMatrix);
    
    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    val.property("PlaneDistance", &Type::m_distance);
    val.property("PlaneNormal", &Type::m_normal);
};

template<typename T>
void RegisterBBox2(const std::string& className) {
    using namespace rttr;
    using namespace Math;
    using Type = BBox2<T>;
    auto val = registration::class_<Type>(className);
    
    val.constructor();
    val.constructor<const Vector2<T>&, const Vector2<T>&>();
    val.constructor<const T, T, T, T>();
    val.constructor<const Type&>();

    val.method("getMin", &Type::getMin);
    val.method("getMax", &Type::getMax);
    val.method("setMin", &Type::setMin);
    val.method("setMax", &Type::setMax);
    val.method("isValid", &Type::isValid);
    val.method("expand", &Type::expand);
    val.method("expanded", &Type::expanded);
    val.method("equals", &Type::equals);
    val.method("getSize", &Type::getSize);
    val.method("getHalfSize", &Type::getHalfSize);
    val.method("getCenter", &Type::getCenter);
    val.method("clearBounds", &Type::clearBounds);
    val.method("insideBounds", &Type::insideBounds);
    val.method("intersects", &Type::intersects);

    val.method("updateBounds", static_cast<void(Type::*)(const Vector2<T>&)>(&Type::updateBounds));
    val.method("updateBounds", static_cast<void(Type::*)(const Type&)>(&Type::updateBounds));   
      
    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

};

template<typename T>
void RegisterBBox3(const std::string& className) {
    using namespace rttr;
    using namespace Math;

    using namespace rttr;
    using namespace Math;
    using Type = BBox3<T>;
    auto val = registration::class_<Type>(className);

    val.constructor();
    val.constructor<const Vector3<T>&, const Vector3<T>&>();
    val.constructor<const T, T, T, T, T, T>();
    val.constructor<const Type&>();

    val.method("getMin", &Type::getMin);
    val.method("getMax", &Type::getMax);
    val.method("setMin", &Type::setMin);
    val.method("setMax", &Type::setMax);
    val.method("isValid", &Type::isValid);
    val.method("expand", &Type::expand);
    val.method("expanded", &Type::expanded);
    val.method("equals", &Type::equals);
    val.method("getSize", &Type::getSize);
    val.method("getHalfSize", &Type::getHalfSize);
    val.method("getCenter", &Type::getCenter);
    val.method("clearBounds", &Type::clearBounds);
    val.method("insideBounds", &Type::insideBounds);
    val.method("intersects", &Type::intersects);

    val.method("updateBounds", static_cast<void(Type::*)(const Vector3<T>&)>(&Type::updateBounds));
    val.method("updateBounds", static_cast<void(Type::*)(const Type&)>(&Type::updateBounds));
       
    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);
};


template<typename T>
void RegisterVector2(const std::string& className) {
    using namespace rttr;
    using namespace Math;
    using Type = Vector2<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T>();
    val.constructor<T, T>();
    val.constructor<const T*>();
    val.constructor<const Type&>();

    val.method("get", &Type::get);
    val.method("getX", &Type::getX);
    val.method("getY", &Type::getY);

    val.method("set",  &Type::set);
    val.method("setX", &Type::setX);
    val.method("setY", &Type::setY);
    val.method("setXY", static_cast<void(Type::*)(T, T)>(&Type::setXY));
    val.method("setXY", static_cast<void(Type::*)(T)>(&Type::setXY));
    val.method("setXY", static_cast<void(Type::*)(const T*)>(&Type::setXY));   

    val.method("getNormalized", &Type::getNormalized);
    val.method("getAbsolute", &Type::getAbsolute);
    val.method("inverted", &Type::inverted);
    val.method("clear", &Type::clear);
    val.method("dot", &Type::dot);
    val.method("cross", &Type::cross);
    val.method("perpendicular", &Type::perpendicular);
    val.method("length", &Type::length);
    val.method("lengthSquared", &Type::lengthSquared);
    val.method("distance", &Type::distance);
    val.method("distanceSquared", &Type::distanceSquared);
    val.method("normalize", &Type::normalize);
    val.method("equals", &Type::equals);
    val.method("lerp", &Type::lerp);
    val.method("getMaxValue", &Type::getMaxValue);
    val.method("getMinValue", &Type::getMinValue);
    val.method("toVector3", &Type::toVector3);
    val.method("randomVector", &Type::randomVector);
    val.method("getPNorm", &Type::getPNormVector);

    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);
  
    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    //operators
    val.method("plus",   static_cast<Type(Type::*)(const Type&) const> (&Type::operator+));
    val.method("minus",  static_cast<Type(Type::*)(const Type&) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(const Type&) const> (&Type::operator/));
    val.method("times",  static_cast<Type(Type::*)(const Type&) const> (&Type::operator*));

    val.method("plus", static_cast<Type(Type::*)(T) const> (&Type::operator+));
    val.method("minus", static_cast<Type(Type::*)(T) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(T) const> (&Type::operator/));
    val.method("times", static_cast<Type(Type::*)(T) const> (&Type::operator*));
    
    val.method("plusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator*=));

    val.method("plusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(T)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(T)> (&Type::operator*=));

    val.method("negate", static_cast<void(Type::*)()> (&Type::operator-));
    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);
    val.method("==", &Type::operator==);
    val.method("!=", &Type::operator!=);

};

template<typename T>
void RegisterVector3(const std::string& className) {
    using namespace rttr;
    using namespace Math;
    using Type = Vector3<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T>();
    val.constructor<const T*>();
    val.constructor<T, T, T>();
    val.constructor<const Type&>();

    val.method("get", &Type::get);
    val.method("getX", &Type::getX);
    val.method("getY", &Type::getY);
    val.method("getZ", &Type::getZ);

    val.method("set", &Type::set);
    val.method("setX", &Type::setX);
    val.method("setY", &Type::setY);
    val.method("setZ", &Type::setZ);
    val.method("setXYZ", static_cast<void(Type::*)(T, T, T)>(&Type::setXYZ));
    val.method("setXYZ", static_cast<void(Type::*)(T)>(&Type::setXYZ));
    val.method("setXYZ", static_cast<void(Type::*)(const T*)>(&Type::setXYZ));

    val.method("getNormalized", &Type::getNormalized);
    val.method("getAbsolute", &Type::getAbsolute);
    val.method("inverted", &Type::inverted);
    val.method("clear", &Type::clear);
    val.method("dot", &Type::dot);
    val.method("cross", &Type::cross);    
    val.method("length", &Type::length);
    val.method("lengthSquared", &Type::lengthSquared);
    val.method("distance", &Type::distance);
    val.method("distanceSquared", &Type::distanceSquared);
    val.method("normalize", &Type::normalize);
    val.method("equals", &Type::equals);
    val.method("lerp", &Type::lerp);
    val.method("getMaxValue", &Type::getMaxValue);
    val.method("getMinValue", &Type::getMinValue);
    val.method("toVector2", &Type::toVector2);
    val.method("toVector4", &Type::toVector4);
    val.method("randomVector", &Type::randomVector);
    val.method("getPNorm", &Type::getPNormVector);
    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    //operators
    val.method("plus", static_cast<Type(Type::*)(const Type&) const> (&Type::operator+));
    val.method("minus", static_cast<Type(Type::*)(const Type&) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(const Type&) const> (&Type::operator/));
    val.method("times", static_cast<Type(Type::*)(const Type&) const> (&Type::operator*));

    val.method("plus", static_cast<Type(Type::*)(T) const> (&Type::operator+));
    val.method("minus", static_cast<Type(Type::*)(T) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(T) const> (&Type::operator/));
    val.method("times", static_cast<Type(Type::*)(T) const> (&Type::operator*));

    val.method("plusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator*=));

    val.method("plusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(T)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(T)> (&Type::operator*=));

    val.method("negate", static_cast<void(Type::*)()> (&Type::operator-));
    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);
    val.method("==", &Type::operator==);
    val.method("!=", &Type::operator!=);
};

template<typename T>
void RegisterVector4(const std::string& className) {
    using namespace rttr;
    using namespace Math;
    using Type = Vector4<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<const T*>();
    val.constructor<T>();
    val.constructor<T, T, T, T>();
    val.constructor<const Vector2<T>&, const Vector2<T>&>();
    val.constructor<const Vector3<T>&, T>();   
    val.constructor<const Type&>();

    val.method("get", &Type::get);
    val.method("getX", &Type::getX);
    val.method("getY", &Type::getY);
    val.method("getZ", &Type::getZ);
    val.method("getW", &Type::getW);

    val.method("set", &Type::set);
    val.method("setX", &Type::setX);
    val.method("setY", &Type::setY);
    val.method("setZ", &Type::setZ);
    val.method("setW", &Type::setW);
    val.method("setXYZW", static_cast<void(Type::*)(T, T, T, T)>(&Type::setXYZW));
    val.method("setXYZW", static_cast<void(Type::*)(T)>(&Type::setXYZW));
    val.method("setXYZW", static_cast<void(Type::*)(const T*)>(&Type::setXYZW));

    val.method("getNormalized", &Type::getNormalized);
    val.method("getAbsolute", &Type::getAbsolute);
    val.method("inverted", &Type::inverted);
    val.method("clear", &Type::clear);
    val.method("dot", &Type::dot);
    val.method("length", &Type::length);
    val.method("lengthSquared", &Type::lengthSquared);
    val.method("distance", &Type::distance);
    val.method("distanceSquared", &Type::distanceSquared);
    val.method("normalize", &Type::normalize);
    val.method("equals", &Type::equals);
    val.method("lerp", &Type::lerp);
    val.method("getMaxValue", &Type::getMaxValue);
    val.method("getMinValue", &Type::getMinValue);
    val.method("toVector3", &Type::toVector3);   
    val.method("randomVector", &Type::randomVector);
    val.method("getPNorm", &Type::getPNormVector);

    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

    //operators
    val.method("plus", static_cast<Type(Type::*)(const Type&) const> (&Type::operator+));
    val.method("minus", static_cast<Type(Type::*)(const Type&) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(const Type&) const> (&Type::operator/));
    val.method("times", static_cast<Type(Type::*)(const Type&) const> (&Type::operator*));

    val.method("plus", static_cast<Type(Type::*)(T) const> (&Type::operator+));
    val.method("minus", static_cast<Type(Type::*)(T) const> (&Type::operator-));
    val.method("divide", static_cast<Type(Type::*)(T) const> (&Type::operator/));
    val.method("times", static_cast<Type(Type::*)(T) const> (&Type::operator*));

    val.method("plusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(const Type&)> (&Type::operator*=));

    val.method("plusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator+=));
    val.method("minusAssign", static_cast<Type&(Type::*)(T)> (&Type::operator-=));
    val.method("divideAssign", static_cast<Type&(Type::*)(T)> (&Type::operator/=));
    val.method("timesAssign", static_cast<Type&(Type::*)(T)> (&Type::operator*=));

    val.method("negate", static_cast<void(Type::*)()> (&Type::operator-));
    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);
    val.method("==", &Type::operator==);
    val.method("!=", &Type::operator!=);
};

template<typename T>
void RegisterVectorX(const std::string& className) {
    using namespace rttr;
    using namespace Math;
};


template<typename T>
void RegisterMatrix2(const std::string& className) {
    using namespace rttr;
    using namespace Math;

    using Type = Matrix2<T>;
    auto val = registration::class_<Type>(className);
    val.constructor();
    val.constructor<T, T, T, T>();
    val.constructor<const Vector2<T>&, const Vector2<T>&>();
    val.constructor<const Type&>();

    val.method("setIdentity", &Type::setIdentity);
    val.method("zero", &Type::zero);
    val.method("invert", &Type::invert);
    val.method("inverted", &Type::inverted);
    val.method("isIdentity", &Type::isIdentity);
    val.method("rotationMatrix", &Type::rotationMatrix);
    val.method("determinant", &Type::determinant);
    val.method("trace", &Type::trace);
    val.method("transpose", static_cast<void(Type::*)(Type&) const>(&Type::transpose));
    val.method("transpose", static_cast<void(Type::*)()>(&Type::transpose));
    val.method("scale", static_cast<void(Type::*)(T)>(&Type::scale));
    val.method("equals", &Type::equals);
    val.method("multiply", &Type::multiply);
    val.method("postMultiply", &Type::postMultiply);
    val.method("preMultiply", &Type::preMultiply);
    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);
    val.method("toMatrix3", &Type::toMatrix3);

    val.method("setElement", static_cast<void(Type::*)(int, T)> (&Type::setElement));
    val.method("setElement", static_cast<void(Type::*)( int, int, T)> (&Type::setElement));

    val.method("getElement", static_cast<T(Type::*)(int)const> (&Type::getElement));
    val.method("getElement", static_cast<T(Type::*)(int, int)const> (&Type::getElement));

    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);
};

template<typename T>
void RegisterMatrix3(const std::string& className) {
    using namespace rttr;
    using namespace Math;

    using Type = Matrix3<T>;
    auto val = registration::class_<Type>(className);

    val.constructor();
    val.constructor<T, T, T, 
                    T, T, T,
                    T, T, T>();
    val.constructor<const Matrix2<T>&>();
    val.constructor<const Vector3<T>&, const Vector3<T>&, const Vector3<T>&>();
    val.constructor<const Type&>();

    val.method("setIdentity", &Type::setIdentity);
    val.method("zero", &Type::zero);
    val.method("invert", &Type::invert);
    val.method("inverted", &Type::inverted);
    val.method("isIdentity", &Type::isIdentity);    
    val.method("determinant", &Type::determinant);
    val.method("trace", &Type::trace);
    val.method("transpose", static_cast<void(Type::*)(Type&) const>(&Type::transpose));
    val.method("transpose", static_cast<void(Type::*)()>(&Type::transpose));
    val.method("transposed", &Type::transposed);
    val.method("scale", static_cast<void(Type::*)(T)>(&Type::scale));
    val.method("equals", &Type::equals);
    val.method("multiply", &Type::multiply);
    val.method("postMultiply", &Type::postMultiply);
    val.method("preMultiply", &Type::preMultiply);
    val.method("adjoint", &Type::adjoint);
    val.method("rotationMatrixX", &Type::rotationMatrixX);
    val.method("rotationMatrixY", &Type::rotationMatrixY);
    val.method("rotationMatrixZ", &Type::rotationMatrixZ);
    val.method("adjoint", &Type::adjoint);
    val.method("toMatrix2", &Type::toMatrix2);
    val.method("toMatrix4", &Type::toMatrix4);

    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);

    val.method("setElement", static_cast<void(Type::*)(int, T)> (&Type::setElement));
    val.method("setElement", static_cast<void(Type::*)(int, int, T)> (&Type::setElement));

    val.method("getElement", static_cast<T(Type::*)(int)const> (&Type::getElement));
    val.method("getElement", static_cast<T(Type::*)(int, int)const> (&Type::getElement));

    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);
};

template<typename T>
void RegisterMatrix4( const std::string& className) {
    using namespace rttr;
    using namespace Math;

    using Type = Matrix4<T>;
    auto val = registration::class_<Type>(className);

    val.constructor();
    val.constructor<T, T, T, T,
                    T, T, T, T,
                    T, T, T, T,
                    T, T, T, T>();
    val.constructor<const Matrix3<T>&>();
    val.constructor<const Vector4<T>&, const Vector4<T>&, const Vector4<T>&, const Vector4<T>&>();
    val.constructor<const Type&>();

    val.method("setIdentity", &Type::setIdentity);
    val.method("zero", &Type::zero);
    val.method("invert", &Type::invert);
    val.method("inverted", &Type::inverted);
    val.method("isIdentity", &Type::isIdentity);   
    val.method("determinant", &Type::determinant);
    val.method("trace", &Type::trace);
    val.method("transpose", static_cast<void(Type::*)(Type&) const>(&Type::transpose));
    val.method("transpose", static_cast<void(Type::*)()>(&Type::transpose));
    val.method("transposed", &Type::transposed);
    val.method("scale", static_cast<void(Type::*)(T)>(&Type::scale));
    val.method("equals", &Type::equals);
    val.method("postMultiply", &Type::postMultiply);
    val.method("preMultiply", &Type::preMultiply);
    val.method("adjoint", &Type::adjoint);
    val.method("rotationMatrixX", &Type::rotationMatrixX);
    val.method("rotationMatrixY", &Type::rotationMatrixY);
    val.method("rotationMatrixZ", &Type::rotationMatrixZ);
    val.method("adjoint", &Type::adjoint);
    val.method("toMatrix3", &Type::toMatrix3);  
    val.method("getDirection", &Type::getDirection);
    val.method("getRight", &Type::getRight);
    val.method("getUp", &Type::getUp);
    val.method("getScale", &Type::getScale);
    val.method("getTranslation", &Type::getTranslation);
    val.method("setDirection", &Type::setDirection);
    val.method("setRight", &Type::setRight);
    val.method("setUp", &Type::setUp);
    val.method("setTranslation", &Type::setTranslation);

    val.method("toPointer", &Type::toPointer);
    val.method("toConstPointer", &Type::toConstPointer);

    val.method("multiply", static_cast<Vector3<T>(Type::*)(const Vector3<T>&) const>(&Type::multiply));
    val.method("multiply", static_cast<Vector4<T>(Type::*)(const Vector4<T>&) const>(&Type::multiply));
    val.method("transformNormal", &Type::transformNormal);
    val.method("lookAt", &Type::lookAt);
    val.method("perspective", &Type::perspectiveMatrix);
    val.method("ortho2d", &Type::orthographicMatrix2d);
    val.method("ortho", &Type::orthographicMatrix);

    val.method("setElement", static_cast<void(Type::*)(int, T)> (&Type::setElement));
    val.method("setElement", static_cast<void(Type::*)(int, int, T)> (&Type::setElement));

    val.method("getElement", static_cast<T(Type::*)(int)const> (&Type::getElement));
    val.method("getElement", static_cast<T(Type::*)(int, int)const> (&Type::getElement));

    val.method("assign", &Type::operator=);
    val.method("=", &Type::operator=);

    val.method("toString", &Type::toString);
    val.method("fromString", &Type::fromString);

};

template<typename T>
void RegisterMatrixX(const std::string& className) {
    using namespace rttr;
    using namespace Math;
};


RTTR_REGISTRATION
{
    using namespace rttr;
    using namespace Math;
    {
        std::cout << "Hello Math Registration\n";
        RegisterGenMath();
        RegisterFloatGenMath<float>("_fp32");
        
        //floating point
        RegisterAxisAngle<float>("AxisAnglef");
        RegisterEulerAngle<float>("EulerAnglef");
        RegisterQuaternion<float>("Quatf");
        RegisterClampedRange<float>("ClampedRangef");
        RegisterRange<float>("Rangef");
        RegisterPlane<float>("Planef");
        RegisterBBox2<float>("BBox2f");
        RegisterBBox3<float>("BBox3f");
        RegisterVector2<float>("Vector2f");
        RegisterVector3<float>("Vector3f");
        RegisterVector4<float>("Vector4f");
        RegisterMatrix2<float>("Matrix2f");
        RegisterMatrix3<float>("Matrix3f");
        RegisterMatrix4<float>("Matrix4f");
        //ui
        RegisterClampedRange<std::uint32_t>("ClampedRangeui");
        RegisterRange<std::uint32_t>("Rangeui");
        //int
        RegisterIntGenMath<int>("_i32");
        RegisterClampedRange<int>("ClampedRangei");
        RegisterRange<int>("Rangei");
        RegisterBBox2<int>("BBox2i");
        RegisterBBox3<int>("BBox3i");
    }
};

namespace Math
{
   
    void Register(Script::ScriptEngine* se)
    {
        using namespace rttr;

        std::vector<std::string> registeredMathTypes =
        {
            "AxisAnglef",
            "EulerAnglef",
            "Quatf",
            "ClampedRangef",
            "Rangef",
            "Planef",
            "BBox2f",
            "BBox3f",
            "Vector2f",
            "Vector3f",
            "Vector4f",
            "Matrix2f",
            "Matrix3f",
            "Matrix4f",

            "ClampedRangeui",
            "Rangeui",

            "ClampedRangei",
            "Rangei",
            "BBox2i",
            "BBox3i"
        };

        bool succeed = true;
        for (const auto& typeName : registeredMathTypes) {
            const auto reflType = type::get_by_name(typeName);
            if (reflType.is_class()) {
                succeed &= se->registerTypeToScript(reflType);
            }
            else {
                se->addLogMessage( "Unknown Type: " + typeName, Log::LOG_LEVEL_WARNING);
            }
        }

    }

}

