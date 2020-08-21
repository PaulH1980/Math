#pragma once

#include "Matrix.h"
#include "Vector.h"
#include "AxisAngle.h"
#include "AABB.h"
#include "Quat.h"
#include "EulerAngle.h"
#include "MathDefs.h"
#include "Plane.h"

using Vec2 = Math::Vector2<float>;
using Vec3 = Math::Vector3<float>;
using Vec4 = Math::Vector4<float>;
using Mat2 = Math::Matrix2<float>;
using Mat3 = Math::Matrix3<float>;
using Mat4 = Math::Matrix4<float>;

using Point2dVector = std::vector<Math::Vector2<float>>;
using Point3dVector = std::vector<Math::Vector3<float>>;
using PlaneVector	= std::vector<Math::Plane<float>>;
using IntVector		= std::vector<int>;
using UIntVector	= std::vector<std::uint32_t>;
using FloatVector	= std::vector<float>;

using HSV32F		= Math::Vector3<float>;
using RGB32F		= Math::Vector3<float>;
using RGBA32F		= Math::Vector4<float>;
using RGB8			= Math::Vector3<std::uint8_t>;
using RGBA8			= Math::Vector4<std::uint8_t>;

