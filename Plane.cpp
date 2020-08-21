#include "Plane.h"

namespace Math
{

    void to_json(IO::JSonObject& obj, const Plane<float>& v)
    {
        obj = v.toJsonObject();
    }

    void from_json(const IO::JSonObject& obj, Plane<float>& v)
    {
        v.fromJsonObject(obj);
    }

}

