#include "Quat.h"

namespace Math
{

    void to_json(IO::JSonObject& obj, const Quat<float>& v)
    {
        obj = v.toJsonObject();
    }

    void from_json(const IO::JSonObject& obj, Quat<float>& v)
    {
        v.fromJsonObject(obj);
    }
}

