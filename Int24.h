#pragma once

#include <limits>
#include "MathDecl.h"

#pragma pack(push)  
#pragma pack(1)   

struct Int24
{
    Int24()
    {
    }

    Int24(const int val)
    {
        *this = val;
    }

    Int24(const Int24& val)
    {
        *this = val;
    }

    operator int() const
    {
        if (m_data[2] & 0x80) // Is this a negative?  Then we need to siingn extend.
        {
            return (0xff << 24) | (m_data[2] << 16) | (m_data[1] << 8) | (m_data[0] << 0);
        }
        else
        {
            return (m_data[2] << 16) | (m_data[1] << 8) | (m_data[0] << 0);
        }
    }

    operator float() const
    {
        return (float)this->operator int();
    }

    Int24& operator =(const Int24& input)
    {
        m_data[0] = input.m_data[0];
        m_data[1] = input.m_data[1];
        m_data[2] = input.m_data[2];

        return *this;
    }

    Int24& operator =(const int input)
    {
        m_data[0] = ((unsigned char*)&input)[0];
        m_data[1] = ((unsigned char*)&input)[1];
        m_data[2] = ((unsigned char*)&input)[2];

        return *this;
    }

    /***********************************************/

    Int24 operator +(const Int24& val) const
    {
        return Int24((int)*this + (int)val);
    }

    Int24 operator -(const Int24& val) const
    {
        return Int24((int)*this - (int)val);
    }

    Int24 operator *(const Int24& val) const
    {
        return Int24((int)*this * (int)val);
    }

    Int24 operator /(const Int24& val) const
    {
        return Int24((int)*this / (int)val);
    }

    /***********************************************/

    Int24 operator +(const int val) const
    {
        return Int24((int)*this + val);
    }

    Int24 operator -(const int val) const
    {
        return Int24((int)*this - val);
    }

    Int24 operator *(const int val) const
    {
        return Int24((int)*this * val);
    }

    Int24 operator /(const int val) const
    {
        return Int24((int)*this / val);
    }

    /***********************************************/
    /***********************************************/


    Int24& operator +=(const Int24& val)
    {
        *this = *this + val;
        return *this;
    }

    Int24& operator -=(const Int24& val)
    {
        *this = *this - val;
        return *this;
    }

    Int24& operator *=(const Int24& val)
    {
        *this = *this * val;
        return *this;
    }

    Int24& operator /=(const Int24& val)
    {
        *this = *this / val;
        return *this;
    }

    /***********************************************/

    Int24& operator +=(const int val)
    {
        *this = *this + val;
        return *this;
    }

    Int24& operator -=(const int val)
    {
        *this = *this - val;
        return *this;
    }

    Int24& operator *=(const int val)
    {
        *this = *this * val;
        return *this;
    }

    Int24& operator /=(const int val)
    {
        *this = *this / val;
        return *this;
    }

    /***********************************************/
    /***********************************************/

    Int24 operator >>(const int val) const
    {
        return Int24((int)*this >> val);
    }

    Int24 operator <<(const int val) const
    {
        return Int24((int)*this << val);
    }

    /***********************************************/

    Int24& operator >>=(const int val)
    {
        *this = *this >> val;
        return *this;
    }

    Int24& operator <<=(const int val)
    {
        *this = *this << val;
        return *this;
    }

    /***********************************************/
    /***********************************************/

    operator bool() const
    {
        return (int)*this != 0;
    }

    bool operator !() const
    {
        return !((int)*this);
    }

    Int24 operator -()
    {
        return Int24(-(int)*this);
    }

    /***********************************************/
    /***********************************************/

    bool operator ==(const Int24& val) const
    {
        return (int)*this == (int)val;
    }

    bool operator !=(const Int24& val) const
    {
        return (int)*this != (int)val;
    }

    bool operator >=(const Int24& val) const
    {
        return (int)*this >= (int)val;
    }

    bool operator <=(const Int24& val) const
    {
        return (int)*this <= (int)val;
    }

    bool operator >(const Int24& val) const
    {
        return (int)*this > (int)val;
    }

    bool operator <(const Int24& val) const
    {
        return (int)*this < (int)val;
    }

    /***********************************************/

    bool operator ==(const int val) const
    {
        return (int)*this == val;
    }

    bool operator !=(const int val) const
    {
        return (int)*this != val;
    }

    bool operator >=(const int val) const
    {
        return (int)*this >= val;
    }

    bool operator <=(const int val) const
    {
        return (int)*this <= val;
    }

    bool operator >(const int val) const
    {
        return ((int)*this) > val;
    }

    bool operator <(const int val) const
    {
        return (int)*this < val;
    }


    char    m_data[3];
};
#pragma pack(pop) 


// CLASS numeric_limits<char>
template<> class std::numeric_limits<Int24>
    : public std::_Num_int_base
{
public:
    
public:
    _NODISCARD static constexpr int min()
    {	// return minimum value
        return - (1 << 23);
    }

    _NODISCARD static constexpr int max()
    {	// return minimum value
        return (1 << 23) - 1;
    }
};