#pragma once
#include <string>
#include <limits>
#include <cmath>
#include <IO/JSonObject.h>
#include "MathDecl.h"

namespace Math
{
	template<typename T>
	class Range
	{
	public:

		Range() 
		{
			reset();
		}
			
		Range( T min, T max)
			: m_minValue(min)
			, m_maxValue(max)
		{}

		Range( const Range<T>& rhs )
			: m_minValue(rhs.m_minValue)
			, m_maxValue(rhs.m_maxValue)
		{}


		Range&		operator = (const Range& rhs)
		{
			m_maxValue = rhs.m_maxValue;
			m_minValue = rhs.m_minValue;
			return *this;
		}

		bool		operator != (const Range& rhs) const
		{
			return  !(*this == rhs);				
		}

		bool		operator == (const Range& rhs) const
		{
			return  m_maxValue == rhs.m_maxValue &&
					m_minValue == rhs.m_minValue;
		}

		bool					intersects(const Range& rhs)
		{
			if (rhs.m_maxValue < m_minValue)
				return false;
			if (rhs.m_maxValue > m_minValue)
				return false;
			return true;
		}

		void					setMin(T val) { m_minValue = val; }
		void					setMax(T val) { m_maxValue = val; }

		void					set(T min, T max)
		{
			m_maxValue = max;
			m_minValue = min;
		}

		bool					valid() const 
		{
			return m_maxValue > m_minValue;
		}

		void					reset()
		{
			m_minValue = static_cast<T>((std::numeric_limits<T>::max)());
			m_maxValue = static_cast<T>((std::numeric_limits<T>::lowest()));
		}

		T						getMin() const { return m_minValue; }
		T						getMax() const { return m_maxValue; }

		std::string				toString() const
		{
			char buf[256] = {'\0'};
#pragma warning( push )
#pragma  warning (disable:4996)
			sprintf(buf, "%g %g", (double)m_minValue, (double) m_maxValue );
#pragma warning( pop )
			return std::string(buf);
		}

		void					update(T val)
		{
			if (val < m_minValue)
				m_minValue = val;
			if (val > m_maxValue)
				m_maxValue = val;
		}


		T				m_minValue,
						m_maxValue;
	};

	using Rangei  = Range<int>;
	using Rangeui = Range<unsigned>;
	using Rangef  = Range<float>;
	using Ranged  = Range<double>;



    template<typename T>
    struct ClampedValue
    {
        ClampedValue()
            : ClampedValue((T)0, T(0))
        {
        }

        ClampedValue(const ClampedValue<T>& copy )
            : ClampedValue( copy.getMin(), copy.getMax(), copy.getValue() )
        {
        }

        ClampedValue(const T& min, const T& max)
            : ClampedValue(min, max, min)
        {
        }

        ClampedValue(const T& min, const T& max, const T& curValue)
            : m_minValue(min)
            , m_maxValue(max)
            , m_curValue(curValue)
        {
        }

        const T&        getMin() const { return m_minValue; }
        const T&        getMax() const { return m_maxValue; }
        const T&        getValue() const { return m_curValue; }               

        
        void            setMin(const T& val) {
            m_minValue = val;
        }

        void            setMax(const T& val) {
            m_maxValue = val;
        }

        void            setValue(const T& val) {
            m_curValue = std::clamp(val, m_minValue, m_maxValue);
        }

        void            setFromPercentage(double percentage) {
            const auto offset = static_cast<T>( (m_maxValue - m_minValue) * percentage );
            m_curValue = m_minValue + offset;
        }

        ClampedValue<T>& operator = (const ClampedValue<T>& rhs) {
            m_minValue = rhs.m_minValue;
            m_maxValue = rhs.m_maxValue;
            m_curValue = rhs.m_curValue;

            return *this;
        }

        ClampedValue<T>& operator = (const T& val)
        {
            setValue(val);
            return *this;
        }

        bool          operator == (const ClampedValue<T>& rhs) const {
            return m_minValue == rhs.m_minValue &&
                   m_maxValue == rhs.m_maxValue && 
                   m_curValue == rhs.m_curValue;
        }

        bool          operator != (const ClampedValue<T>& rhs) const{
            return !(*this == rhs);
        }

        bool        equals( const ClampedValue<T>& rhs, const T& epsilon = (T) 0.0 ) const 
        {
            return std::abs( m_minValue - rhs.m_minValue ) < epsilon  && 
                   std::abs( m_maxValue - rhs.m_maxValue ) < epsilon  &&
                   std::abs( m_curValue - rhs.m_curValue ) < epsilon;
        }



        T               m_minValue;
        T               m_maxValue;
        T               m_curValue;
    };

    template <typename T>
    void to_json(IO::JSonObject& obj, const ClampedValue<T>& v) {
        IO::JSonObject result;
        result["Min"]   = v.getMin();
        result["Max"]   = v.getMax();
        result["Value"] = v.getValue();   
        obj = result;
    }
    template <typename T>
    void from_json(const IO::JSonObject& obj, ClampedValue<T>& v) {
        v.setMin(obj["Min"].get<T>());
        v.setMax(obj["Max"].get<T>());
        v.setValue(obj["Value"].get<T>());
    }

    template <typename T>
    void to_json(IO::JSonObject& obj, const Range<T>& v) {
        IO::JSonObject result;
        result["Min"] = v.getMin();
        result["Max"] = v.getMax();      
        obj = result;
    }
    template <typename T>
    void from_json(const IO::JSonObject& obj, Range<T>& v) {
        v.setMin(obj["Min"].get<T>());
        v.setMax(obj["Max"].get<T>());
    }

    using FloatRangedValue  = ClampedValue<float>;
    using IntRangedValue    = ClampedValue<int>;


}
