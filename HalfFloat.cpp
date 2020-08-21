#include "HalfFloat.h"

namespace Math
{

	HalfFloatLUT HalfFloat::g_LUTs;		//initialize look-up tables

	static std::uint32_t	convertMantissa(std::uint32_t index)
	{
		std::uint32_t m = index << 13;
		std::uint32_t e = 0;
		while (!(m & 0x00800000))
		{
			e -= 0x00800000;
			m <<= 1;
		}
		m &= ~0x00800000;
		e += 0x38800000;

		return (m | e);
	}

	HalfFloat::HalfFloat()
	{
		fromFloat(0.0f);
	}

	HalfFloat::HalfFloat(float val)
	{
		fromFloat(val);
	}

	HalfFloat::HalfFloat(const HalfFloat& other)
	{
		m_floatData = other.m_floatData;
	}




	HalfFloatLUT::HalfFloatLUT()
	{
		m_mantissaTable[0] = 0;
		for (int i = 1; i < 1024; ++i)
			m_mantissaTable[i] = convertMantissa(i);
		for (int i = 1024; i < 2048; ++i)
			m_mantissaTable[i] = 0x38000000 + ((i - 1024) << 13);

		for (int i = 1; i < 31; ++i)
			m_exponentTable[i] = i << 23;
		for (int i = 33; i < 63; ++i)
			m_exponentTable[i] = 0x80000000 + ((i - 32) << 23);

		m_exponentTable[0] = 0;
		m_exponentTable[31] = 0x47800000;
		m_exponentTable[32] = 0x80000000;
		m_exponentTable[63] = 0xC7800000;


		for (int i = 0; i < 64; ++i)
			m_offsetTable[i] = 1024;

		m_offsetTable[0] = 0;
		m_offsetTable[32] = 0;

		int e;

		for (int i = 0; i<256; ++i)
		{
			e = i - 127;
			if (e<-24) {
				m_baseTable[i | 0x000] = 0x0000;
				m_baseTable[i | 0x100] = 0x8000;
				m_shiftTable[i | 0x000] = 24;
				m_shiftTable[i | 0x100] = 24;
			}
			else if (e<-14) {
				m_baseTable[i | 0x000] = (0x0400 >> (-e - 14));
				m_baseTable[i | 0x100] = (0x0400 >> (-e - 14)) | 0x8000;
				m_shiftTable[i | 0x000] = -e - 1;
				m_shiftTable[i | 0x100] = -e - 1;
			}
			else if (e <= 15) {
				m_baseTable[i | 0x000] = ((e + 15) << 10);
				m_baseTable[i | 0x100] = ((e + 15) << 10) | 0x8000;
				m_shiftTable[i | 0x000] = 13;
				m_shiftTable[i | 0x100] = 13;
			}
			else if (e<128) {
				m_baseTable[i | 0x000] = 0x7C00;
				m_baseTable[i | 0x100] = 0xFC00;
				m_shiftTable[i | 0x000] = 24;
				m_shiftTable[i | 0x100] = 24;
			}
			else {
				m_baseTable[i | 0x000] = 0x7C00;
				m_baseTable[i | 0x100] = 0xFC00;
				m_shiftTable[i | 0x000] = 13;
				m_shiftTable[i | 0x100] = 13;
			}
		}
	}
}

