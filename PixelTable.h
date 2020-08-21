#pragma once
#include "MathDecl.h"
#include "MathDefs.h"
#include "GenMath.h"


namespace Math
{
	
	///////////////////////////////////////////////////////////////////////////
	//\Brief:Subdividves a rectangular( pixel ) grid in smaller grid sized portions 
	//////////////////////////////////////////////////////////////////////////
	class PixTable
	{
	public:

		PixTable()
			: m_packetSize(8)
			, m_width(0)
			, m_height(0)
			, m_roundWidth(0)
			, m_roundHeight( 0 )
			, m_numRowsX( 0 )
			, m_numRowsY( 0 )
		{

		}


		void	resize(int width, int height)
		{
			if (width == m_width && height == m_height)
				return;
			
			bool recalc = width != m_width || height != m_height;
			m_width = width;
			m_height = height;

			m_roundWidth  = Math::SnapToGrid(m_width + m_packetSize, m_packetSize);
			m_roundHeight = Math::SnapToGrid(m_height + m_packetSize, m_packetSize);

			m_numRowsX = m_roundWidth / m_packetSize;
			m_numRowsY = m_roundHeight / m_packetSize;

			rebuild();

		}
		void    rebuild()
		{
			m_pixelsX.reserve(m_roundWidth * m_roundHeight);
			m_pixelsY.reserve(m_roundWidth * m_roundHeight);
			m_pixelsX.resize(m_roundWidth * m_roundHeight);
			m_pixelsY.resize(m_roundWidth * m_roundHeight);

			int count = 0;

			for (int h = 0; h < m_roundHeight; h += m_packetSize) //loop through height
			{
				for (int w = 0; w < m_roundWidth; w += m_packetSize ) //loop through width
				{
					for (int i = 0; i < m_packetSize; ++i) //inner height loop
					{
						int startHeight = h + i;
						for (int j = 0; j < m_packetSize; ++j) //inner width loop
						{
							int startWidth = w + j;
							m_pixelsX[count] = (float)startWidth;
							m_pixelsY[count] = (float)startHeight;
							count++;
						}
					}
				}
			}
		}


		std::vector<float>			  m_pixelsX;
		std::vector<float>		  	  m_pixelsY;
		//Common::Array<Math::Vector2f> m_points;
		int							  m_packetSize;
		int							  m_width, m_height;
		int							  m_roundWidth, m_roundHeight;
		int							  m_numRowsX, m_numRowsY;



	};
}