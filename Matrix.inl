
	template< class T>
	Matrix2<T>::Matrix2()
	{
		setIdentity();
	}

	template< class T>
	Matrix2<T>::Matrix2( T _11, T _12, T _21, T _22 )
	{
		m_rows[0].setXY( _11, _12 );
		m_rows[1].setXY( _21, _22 );	
	}


	template< class T>
	Matrix2<T>::Matrix2( const T* dataPtr )
	{
		m_rows[0].setXY( dataPtr[0], dataPtr[1] );
		m_rows[1].setXY( dataPtr[2], dataPtr[3] );
	}

	template< class T>
	Matrix2<T>::Matrix2( const Vector2<T>& row1, const Vector2<T>& row2 )
	{
		m_rows[0] = row1;
		m_rows[1] = row2;
	}

	template< class T>
	Matrix2<T>::Matrix2( const Matrix2<T>& other )
	{
		memcpy( m_data, other.m_data, sizeof( T ) * 4 );
	}


	template< class T>
	void Matrix2<T>::setIdentity()
	{
		m_rows[0].setXY( (T)1, (T)0 );
		m_rows[1].setXY( (T)0, (T)1 );
	}

	template< class T>
	void Matrix2<T>::zero()
	{
		m_rows[0].setXY( (T)0, (T)0 );
		m_rows[1].setXY( (T)0, (T)0 );
	}

    template <typename T>
    Matrix3<T> Matrix2<T>::toMatrix3() const
    {
        return  Matrix3<T>( m_data[0], m_data[1], (T)(0),
                            m_data[2], m_data[3], (T)(0),
                            (T)(0)   , (T)(0)   , (T)(1));
    }

	template< class T>
	Matrix2<T> Matrix2<T>::inverted() const
	{
		Matrix2<T> inverted;
		this->invert(inverted);
		return inverted;
	}

	

	template< class T>
	bool Matrix2<T>::invert( Matrix2<T> & invertedMat ) const
	{
		T det   = determinant();
		if( std::abs(det) < (T)MATRIX_INVERT_EPSILON)
			return false;
		
		T detInv = (T)1.0 / det;
		invertedMat.m_rows[0].setXY( m_data[3] * detInv, -( m_data[1] * detInv )  );
		invertedMat.m_rows[1].setXY( -( m_data[2] * detInv ), m_data[0] * detInv );
		return true;
	}

	template<class T>
	Matrix2<T> Matrix2<T>::postMultiply( const Matrix2<T>& rhs ) const
	{
		T resultMat[4];
		//first row
		resultMat[0]    = ( m_data[0] * rhs.m_data[0] ) + ( m_data[1] * rhs.m_data[2] );
		resultMat[1]    = ( m_data[0] * rhs.m_data[1] ) + ( m_data[1] * rhs.m_data[3] );
		//second row
		resultMat[2]    = ( m_data[2] * rhs.m_data[0] ) + ( m_data[3] * rhs.m_data[2] );
		resultMat[3]    = ( m_data[2] * rhs.m_data[1] ) + ( m_data[3] * rhs.m_data[3] );
		return Matrix2<T>( resultMat );
	}

	template<class T>
	Matrix2<T> Matrix2<T>::preMultiply( const Matrix2<T>& rhs ) const
	{
		T resultMat[4];
		//first row
		resultMat[0]    = ( rhs.m_data[0] * m_data[0] ) + ( rhs.m_data[1] * m_data[2] );
		resultMat[1]    = ( rhs.m_data[0] * m_data[1] ) + ( rhs.m_data[1] * m_data[3] );
		//second row
		resultMat[2]    = ( rhs.m_data[2] * m_data[0] ) + ( rhs.m_data[3] * m_data[2] );
		resultMat[3]    = ( rhs.m_data[2] * m_data[1] ) + ( rhs.m_data[3] * m_data[3] );
		return Matrix2<T>( resultMat );
	}


	template< class T>
	T Matrix2<T>::determinant() const
	{
		return( ( m_rows[0].getX() * m_rows[1].getY() )  - ( m_rows[0].getY() * m_rows[1].getX() ) );
	}

	template< class T>
	T Matrix2<T>::trace() const
	{
		return ( m_data[0] + m_data[3] );
	}

    template <typename T>
    void Matrix2<T>::setScale(const Vector2<T>& scale)
    {
        m_data[0] = scale[0];
        m_data[3] = scale[1];
    }

	template< class T>
	void Matrix2<T>::transpose()
	{
		T tmp = m_data[2];
		m_data[2]  = m_data[1];
		m_data[1]  = tmp;
	}

	template< class T>
	void Matrix2<T>::transpose( Matrix2<T>& transposeMat ) const
	{
		transposeMat.m_rows[0].setXY( m_data[0], m_data[2] );
		transposeMat.m_rows[1].setXY( m_data[1], m_data[3] );
	}



	template< typename T>
	std::string Matrix2<T>::toString() const
	{
		char buf[256];
#pragma warning( push )
#pragma  warning (disable:4996)
		sprintf(buf, "%g %g %g %g", m_data[0], m_data[1], m_data[2], m_data[3] );
#pragma warning(pop)
		return std::string(buf);
	}

	template< typename T>
	bool Matrix2<T>::fromString(const std::string& val)
	{
		const auto tokens = Common::SplitString(val);
        if (tokens.size() != 4)
            return false;
        for (unsigned i = 0; i < tokens.size(); ++i)
            this->m_data[i] = (T)std::stod(tokens[i]);		
		return true;
	}
	
	template< class T>
	void Matrix2<T>::rotationMatrix( T angleRadians )
	{
		T cosA = cos( angleRadians );
		T sinA = sin( angleRadians );

		m_rows[0].setXY( cosA, -sinA );
		m_rows[1].setXY( sinA, cosA );
	}

	template< class T>
	T*	Matrix2<T>::toPointer()
	{
		return &m_data[0];
	}

	template< class T>
	const T*	Matrix2<T>::toConstPointer() const
	{
		return &m_data[0];
	}

	template< class T>
	Vector2<T> Matrix2<T>::multiply( const Vector2<T>& val ) const
	{
		T newX = m_rows[0].dot( val );
		T newY = m_rows[1].dot( val );
		return Vector2<T>( newX, newY );
	}
	
	template< class T>
	void Matrix2<T>::setElement( int row, int col, T val )
	{
		m_rows[row][col] = val;
	}

	template< class T>
	void Matrix2<T>::setElement( int index, T val )
	{
		m_data[index] = val;
	}

	template< typename T>
	Matrix2<T>& Matrix2<T>::operator=(const Matrix2<T>&rhs)
	{
		m_rows[0] = rhs.m_rows[0];
		m_rows[1] = rhs.m_rows[1];
		return *this;
	}
	

	template< class T>
	T Matrix2<T>::getElement( int row, int col ) const
	{
		return m_rows[row][col];
	}

	template< class T>
	T Matrix2<T>::getElement( int index ) const
	{
		return m_data[index];
	}


	template< typename T>
	void Matrix2<T>::scale( T scaleVal )
	{
		for( int i = 0; i < 4; ++i )
			m_data[i] *= scaleVal;
	}

	template< typename T>
	bool Matrix2<T>::equals( const Matrix2<T>& other, T epsilon /*=(T) 0.0001 */ ) const
	{
		for( int i = 0; i < 4; ++i )
			if( std::abs( m_data[i] - other.m_data[i]) > epsilon )
				return false;
		return true;
	}


	template< typename T>
	T& Matrix2<T>::operator[](Mat2::MATRIX2_ELEMS idx)
	{
		return m_data[idx];
	}

	template< typename T>
	T Matrix2<T>::operator[](Mat2::MATRIX2_ELEMS idx) const
	{
		return m_data[idx];
	}

	template< class T>
	const Vector2<T>& Matrix2<T>::operator[]( int index ) const
	{
		assert( index < 2 && index >= 0 );
		return m_rows[index];
	}
	
	template< class T>
	Vector2<T>& Matrix2<T>::operator[]( int index )
	{
		assert( index < 2 && index >= 0 );
		return m_rows[index];
	}
	
	template< class T>
	Vector2<T> Matrix2<T>::operator * ( const Vector2<T>& rhs ) const
	{
		return multiply( rhs );
	}

	template< class T>
	Matrix2<T>&  Matrix2<T>::operator *= ( const Matrix2<T>& rhs ) 
	{
		Matrix2<T> resultMat = postMultiply( rhs );
		*this = resultMat;
		return *this;
	}

	template< class T>
	Matrix2<T>  Matrix2<T>::operator * ( const Matrix2<T>& rhs ) const
	{
		return postMultiply( rhs );
	}

	template< typename T>
	Matrix2<T>& Matrix2<T>::operator*=( T val )
	{
		for( int i = 0; i < 4; ++i )
			m_data[i] *= val;
		return *this;
	}

	template< typename T>
	Matrix2<T> Matrix2<T>::operator-( const Matrix2<T>& rhs ) const
	{
		return Matrix2( m_rows[0] - rhs.m_rows[0], m_rows[1] - rhs.m_rows[1] );
	}

	template< typename T>
	Matrix2<T> Matrix2<T>::operator+( const Matrix2<T>& rhs ) const
	{
		return Matrix2( m_rows[0] + rhs.m_rows[0], m_rows[1] + rhs.m_rows[1] );
	}

	template< typename T>
	Matrix2<T>& Matrix2<T>::operator-=( const Matrix2<T>& rhs )
	{
		for (int i = 0; i < 2; ++i)
			m_data[i] -= rhs.m_data[i];
		return *this;
	}


	template< typename T>
	Matrix2<T>& Matrix2<T>::operator+=( const Matrix2<T>& rhs )
	{
		for( int i = 0; i < 2; ++i)
			m_data[i] += rhs.m_data[i];
		return *this;
	}

	template< typename T>
	Matrix2<T> Matrix2<T>::operator*( T val ) const
	{
		return Matrix2( m_rows[0] * val, m_rows[1] *val );
	}


	template< typename T>
	bool Matrix2<T>::operator==(const Matrix2<T>& rhs) const
	{
		return equals( rhs );
	}

	template <typename T>
	template<typename Func>	
	void Matrix2<T>::apply(Func fun)
	{
		for (auto &val : m_data)
			val = fun(val);
	}




	template< typename T>
	Matrix3<T>::Matrix3()
	{
		setIdentity();
	}

	template< typename T>
	Matrix3<T>::Matrix3( const Vector3<T>& row1, const Vector3<T>& row2, const Vector3<T>& row3 )
	{
		m_rows[0] = row1;
		m_rows[1] = row2;
		m_rows[2] = row3;
	}

	template< typename T>
	Matrix3<T>::Matrix3( T _11, T _12, T _13, T _21, T _22, T _23, T _31, T _32, T _33 )
	{
		m_rows[0].setXYZ( _11, _12, _13 );
		m_rows[1].setXYZ( _21, _22, _23 );
		m_rows[2].setXYZ( _31, _32, _33 );
	}

	template< typename T>
	Matrix3<T>::Matrix3( const T* dataPtr )
	{
		memcpy( m_data, dataPtr, sizeof(T) * 9 );
	}

	template< typename T>
	Matrix3<T>::Matrix3( const Matrix3<T>& other )
	{
		memcpy( m_data, other.m_data, sizeof(T) * 9 );
	}


	template< typename T>
	void Matrix3<T>::setIdentity()
	{
		m_rows[0].setXYZ( 1, 0, 0 );
		m_rows[1].setXYZ( 0, 1, 0 );
		m_rows[2].setXYZ( 0, 0, 1 );
	}

	template< typename T>
	void Matrix3<T>::zero()
	{
		m_rows[0].setXYZ( 0, 0, 0 );
		m_rows[1].setXYZ( 0, 0, 0 );
		m_rows[2].setXYZ( 0, 0, 0 );
	}


	template< typename T>
	std::string Matrix3<T>::toString() const
	{
		char buf[256];
#pragma warning( push )
#pragma  warning (disable:4996)
		sprintf(buf, "%g %g %g %g %g %g %g %g %g", m_data[0], m_data[1], m_data[2], m_data[3], m_data[4], m_data[5], m_data[6], m_data[7], m_data[8] );
		return std::string(buf);
#pragma  warning( pop )
	}


	template< typename T>
	bool Matrix3<T>::fromString(const std::string& val)
	{
        const auto tokens = Common::SplitString(val);
        if (tokens.size() != 9)
            return false;
        for (unsigned i = 0; i < tokens.size(); ++i)
            this->m_data[i] = (T)std::stod(tokens[i]);
        return true;
	}


	template< typename T>
	T Matrix3<T>::trace() const
	{
		return ( m_rows[0][0] + m_rows[1][1] + m_rows[2][2] );
	}

	template< typename T>
	T Matrix3<T>::determinant() const
	{
		T val = m_rows[0][0] * m_rows[1][1] * m_rows[2][2];
		val +=  m_rows[0][1] * m_rows[1][2] * m_rows[2][0];
		val +=  m_rows[0][2] * m_rows[1][0] * m_rows[2][1];

		val -=  m_rows[0][2] * m_rows[1][1] * m_rows[2][0]; 
		val -=  m_rows[0][1] * m_rows[1][0] * m_rows[2][2];
		val -=  m_rows[0][0] * m_rows[1][2] * m_rows[2][1];

		return val;
	}

	template< typename T>
	void Matrix3<T>::rotationMatrixX( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );
		//rotation over YZ plane
		m_rows[0].setXYZ( 1.0,	0.0,	0.0 );
		m_rows[1].setXYZ( 0.0,	  c,  	 s  );
		m_rows[2].setXYZ( 0.0,	  -s,	  c );
	}

	template< typename T>
	void Matrix3<T>::rotationMatrixY( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );
		//rotation XZ plane
		m_rows[0].setXYZ(   c,	0.0,	  -s );
		m_rows[1].setXYZ( 0.0,	1.0,	0.0  );
		m_rows[2].setXYZ(  s,	0.0,	  c  );
	}


	template< typename T>
	void Matrix3<T>::rotationMatrixZ( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );
		//rotation over XY plane
		m_rows[0].setXYZ(   c,	  s,	0.0	);
		m_rows[1].setXYZ(  -s,	  c,	0.0 );
		m_rows[2].setXYZ( 0.0,	0.0,	1.0 );
	}


	template< typename T>
	void  Matrix3<T>::adjoint( Matrix3<T> & adjointMat) const
	{
		adjointMat.m_data[0] = m_data[4]*m_data[8] - m_data[5]*m_data[7];
		adjointMat.m_data[1] = m_data[2]*m_data[7] - m_data[1]*m_data[8];
		adjointMat.m_data[2] = m_data[1]*m_data[5] - m_data[2]*m_data[4];

		adjointMat.m_data[3] = m_data[5]*m_data[6] - m_data[3]*m_data[8];
		adjointMat.m_data[4] = m_data[0]*m_data[8] - m_data[2]*m_data[6];
		adjointMat.m_data[5] = m_data[2]*m_data[3] - m_data[0]*m_data[5];

		adjointMat.m_data[6] = m_data[3]*m_data[7] - m_data[4]*m_data[6];
		adjointMat.m_data[7] = m_data[1]*m_data[6] - m_data[0]*m_data[7];
		adjointMat.m_data[8] = m_data[0]*m_data[4] - m_data[1]*m_data[3];
	}


	template< typename T>
	bool Matrix3<T>::equals( const Matrix3<T>& other, T epsilon /*=(T) 0.0001 */ ) const
	{
		for( int i = 0; i < 9; ++i )
			if( std::abs( m_data[i] - other.m_data[i]) > epsilon )
				return false;

		return true;
	}


	template< typename T>
	bool Matrix3<T>::invert( Matrix3<T> & invertedMat ) const
	{
		T det = determinant();
		if( det < (T)MATRIX_INVERT_EPSILON)
			return false;

		Matrix3<T>adjointMat;
		adjoint( adjointMat );

		T invDet = ((T)1 ) / det;

		for( int i = 0; i < 9; ++i )
			invertedMat.m_data[i] = adjointMat.m_data[i] * invDet;

		return true;
	}



	template< typename T>
	Matrix3<T> Matrix3<T>::inverted() const
	{
		Matrix3<T> invertMat;
		this->invert(invertMat);
		return invertMat;
	}

	template< typename T>
	void Matrix3<T>::transpose( Matrix3<T> & transMat ) const
	{

		for( int row = 0; row < 3; row++ )
			for( int col = 0; col < 3; col++ )
				transMat[col][row] = m_rows[row][col];

	}

	template< typename T>
	void Matrix3<T>::transpose()
	{
		for( int row = 0; row < 2; ++row )
		{	
			for( int col = row + 1; col < 3; ++col )
			{
				T tmp			 = m_rows[row][col];
				m_rows[row][col] = m_rows[col][row];
				m_rows[col][row] = tmp;
			}
		}
	}

	template< typename T>
	Matrix3<T> Matrix3<T>::transposed() const
	{
		Matrix3<T> tmp = *this;
		tmp.transpose();
		return tmp;
	}

	template< typename T>
	const T* Matrix3<T>::toConstPointer() const
	{
		return &m_data[0];
	}

	template< typename T>
	T* Matrix3<T>::toPointer()
	{
		return &m_data[0];
	}

	template< typename T>
	void Matrix3<T>::scale( T scaleVal )
	{
		for( int i = 0; i < 9; ++i )
			m_data[i] *= scaleVal; 
	}

    template <typename T>
    void Matrix3<T>::setScale(const Vector3<T>& val)
    {
        m_data[0] = val.getX();
        m_data[4] = val.getY();
        m_data[8] = val.getZ();
    }


	template< typename T>
	T Matrix3<T>::getElement( int index ) const
	{
		return m_data[index];
	}

	template< typename T>
	T Matrix3<T>::getElement( int row, int col ) const
	{
		return m_rows[row][col];
	}

	template< typename T>
	void Matrix3<T>::setElement( int index, T val )
	{
		m_data[index] = val;
	}

	template< typename T>
	void Matrix3<T>::setElement( int row, int col, T val )
	{
		m_rows[row][col] = val;
	}

	template< typename T>
	Matrix3<T> Matrix3<T>::preMultiply( const Matrix3<T>& rhs ) const
	{
		T resultMat[9];			

		//first row
		resultMat[0] = ( rhs.m_data[0] * m_data[0] ) + ( rhs.m_data[1] * m_data[3] ) + ( rhs.m_data[2] * m_data[6] );
		resultMat[1] = ( rhs.m_data[0] * m_data[1] ) + ( rhs.m_data[1] * m_data[4] ) + ( rhs.m_data[2] * m_data[7] );
		resultMat[2] = ( rhs.m_data[0] * m_data[2] ) + ( rhs.m_data[1] * m_data[5] ) + ( rhs.m_data[2] * m_data[8] );
		//second row																								 
		resultMat[3] = ( rhs.m_data[3] * m_data[0] ) + ( rhs.m_data[4] * m_data[3] ) + ( rhs.m_data[5] * m_data[6] );
		resultMat[4] = ( rhs.m_data[3] * m_data[1] ) + ( rhs.m_data[4] * m_data[4] ) + ( rhs.m_data[5] * m_data[7] );
		resultMat[5] = ( rhs.m_data[3] * m_data[2] ) + ( rhs.m_data[4] * m_data[5] ) + ( rhs.m_data[5] * m_data[8] );
		//third																									 
		resultMat[6] = ( rhs.m_data[6] * m_data[0] ) + ( rhs.m_data[7] * m_data[3] ) + ( rhs.m_data[8] * m_data[6] );
		resultMat[7] = ( rhs.m_data[6] * m_data[1] ) + ( rhs.m_data[7] * m_data[4] ) + ( rhs.m_data[8] * m_data[7] );
		resultMat[8] = ( rhs.m_data[6] * m_data[2] ) + ( rhs.m_data[7] * m_data[5] ) + ( rhs.m_data[8] * m_data[8] );		

		return Matrix3<T>( resultMat );
	}

	template< typename T>
	Matrix3<T> Matrix3<T>::postMultiply( const Matrix3<T>& rhs ) const
	{
		T resultMat[9];			

		//first row
		resultMat[0] = ( m_data[0] * rhs.m_data[0] ) + ( m_data[1] * rhs.m_data[3] ) + ( m_data[2] * rhs.m_data[6] );
		resultMat[1] = ( m_data[0] * rhs.m_data[1] ) + ( m_data[1] * rhs.m_data[4] ) + ( m_data[2] * rhs.m_data[7] );
		resultMat[2] = ( m_data[0] * rhs.m_data[2] ) + ( m_data[1] * rhs.m_data[5] ) + ( m_data[2] * rhs.m_data[8] );
		//second row																								 
		resultMat[3] = ( m_data[3] * rhs.m_data[0] ) + ( m_data[4] * rhs.m_data[3] ) + ( m_data[5] * rhs.m_data[6] );
		resultMat[4] = ( m_data[3] * rhs.m_data[1] ) + ( m_data[4] * rhs.m_data[4] ) + ( m_data[5] * rhs.m_data[7] );
		resultMat[5] = ( m_data[3] * rhs.m_data[2] ) + ( m_data[4] * rhs.m_data[5] ) + ( m_data[5] * rhs.m_data[8] );
		//third																									 
		resultMat[6] = ( m_data[6] * rhs.m_data[0] ) + ( m_data[7] * rhs.m_data[3] ) + ( m_data[8] * rhs.m_data[6] );
		resultMat[7] = ( m_data[6] * rhs.m_data[1] ) + ( m_data[7] * rhs.m_data[4] ) + ( m_data[8] * rhs.m_data[7] );
		resultMat[8] = ( m_data[6] * rhs.m_data[2] ) + ( m_data[7] * rhs.m_data[5] ) + ( m_data[8] * rhs.m_data[8] );		

		return Matrix3<T>( resultMat );	

	}



	template <typename T>
	template<typename Func>
	void Matrix3<T>::apply(Func fun)
	{
		for (auto &val : m_data)
			val = fun(val);
	}


	template< typename T>
	bool Matrix3<T>::operator==(const Matrix3<T>& rhs) const
	{
		return equals( rhs );
	}

	template< typename T>
	Vector3<T> Matrix3<T>::multiply( const Vector3<T>& val ) const
	{
		T newX, newY, newZ;		

		newX = m_rows[0][0] * val[0] + m_rows[1][0] * val[1] + m_rows[2][0] * val[2];//m_rows[0].dot( val ); 
		newY = m_rows[0][1] * val[0] + m_rows[1][1] * val[1] + m_rows[2][1] * val[2];//m_rows[1].dot( val ); 
		newZ = m_rows[0][2] * val[0] + m_rows[1][2] * val[1] + m_rows[2][2] * val[2];//m_rows[2].dot( val ); 

		return Vector3<T>( newX, newY, newZ );
	}

	template< typename T>
	Vector3<T> Matrix3<T>::multiplyTrans(const Vector3<T>& val) const
	{
		T newX, newY, newZ;

		newX = m_rows[0].dot(val);
		newY = m_rows[1].dot(val);
		newZ = m_rows[2].dot(val);

		//newX = m_rows[0][0] * val[0] + m_rows[1][0] * val[1] + m_rows[2][0] * val[2];//m_rows[0].dot( val ); 
		//newY = m_rows[0][1] * val[0] + m_rows[1][1] * val[1] + m_rows[2][1] * val[2];//m_rows[1].dot( val ); 
		//newZ = m_rows[0][2] * val[0] + m_rows[1][2] * val[1] + m_rows[2][2] * val[2];//m_rows[2].dot( val ); 

		return Vector3<T>(newX, newY, newZ);
	}

    template <typename T>
    Matrix2<T> Matrix3<T>::toMatrix2() const
    {
        return Matrix2<T>(m_data[0], m_data[1], m_data[3], m_data[4]);
    }

    template <typename T>
    Matrix4<T> Matrix3<T>::toMatrix4() const
    {

        return Matrix4<T>(m_data[0], m_data[1], m_data[2], T(0),
                          m_data[3], m_data[4], m_data[5], T(0),
                          m_data[6], m_data[7], m_data[8], T(0),
                          T(0), T(0), T(0), T(1)
            );
    }



	template< typename T>
	Vector3<T> Matrix3<T>::operator*( const Vector3<T>& rhs ) const
	{
		return multiply( rhs );
	}

	template< typename T>
	Matrix3<T>& Matrix3<T>::operator*=( const Matrix3<T>& rhs )
	{
		Matrix3<T> resultMat = postMultiply( rhs );
		*this = resultMat;
		return *this;
	}

	template< typename T>
	Matrix3<T> Matrix3<T>::operator*( const Matrix3<T>& rhs ) const
	{
		return postMultiply( rhs );
	}

	template< typename T>
	T& Matrix3<T>::operator[](Mat3::MATRIX3_ELEMS idx)
	{
		return m_data[idx];
	}

	template< typename T>
	T Matrix3<T>::operator[](Mat3::MATRIX3_ELEMS idx) const
	{
		return m_data[idx];
	}

	template< typename T>
	Vector3<T>& Matrix3<T>::operator[]( int rowIndex )
	{
		assert( rowIndex < 3 && rowIndex >= 0 );
		return m_rows[rowIndex];
	}

	template< typename T>
	const Vector3<T>& Matrix3<T>::operator[]( int rowIndex ) const
	{
		assert( rowIndex < 3 && rowIndex >= 0 );
		return m_rows[rowIndex];
	}


	template< typename T>
	Matrix3<T>& Matrix3<T>::operator*=( T val )
	{
		for( int i = 0; i < 9; ++i )
			m_data[i] *= val;
		return *this;
	}

	template< typename T>
	Matrix3<T>& Matrix3<T>::operator=(const Matrix3<T>&rhs)
	{
		m_rows[0] = rhs.m_rows[0];
		m_rows[1] = rhs.m_rows[1];
		m_rows[2] = rhs.m_rows[2];

		return *this;

	}


	template< typename T>
	Matrix3<T>& Matrix3<T>::operator-=( const Matrix3<T>& rhs )
	{
		for( int i = 0; i < 9; ++i )
			m_data[i] -= rhs.m_data[i];
		return *this;
	}

	template< typename T>
	Matrix3<T>& Matrix3<T>::operator+=( const Matrix3<T>& rhs )
	{
		for( int i = 0; i < 9; ++i )
			m_data[i] += rhs.m_data[i];
		return *this;
	}

	template< typename T>
	Matrix3<T> Matrix3<T>::operator+( const Matrix3<T>& rhs ) const
	{
		return Matrix3<T>( m_rows[0] + rhs.m_rows[0], m_rows[1] + rhs.m_rows[1], m_rows[2] + rhs.m_rows[2] );
	}


	template< typename T>
	Matrix3<T> Matrix3<T>::operator*( T val ) const
	{
		return Matrix3<T>( m_rows[0] * val, m_rows[1] * val, m_rows[2] * val );
	}


	template< typename T>
	Matrix3<T> Matrix3<T>::operator-( const Matrix3<T>& rhs ) const
	{
		return Matrix3<T>( m_rows[0] - rhs.m_rows[0], m_rows[1] - rhs.m_rows[1], m_rows[2] - rhs.m_rows[2] );
	}



	template< typename T>
	Matrix4<T>::Matrix4()
	{
		setIdentity();
	}

	template< typename T>
	Matrix4<T>::Matrix4( const Vector4<T>& row1, const Vector4<T>& row2, const Vector4<T>& row3, const Vector4<T>& row4 )
	{
		m_rows[0] = row1;
		m_rows[1] = row2;
		m_rows[2] = row3;
		m_rows[3] = row4;
	}

	template< typename T>
	Matrix4<T>::Matrix4( T _11, T _12, T _13, T _14, T _21, T _22, T _23, T _24, T _31, T _32, T _33, T _34, T _41, T _42, T _43, T _44 )
	{
		m_rows[0].setXYZW( _11, _12, _13, _14 );
		m_rows[1].setXYZW( _21, _22, _23, _24 );
		m_rows[2].setXYZW( _31, _32, _33, _34 );
		m_rows[3].setXYZW( _41, _42, _43, _44 );
	}

	template< typename T>
	Matrix4<T>::Matrix4( const T* dataPtr )
	{
		memcpy( m_data, dataPtr, sizeof( T) * 16 );
	}

	template< typename T>
	Matrix4<T>::Matrix4( const Matrix4<T>& other )
	{
        /*m_rows[0] = other.m_rows[0];
        m_rows[1] = other.m_rows[1];
        m_rows[2] = other.m_rows[2];
        m_rows[3] = other.m_rows[3];*/
        memcpy( m_data, other.m_data, sizeof(T) * 16 );
	}


    template< typename T>
    Matrix4<T> Matrix4<T>::transposed() const
    {
        Matrix4<T> result;
        transpose( result );
        return result;
    }


	


	template< typename T>
	std::string Matrix4<T>::toString() const
	{
		char buf[256];
#pragma warning( push )
#pragma  warning (disable:4996)
		sprintf(buf, "%g %g %g %g \
					 %g %g %g %g \
					 %g %g %g %g \
					 %g %g %g %g", 
					 m_data[0], m_data[1], m_data[2], m_data[3], 
					 m_data[4], m_data[5], m_data[6], m_data[7], 
					 m_data[8], m_data[9], m_data[10], m_data[11],
					 m_data[12], m_data[13], m_data[14], m_data[15] );
#pragma  warning ( pop )
		return std::string(buf);
	}


	template< typename T>
	bool Matrix4<T>::fromString(const std::string& val)
	{
        const auto tokens = Common::SplitString(val);
        if (tokens.size() != 16)
            return false;
        for (unsigned i = 0; i < tokens.size(); ++i)
            this->m_data[i] = (T)std::stod(tokens[i]);
        return true;
	}

	template< typename T>
	void Matrix4<T>::setIdentity()
	{
		m_rows[0].setXYZW( 1, 0, 0, 0 );
		m_rows[1].setXYZW( 0, 1, 0, 0 );
		m_rows[2].setXYZW( 0, 0, 1, 0 );
		m_rows[3].setXYZW( 0, 0, 0, 1 );
	}

	template< typename T>
	void Matrix4<T>::zero()
	{
		m_rows[0].setXYZW( 0, 0, 0, 0 );
		m_rows[1].setXYZW( 0, 0, 0, 0 );
		m_rows[2].setXYZW( 0, 0, 0, 0 );
		m_rows[3].setXYZW( 0, 0, 0, 0 );
	}


    template< typename T>
    void Matrix4<T>::setTranslation( const Vector3<T>& translation )
    {
		m_rows[3][0] = translation[0];
		m_rows[3][1] = translation[1];
		m_rows[3][2] = translation[2];
    }

    template< typename T>
    Vector3<T>Matrix4<T>::getTranslation() const
    {
		return m_rows[3].toVector3();
    }

    template< typename T>
    void Matrix4<T>::setScale( const Vector3<T>& scale )
    {
        m_rows[0][0] = scale[0];
        m_rows[1][1] = scale[1];
        m_rows[2][2] = scale[2];
    }

    template< typename T>
    void Matrix4<T>::setDirection( const Vector3<T>& direction )
    {
        m_rows[2][0] = direction[0];
        m_rows[2][1] = direction[1];
        m_rows[2][2] = direction[2];
    }

    template< typename T>
    void Matrix4<T>::setUp( const Vector3<T>& up )
    {
        m_rows[1][0] = up[0];
        m_rows[1][1] = up[1];
        m_rows[1][2] = up[2];
    }

    template< typename T>
    void Matrix4<T>::setRight( const Vector3<T>& right )
    {
        m_rows[0][0] = right[0];
        m_rows[0][1] = right[1];
        m_rows[0][2] = right[2];
    }

    template< typename T>
    Vector3<T> Matrix4<T>::getRight() const
    {
        return Vector3<T>( m_rows[0][0], m_rows[0][1], m_rows[0][2] );
    }

    template< typename T>
    Vector3<T> Matrix4<T>::getUp() const
    {
        return Vector3<T>( m_rows[1][0], m_rows[1][1], m_rows[1][2] );
    }

    template< typename T>
    Vector3<T> Matrix4<T>::getDirection() const
    {
        return Vector3<T>( m_rows[2][0], m_rows[2][1], m_rows[2][2] );
    }

	template< typename T>
	Vector3<T> Matrix4<T>::getScale(bool normalize 	) const
	{
		if( normalize )
		{
			return Vector3<T>
				(  (T)sqrt( m_rows[0][0] * m_rows[0][0] + m_rows[0][1] * m_rows[0][1]  + m_rows[0][2] * m_rows[0][2] ),
				   (T)sqrt( m_rows[1][0] * m_rows[1][0] + m_rows[1][1] * m_rows[1][1]  + m_rows[1][2] * m_rows[1][2] ),
				   (T)sqrt( m_rows[2][0] * m_rows[2][0] + m_rows[2][1] * m_rows[2][1]  + m_rows[2][2] * m_rows[2][2] ) );
		}
		
		return Vector3<T>( m_rows[0][0], m_rows[1][1], m_rows[2][2] );
	}


    template <typename T>
    Matrix4<T> Matrix4<T>::CreateScaleMatrix(const Vector3<T>& xyz)
    {
        Matrix4<T> result;
        result.setScale(xyz);
        return result;
    }
	
    template <typename T>
    Matrix4<T> Matrix4<T>::CreateTranslationMatrix(const Vector3<T>& xyz)
    {
        Matrix4<T> result;
        result.setTranslation(xyz);
        return result;
    }



	template< typename T>
	void Matrix4<T>::rotationMatrixX( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );

		//rotation over YZ plane
		m_rows[0].setXYZW( 1.0,		0.0,	0.0,	0.0 );
		m_rows[1].setXYZW( 0.0,		  c,	  s,	0.0 );
		m_rows[2].setXYZW( 0.0,		  -s,	  c,	0.0 );
		m_rows[3].setXYZW( 0.0,		0.0,	0.0 ,	1.0 );

	}

	template< typename T>
	void Matrix4<T>::rotationMatrixY( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );

		//rotation over YZ plane
		m_rows[0].setXYZW(   c,		0.0,	  -s,	0.0 );
		m_rows[1].setXYZW( 0.0,		1.0,	0.0,	0.0 );
		m_rows[2].setXYZW(  s,		0.0,	  c,	0.0 );
		m_rows[3].setXYZW( 0.0,		0.0,	0.0 ,	1.0 );
	}

	template< typename T>
	void Matrix4<T>::rotationMatrixZ( T angleRadians )
	{
		T s, c;
		sinCos( angleRadians, s, c );

		//rotation over XY plane
		m_rows[0].setXYZW(   c,		  s,	0.0,	0.0 );
		m_rows[1].setXYZW(  -s,		  c,	0.0,	0.0 );
		m_rows[2].setXYZW( 0.0,		0.0,	1.0,	0.0 );
		m_rows[3].setXYZW( 0.0,		0.0,	0.0 ,	1.0 );
	}


	template< typename T>
	Matrix4<T> Matrix4<T>::inverted() const
	{
		Matrix4<T> ret;
		invert( ret );
		return ret;
	}



	template <class T>
	template <class U>
	Matrix4<U>	Matrix4<T>::convertTo() const
	{

		Matrix4<U> newMat;
		for( int i = 0; i < 16; ++i )
			newMat.setElement( i, (U) m_data[i] );
		return newMat;
	}


	template< typename T>
	T Matrix4<T>::trace() const
	{
		return ( m_data[0] + m_data[5] + m_data[10] + m_data[15] );
	}

	template< typename T>
	void Matrix4<T>::transpose()
	{
		for( int row = 0; row < 3; row++ )
		{	
			for( int col = row + 1; col < 4; col++ )
			{
				T tmp			 = m_rows[row][col];
				m_rows[row][col] = m_rows[col][row];
				m_rows[col][row] = tmp;
			}
		}
	}



	template< typename T>
	void Matrix4<T>::transpose( Matrix4<T>& transposeOut ) const
	{
		for( int row = 0; row < 4; row++ )
			for( int col = 0; col < 4; col++ )
				transposeOut[col][row] = m_rows[row][col];
	}

	template< typename T>
	Matrix4<T>& Matrix4<T>::operator=(const Matrix4<T>&rhs)
	{
	    memcpy( m_data, rhs.m_data, sizeof(T) * 16 );
		return *this;
	}




    template< typename T>
   void Matrix4<T>::lookAt( const Vector3<T>& eye, const Vector3<T>& target, const Vector3<T>& up )
   {
        
        
        auto dir = (target - eye).getNormalized();
        auto side = dir.cross( up ).getNormalized();
        auto newUp = side.cross( dir ).getNormalized();

        setDirection( dir *-1.0f );
        setRight( side );
        setUp( newUp );

        auto eyeOrig = eye * -1.0f;

        setTranslation( Math::Vector3f( eyeOrig.dot( side ), eyeOrig.dot( newUp ), eyeOrig.dot( dir * -1.0f ) ) );



    }


   template<typename T>
   template<typename Func>
   void Matrix4<T>::apply(Func fun)
   {
	   for (auto &val : m_data)
		   val = fun(val);
   }




	template< typename T>
	bool Matrix4<T>::invert( Matrix4<T>& inverseOut ) const
	{
        T inv[16];
        T det;
        int i;

        inv[0] = m_data[5] * m_data[10] * m_data[15] -
            m_data[5] * m_data[11] * m_data[14] -
            m_data[9] * m_data[6] * m_data[15] +
            m_data[9] * m_data[7] * m_data[14] +
            m_data[13] * m_data[6] * m_data[11] -
            m_data[13] * m_data[7] * m_data[10];

        inv[4] = -m_data[4] * m_data[10] * m_data[15] +
            m_data[4] * m_data[11] * m_data[14] +
            m_data[8] * m_data[6] * m_data[15] -
            m_data[8] * m_data[7] * m_data[14] -
            m_data[12] * m_data[6] * m_data[11] +
            m_data[12] * m_data[7] * m_data[10];

        inv[8] = m_data[4] * m_data[9] * m_data[15] -
            m_data[4] * m_data[11] * m_data[13] -
            m_data[8] * m_data[5] * m_data[15] +
            m_data[8] * m_data[7] * m_data[13] +
            m_data[12] * m_data[5] * m_data[11] -
            m_data[12] * m_data[7] * m_data[9];

        inv[12] = -m_data[4] * m_data[9] * m_data[14] +
            m_data[4] * m_data[10] * m_data[13] +
            m_data[8] * m_data[5] * m_data[14] -
            m_data[8] * m_data[6] * m_data[13] -
            m_data[12] * m_data[5] * m_data[10] +
            m_data[12] * m_data[6] * m_data[9];

        inv[1] = -m_data[1] * m_data[10] * m_data[15] +
            m_data[1] * m_data[11] * m_data[14] +
            m_data[9] * m_data[2] * m_data[15] -
            m_data[9] * m_data[3] * m_data[14] -
            m_data[13] * m_data[2] * m_data[11] +
            m_data[13] * m_data[3] * m_data[10];

        inv[5] = m_data[0] * m_data[10] * m_data[15] -
            m_data[0] * m_data[11] * m_data[14] -
            m_data[8] * m_data[2] * m_data[15] +
            m_data[8] * m_data[3] * m_data[14] +
            m_data[12] * m_data[2] * m_data[11] -
            m_data[12] * m_data[3] * m_data[10];

        inv[9] = -m_data[0] * m_data[9] * m_data[15] +
            m_data[0] * m_data[11] * m_data[13] +
            m_data[8] * m_data[1] * m_data[15] -
            m_data[8] * m_data[3] * m_data[13] -
            m_data[12] * m_data[1] * m_data[11] +
            m_data[12] * m_data[3] * m_data[9];

        inv[13] = m_data[0] * m_data[9] * m_data[14] -
            m_data[0] * m_data[10] * m_data[13] -
            m_data[8] * m_data[1] * m_data[14] +
            m_data[8] * m_data[2] * m_data[13] +
            m_data[12] * m_data[1] * m_data[10] -
            m_data[12] * m_data[2] * m_data[9];

        inv[2] = m_data[1] * m_data[6] * m_data[15] -
            m_data[1] * m_data[7] * m_data[14] -
            m_data[5] * m_data[2] * m_data[15] +
            m_data[5] * m_data[3] * m_data[14] +
            m_data[13] * m_data[2] * m_data[7] -
            m_data[13] * m_data[3] * m_data[6];

        inv[6] = -m_data[0] * m_data[6] * m_data[15] +
            m_data[0] * m_data[7] * m_data[14] +
            m_data[4] * m_data[2] * m_data[15] -
            m_data[4] * m_data[3] * m_data[14] -
            m_data[12] * m_data[2] * m_data[7] +
            m_data[12] * m_data[3] * m_data[6];

        inv[10] = m_data[0] * m_data[5] * m_data[15] -
            m_data[0] * m_data[7] * m_data[13] -
            m_data[4] * m_data[1] * m_data[15] +
            m_data[4] * m_data[3] * m_data[13] +
            m_data[12] * m_data[1] * m_data[7] -
            m_data[12] * m_data[3] * m_data[5];

        inv[14] = -m_data[0] * m_data[5] * m_data[14] +
            m_data[0] * m_data[6] * m_data[13] +
            m_data[4] * m_data[1] * m_data[14] -
            m_data[4] * m_data[2] * m_data[13] -
            m_data[12] * m_data[1] * m_data[6] +
            m_data[12] * m_data[2] * m_data[5];

        inv[3] = -m_data[1] * m_data[6] * m_data[11] +
            m_data[1] * m_data[7] * m_data[10] +
            m_data[5] * m_data[2] * m_data[11] -
            m_data[5] * m_data[3] * m_data[10] -
            m_data[9] * m_data[2] * m_data[7] +
            m_data[9] * m_data[3] * m_data[6];

        inv[7] = m_data[0] * m_data[6] * m_data[11] -
            m_data[0] * m_data[7] * m_data[10] -
            m_data[4] * m_data[2] * m_data[11] +
            m_data[4] * m_data[3] * m_data[10] +
            m_data[8] * m_data[2] * m_data[7] -
            m_data[8] * m_data[3] * m_data[6];

        inv[11] = -m_data[0] * m_data[5] * m_data[11] +
            m_data[0] * m_data[7] * m_data[9] +
            m_data[4] * m_data[1] * m_data[11] -
            m_data[4] * m_data[3] * m_data[9] -
            m_data[8] * m_data[1] * m_data[7] +
            m_data[8] * m_data[3] * m_data[5];

        inv[15] = m_data[0] * m_data[5] * m_data[10] -
            m_data[0] * m_data[6] * m_data[9] -
            m_data[4] * m_data[1] * m_data[10] +
            m_data[4] * m_data[2] * m_data[9] +
            m_data[8] * m_data[1] * m_data[6] -
            m_data[8] * m_data[2] * m_data[5];

        det = m_data[0] * inv[0] + m_data[1] * inv[4] + m_data[2] * inv[8] + m_data[3] * inv[12];

		if (std::abs(det) < (T)MATRIX_INVERT_EPSILON)
            return false;

        det = (T)1.0 / det;

        for(i = 0; i < 16; i++)
            inverseOut.m_data[i] = inv[i] * det;

        return true;
        
        
        
        /*T a0 = m_data[ 0]*m_data[ 5] - m_data[ 1]*m_data[ 4];
		T a1 = m_data[ 0]*m_data[ 6] - m_data[ 2]*m_data[ 4];
		T a2 = m_data[ 0]*m_data[ 7] - m_data[ 3]*m_data[ 4];
		T a3 = m_data[ 1]*m_data[ 6] - m_data[ 2]*m_data[ 5];
		T a4 = m_data[ 1]*m_data[ 7] - m_data[ 3]*m_data[ 5];
		T a5 = m_data[ 2]*m_data[ 7] - m_data[ 3]*m_data[ 6];
		T b0 = m_data[ 8]*m_data[13] - m_data[ 9]*m_data[12];
		T b1 = m_data[ 8]*m_data[14] - m_data[10]*m_data[12];
		T b2 = m_data[ 8]*m_data[15] - m_data[11]*m_data[12];
		T b3 = m_data[ 9]*m_data[14] - m_data[10]*m_data[13];
		T b4 = m_data[ 9]*m_data[15] - m_data[11]*m_data[13];
		T b5 = m_data[10]*m_data[15] - m_data[11]*m_data[14];

		T det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;

		if( abs( det) < (T) 0.00000000001 )
			return false;

		inverseOut.m_data[0 ] =	+ m_data[ 5]*b5 - m_data[ 6]*b4 + m_data[ 7]*b3;
		inverseOut.m_data[1 ] =	- m_data[ 1]*b5 + m_data[ 2]*b4 - m_data[ 3]*b3;
		inverseOut.m_data[2 ] =	+ m_data[13]*a5 - m_data[14]*a4 + m_data[15]*a3;
		inverseOut.m_data[3 ] =	- m_data[ 9]*a5 + m_data[10]*a4 - m_data[11]*a3;

		inverseOut.m_data[4 ] =	- m_data[ 4]*b5 + m_data[ 6]*b2 - m_data[ 7]*b1;
		inverseOut.m_data[5 ] =	+ m_data[ 0]*b5 - m_data[ 2]*b2 + m_data[ 3]*b1;
		inverseOut.m_data[6 ] =	- m_data[12]*a5 + m_data[14]*a2 - m_data[15]*a1;
		inverseOut.m_data[7 ] =	+ m_data[ 8]*a5 - m_data[10]*a2 + m_data[11]*a1;

		inverseOut.m_data[8 ] = + m_data[ 4]*b4 - m_data[ 5]*b2 + m_data[ 7]*b0;
		inverseOut.m_data[9 ] = - m_data[ 0]*b4 + m_data[ 1]*b2 - m_data[ 3]*b0;
		inverseOut.m_data[10] =	+ m_data[12]*a4 - m_data[13]*a2 + m_data[15]*a0;
		inverseOut.m_data[11] =	- m_data[ 8]*a4 + m_data[ 9]*a2 - m_data[11]*a0;

		inverseOut.m_data[12] = - m_data[ 4]*b3 + m_data[ 5]*b1 - m_data[ 6]*b0;
		inverseOut.m_data[13] = + m_data[ 0]*b3 - m_data[ 1]*b1 + m_data[ 2]*b0;
		inverseOut.m_data[14] = - m_data[12]*a3 + m_data[13]*a1 - m_data[14]*a0;
		inverseOut.m_data[15] = + m_data[ 8]*a3 - m_data[ 9]*a1 + m_data[10]*a0;

		T invDet = (T)(1.0/det);

		inverseOut *= invDet;

		return true;*/
	}

    template <typename T>
    void Matrix4<T>::setScale(const Vector4<T>& val)
    {
        m_data[0] = val.getX();
        m_data[5] = val.getY();
        m_data[10] = val.getZ();
        m_data[15] = val.getW();
    }



	template< typename T>
	void Matrix4<T>::perspectiveMatrix( T angle, T viewportRatio, T nearPlane, T farPlane )
	{
		setIdentity();				

		T xmin, xmax, ymin, ymax;
		T one_deltax, one_deltay, one_deltaz, doubleznear;

		ymax = nearPlane * tan( angle );
		ymin = -ymax;

		xmax = ymax * viewportRatio;
		xmin = -xmax;

		doubleznear = (T)2.0 * nearPlane;
		one_deltax	= (T)1.0 / ( xmax - xmin );
		one_deltay	= (T)1.0 / ( ymax - ymin );
		one_deltaz	= (T)1.0 / ( farPlane - nearPlane );

		setElement(  0,  ( doubleznear * one_deltax ) );
		setElement(  2,  ( ( xmax + xmin ) * one_deltax ) );
		setElement(  5,  ( doubleznear * one_deltay ) );
		setElement(  6, ( ( ymax + ymin ) * one_deltay ) );
		setElement( 10, (-( farPlane+ nearPlane ) * one_deltaz ) );
		setElement( 11, (- ( farPlane * doubleznear ) * one_deltaz ) );
		setElement( 14, (T)-1.0 );
		setElement( 15, (T) 0.0 );
		//transpose();
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::Perspective(T fovy, T aspect, T zNear, T zFar)
	{
		const T bottom = -zNear * (T)tanf(0.5f * fovy * 3.14159265358 / 180.0f);
		const T top = -bottom;

		const T left = aspect * bottom;
		const T right = -left;

		return Frustum(left, right, bottom, top, zNear, zFar);
	}


	template< typename T>
	Matrix4<T> Matrix4<T>::Ortho2DMatrix(T width, T height)
	{
		Matrix4<T> result;
		result.orthographicMatrix2d(width, height);
		return result;
	}


	template< typename T>
	Matrix4<T> Matrix4<T>::Frustum(T left, T right, T bottom, T top, T zNear, T zFar)
	{
		const T dx = right - left;
		const T dy = bottom - top;
		const T dz = zFar - zNear;

		const T mx = (T)0.5 * (left + right);
		const T my = (T)0.5 * (bottom + top);
		const T mz = (T)0.5 * (zNear + zFar);

		const T n = zNear;
		const T nf = zNear * zFar;

		return Matrix4<T>(2.0f * n / dx, 0.0f, 2.0f * mx / dx, 0.0f,
						  0.0f, 2.0f * n / dy, 2.0f * my / dy, 0.0f,
						  0.0f, 0.0f, -2.0f * mz / dz, -2.0f * nf / dz,
						  0.0f, 0.0f, -1.0f, 0.0f);
	}

	template< typename T>
	void Matrix4<T>::orthographicMatrix( T left, T right, T bottom, T top, T znear, T zfar )
	{
		
		const T invX = (T)1.0 / ( right - left );
		const T invY = (T)1.0 / ( top - bottom );
		const T invZ = (T)1.0 / ( zfar - znear );
		
		setElement(  0, (T) 2.0 * invX );
		setElement(  5, (T) 2.0 * invY );
		setElement( 10, (T)-2.0 * invZ );

		setElement( 3, ( -( right + left ) * invX ) );
		setElement( 7, ( -( top + bottom ) * invY ) );
		setElement( 11,( -( zfar + znear ) * invZ ) );
		setElement( 15, (T) 1.0 );
	}



	template< typename T>
	Matrix3<T> Matrix4<T>::toMatrix3() const
	{
		Vector3<T> rows[3];
		for( int i = 0; i < 3; i++ )
			rows[i].setXYZ( m_rows[i][0], m_rows[i][1], m_rows[i][2] );

		return Matrix3<T>( rows[0], rows[1], rows[2] );
	}




	template< typename T>
	T Matrix4<T>::determinant() const
	{
		T a0 = m_data[ 0]*m_data[ 5] - m_data[ 1]*m_data[ 4];
		T a1 = m_data[ 0]*m_data[ 6] - m_data[ 2]*m_data[ 4];
		T a2 = m_data[ 0]*m_data[ 7] - m_data[ 3]*m_data[ 4];
		T a3 = m_data[ 1]*m_data[ 6] - m_data[ 2]*m_data[ 5];
		T a4 = m_data[ 1]*m_data[ 7] - m_data[ 3]*m_data[ 5];
		T a5 = m_data[ 2]*m_data[ 7] - m_data[ 3]*m_data[ 6];
		T b0 = m_data[ 8]*m_data[13] - m_data[ 9]*m_data[12];
		T b1 = m_data[ 8]*m_data[14] - m_data[10]*m_data[12];
		T b2 = m_data[ 8]*m_data[15] - m_data[11]*m_data[12];
		T b3 = m_data[ 9]*m_data[14] - m_data[10]*m_data[13];
		T b4 = m_data[ 9]*m_data[15] - m_data[11]*m_data[13];
		T b5 = m_data[10]*m_data[15] - m_data[11]*m_data[14];

		T det = a0*b5 - a1*b4 + a2*b3 + a3*b2 - a4*b1 + a5*b0;
		return det;
	}


	template< typename T>
	void  Matrix4<T>::adjoint( Matrix4<T>& adjointMat ) 
	{
		T a0 = m_data[ 0]*m_data[ 5] - m_data[ 1]*m_data[ 4];
		T a1 = m_data[ 0]*m_data[ 6] - m_data[ 2]*m_data[ 4];
		T a2 = m_data[ 0]*m_data[ 7] - m_data[ 3]*m_data[ 4];
		T a3 = m_data[ 1]*m_data[ 6] - m_data[ 2]*m_data[ 5];
		T a4 = m_data[ 1]*m_data[ 7] - m_data[ 3]*m_data[ 5];
		T a5 = m_data[ 2]*m_data[ 7] - m_data[ 3]*m_data[ 6];
		T b0 = m_data[ 8]*m_data[13] - m_data[ 9]*m_data[12];
		T b1 = m_data[ 8]*m_data[14] - m_data[10]*m_data[12];
		T b2 = m_data[ 8]*m_data[15] - m_data[11]*m_data[12];
		T b3 = m_data[ 9]*m_data[14] - m_data[10]*m_data[13];
		T b4 = m_data[ 9]*m_data[15] - m_data[11]*m_data[13];
		T b5 = m_data[10]*m_data[15] - m_data[11]*m_data[14];

		adjointMat.m_data[0 ] =	+ m_data[ 5]*b5 - m_data[ 6]*b4 + m_data[ 7]*b3;
		adjointMat.m_data[1 ] =	- m_data[ 1]*b5 + m_data[ 2]*b4 - m_data[ 3]*b3;
		adjointMat.m_data[2 ] =	+ m_data[13]*a5 - m_data[14]*a4 + m_data[15]*a3;
		adjointMat.m_data[3 ] =	- m_data[ 9]*a5 + m_data[10]*a4 - m_data[11]*a3;

		adjointMat.m_data[4 ] =	- m_data[ 4]*b5 + m_data[ 6]*b2 - m_data[ 7]*b1;
		adjointMat.m_data[5 ] =	+ m_data[ 0]*b5 - m_data[ 2]*b2 + m_data[ 3]*b1;
		adjointMat.m_data[6 ] =	- m_data[12]*a5 + m_data[14]*a2 - m_data[15]*a1;
		adjointMat.m_data[7 ] =	+ m_data[ 8]*a5 - m_data[10]*a2 + m_data[11]*a1;

		adjointMat.m_data[8 ] = + m_data[ 4]*b4 - m_data[ 5]*b2 + m_data[ 7]*b0;
		adjointMat.m_data[9 ] = - m_data[ 0]*b4 + m_data[ 1]*b2 - m_data[ 3]*b0;
		adjointMat.m_data[10] =	+ m_data[12]*a4 - m_data[13]*a2 + m_data[15]*a0;
		adjointMat.m_data[11] =	- m_data[ 8]*a4 + m_data[ 9]*a2 - m_data[11]*a0;

		adjointMat.m_data[12] = - m_data[ 4]*b3 + m_data[ 5]*b1 - m_data[ 6]*b0;
		adjointMat.m_data[13] = + m_data[ 0]*b3 - m_data[ 1]*b1 + m_data[ 2]*b0;
		adjointMat.m_data[14] = - m_data[12]*a3 + m_data[13]*a1 - m_data[14]*a0;
		adjointMat.m_data[15] = + m_data[ 8]*a3 - m_data[ 9]*a1 + m_data[10]*a0;

	}


	template< typename T>
	void Matrix4<T>::scale( T scaleVal )
	{
		for( int i = 0; i < 16; ++i )
			m_data[i] *= scaleVal;
	}

	template< typename T>
	bool Matrix4<T>::equals( const Matrix4<T>& other, T epsilon /*=(T) 0.0001 */ ) const
	{
		for( int i= 0; i < 16; ++i )
			if( std::abs( m_data[i] - other.m_data[i]) > epsilon )
				return false;
		return true;
	}

		
	template< typename T>
	Vector3<T> Matrix4<T>::multiply( const Vector3<T>& val ) const
	{
		Vector3<T> result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = m_data[i +  0]	* val[0] + 
						m_data[i +  4]	* val[1] + 
						m_data[i +  8]	* val[2] + 
						m_data[i + 12];
		}	
		return  result;
	}

	template< typename T>
	Vector3<T> Matrix4<T>::transformNormal(const Vector3<T>& val) const
	{
		Vector3<T> result;
		for (int i = 0; i < 3; ++i)
		{
			result[i] = m_data[i + 0] * val[0] +
						m_data[i + 4] * val[1] +
						m_data[i + 8] * val[2];
		}
		return  result;
	}


	template< typename T>
	Vector4<T> Matrix4<T>::multiply( const Vector4<T>& val ) const
	{
		Vector4<T> result;
		for (int i = 0; i < 4; ++i)
		{
			result[i] = m_data[i +  0] * val[0] +
						m_data[i +  4] * val[1] +
						m_data[i +  8] * val[2] +
						m_data[i + 12] * val[3];
		}
		return result;
	}


	template< typename T>
	Matrix4<T> Matrix4<T>::postMultiply( const Matrix4<T>& rhs ) const
	{
		T resultMat[16];

		resultMat[ 0] = ( m_data[ 0] * rhs.m_data[0] ) + ( m_data[ 1] * rhs.m_data[4] ) + ( m_data[ 2] * rhs.m_data[ 8] ) + ( m_data[ 3] * rhs.m_data[12] );
		resultMat[ 1] = ( m_data[ 0] * rhs.m_data[1] ) + ( m_data[ 1] * rhs.m_data[5] ) + ( m_data[ 2] * rhs.m_data[ 9] ) + ( m_data[ 3] * rhs.m_data[13] );
		resultMat[ 2] = ( m_data[ 0] * rhs.m_data[2] ) + ( m_data[ 1] * rhs.m_data[6] ) + ( m_data[ 2] * rhs.m_data[10] ) + ( m_data[ 3] * rhs.m_data[14] );
		resultMat[ 3] = ( m_data[ 0] * rhs.m_data[3] ) + ( m_data[ 1] * rhs.m_data[7] ) + ( m_data[ 2] * rhs.m_data[11] ) + ( m_data[ 3] * rhs.m_data[15] );

		resultMat[ 4] = ( m_data[ 4] * rhs.m_data[0] ) + ( m_data[ 5] * rhs.m_data[4] ) + ( m_data[ 6] * rhs.m_data[ 8] ) + ( m_data[ 7] * rhs.m_data[12] );
		resultMat[ 5] = ( m_data[ 4] * rhs.m_data[1] ) + ( m_data[ 5] * rhs.m_data[5] ) + ( m_data[ 6] * rhs.m_data[ 9] ) + ( m_data[ 7] * rhs.m_data[13] );
		resultMat[ 6] = ( m_data[ 4] * rhs.m_data[2] ) + ( m_data[ 5] * rhs.m_data[6] ) + ( m_data[ 6] * rhs.m_data[10] ) + ( m_data[ 7] * rhs.m_data[14] );
		resultMat[ 7] = ( m_data[ 4] * rhs.m_data[3] ) + ( m_data[ 5] * rhs.m_data[7] ) + ( m_data[ 6] * rhs.m_data[11] ) + ( m_data[ 7] * rhs.m_data[15] );

		resultMat[ 8] = ( m_data[ 8] * rhs.m_data[0] ) + ( m_data[ 9] * rhs.m_data[4] ) + ( m_data[10] * rhs.m_data[8 ] ) + ( m_data[11] * rhs.m_data[12] );
		resultMat[ 9] = ( m_data[ 8] * rhs.m_data[1] ) + ( m_data[ 9] * rhs.m_data[5] ) + ( m_data[10] * rhs.m_data[9 ] ) + ( m_data[11] * rhs.m_data[13] );
		resultMat[10] = ( m_data[ 8] * rhs.m_data[2] ) + ( m_data[ 9] * rhs.m_data[6] ) + ( m_data[10] * rhs.m_data[10] ) + ( m_data[11] * rhs.m_data[14] );
		resultMat[11] = ( m_data[ 8] * rhs.m_data[3] ) + ( m_data[ 9] * rhs.m_data[7] ) + ( m_data[10] * rhs.m_data[11] ) + ( m_data[11] * rhs.m_data[15] );

		resultMat[12] = ( m_data[12] * rhs.m_data[0] ) + ( m_data[13] * rhs.m_data[4] ) + ( m_data[14] * rhs.m_data[8 ] ) + ( m_data[15] * rhs.m_data[12] );
		resultMat[13] = ( m_data[12] * rhs.m_data[1] ) + ( m_data[13] * rhs.m_data[5] ) + ( m_data[14] * rhs.m_data[9 ] ) + ( m_data[15] * rhs.m_data[13] );
		resultMat[14] = ( m_data[12] * rhs.m_data[2] ) + ( m_data[13] * rhs.m_data[6] ) + ( m_data[14] * rhs.m_data[10] ) + ( m_data[15] * rhs.m_data[14] );
		resultMat[15] = ( m_data[12] * rhs.m_data[3] ) + ( m_data[13] * rhs.m_data[7] ) + ( m_data[14] * rhs.m_data[11] ) + ( m_data[15] * rhs.m_data[15] );

		return Matrix4<T>( resultMat );
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::preMultiply( const Matrix4<T>& rhs ) const
	{
		T resultMat[16];

		resultMat[ 0] = ( rhs.m_data[ 0] * m_data[0] ) + ( rhs.m_data[ 1] * m_data[4] ) + ( rhs.m_data[ 2] * m_data[ 8] ) + (  rhs.m_data[ 3] *m_data[12] );
		resultMat[ 1] = ( rhs.m_data[ 0] * m_data[1] ) + ( rhs.m_data[ 1] * m_data[5] ) + ( rhs.m_data[ 2] * m_data[ 9] ) + (  rhs.m_data[ 3] *m_data[13] );
		resultMat[ 2] = ( rhs.m_data[ 0] * m_data[2] ) + ( rhs.m_data[ 1] * m_data[6] ) + ( rhs.m_data[ 2] * m_data[10] ) + (  rhs.m_data[ 3] *m_data[14] );
		resultMat[ 3] = ( rhs.m_data[ 0] * m_data[3] ) + ( rhs.m_data[ 1] * m_data[7] ) + ( rhs.m_data[ 2] * m_data[11] ) + (  rhs.m_data[ 3] *m_data[15] );

		resultMat[ 4] = ( rhs.m_data[ 4] * m_data[0] ) + ( rhs.m_data[ 5] * m_data[4] ) + ( rhs.m_data[ 6] * m_data[ 8] ) + (  rhs.m_data[ 7] *m_data[12] );
		resultMat[ 5] = ( rhs.m_data[ 4] * m_data[1] ) + ( rhs.m_data[ 5] * m_data[5] ) + ( rhs.m_data[ 6] * m_data[ 9] ) + (  rhs.m_data[ 7] *m_data[13] );
		resultMat[ 6] = ( rhs.m_data[ 4] * m_data[2] ) + ( rhs.m_data[ 5] * m_data[6] ) + ( rhs.m_data[ 6] * m_data[10] ) + (  rhs.m_data[ 7] *m_data[14] );
		resultMat[ 7] = ( rhs.m_data[ 4] * m_data[3] ) + ( rhs.m_data[ 5] * m_data[7] ) + ( rhs.m_data[ 6] * m_data[11] ) + (  rhs.m_data[ 7] *m_data[15] );

		resultMat[ 8] = ( rhs.m_data[ 8] * m_data[0] ) + ( rhs.m_data[ 9] * m_data[4] ) + ( rhs.m_data[10] * m_data[8 ] ) + (  rhs.m_data[11] *m_data[12] );
		resultMat[ 9] = ( rhs.m_data[ 8] * m_data[1] ) + ( rhs.m_data[ 9] * m_data[5] ) + ( rhs.m_data[10] * m_data[9 ] ) + (  rhs.m_data[11] *m_data[13] );
		resultMat[10] = ( rhs.m_data[ 8] * m_data[2] ) + ( rhs.m_data[ 9] * m_data[6] ) + ( rhs.m_data[10] * m_data[10] ) + (  rhs.m_data[11] *m_data[14] );
		resultMat[11] = ( rhs.m_data[ 8] * m_data[3] ) + ( rhs.m_data[ 9] * m_data[7] ) + ( rhs.m_data[10] * m_data[11] ) + (  rhs.m_data[11] *m_data[15] );

		resultMat[12] = ( rhs.m_data[12] * m_data[0] ) + ( rhs.m_data[13] * m_data[4] ) + ( rhs.m_data[14] * m_data[8 ] ) + (  rhs.m_data[15] *m_data[12] );
		resultMat[13] = ( rhs.m_data[12] * m_data[1] ) + ( rhs.m_data[13] * m_data[5] ) + ( rhs.m_data[14] * m_data[9 ] ) + (  rhs.m_data[15] *m_data[13] );
		resultMat[14] = ( rhs.m_data[12] * m_data[2] ) + ( rhs.m_data[13] * m_data[6] ) + ( rhs.m_data[14] * m_data[10] ) + (  rhs.m_data[15] *m_data[14] );
		resultMat[15] = ( rhs.m_data[12] * m_data[3] ) + ( rhs.m_data[13] * m_data[7] ) + ( rhs.m_data[14] * m_data[11] ) + (  rhs.m_data[15] *m_data[15] );

		return Matrix4<T>( resultMat );
	}

	template< typename T>
	T Matrix4<T>::getElement( int row, int col ) const
	{
		return m_rows[row][col];
	}

	template< typename T>
	T Matrix4<T>::getElement( int index ) const
	{
		return m_data[index];
	}

	template< typename T>
	T* Matrix4<T>::toPointer()
	{
		return &m_data[0];
	}

	template< typename T>
	const T* Matrix4<T>::toConstPointer() const
	{
		return &m_data[0];
	}

	template< typename T>
	const Vector4<T>& Matrix4<T>::operator[]( int rowIndex ) const
	{
#if DEBUG
		assert( rowIndex < 4 && rowIndex >= 0 );
#endif

		return m_rows[rowIndex];
	}

	template< typename T>
	Vector4<T>& Matrix4<T>::operator[]( int rowIndex )
	{
#if DEBUG
		assert( rowIndex < 4 && rowIndex >= 0 );
#endif	
		return m_rows[rowIndex];
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::operator*( const Matrix4<T>& rhs ) const
	{
		return postMultiply( rhs );
	}

	template< typename T>
	Matrix4<T>& Matrix4<T>::operator*=( const Matrix4<T>& rhs )
	{
		Matrix4<T> resultMat = postMultiply( rhs );
		*this = resultMat;
		return *this;
	}

	template< typename T>
	Vector4<T> Matrix4<T>::operator*( const Vector4<T>& rhs ) const
	{
		return multiply( rhs );
	}


	template< typename T>
	Vector3<T> Matrix4<T>::operator*(const Vector3<T>& rhs) const
	{
		return multiply( rhs );
	}

	template< typename T>
	Matrix4<T>& Matrix4<T>::operator*=( T val )
	{
		for( int i = 0; i < 16; ++i )
			m_data[i] *= val;

		return *this;
	}



	template< typename T>
	T Matrix4<T>::operator[](Mat4::MATRIX4_ELEMS idx) const
	{
		return m_data[idx];
	}

	template< typename T>
	T& Matrix4<T>::operator[](Mat4::MATRIX4_ELEMS idx)
	{
		return m_data[idx];
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::operator+( const Matrix4<T>& rhs ) const
	{
        Matrix4<T> res;
        for(int i = 0; i < 16; ++i)
            res.m_data[i] += m_data[i] + rhs.m_data[i];

        return res;
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::operator-( const Matrix4<T>& rhs ) const
	{
		
        Matrix4<T> res;
        for(int i = 0; i < 16; ++i)
            res.m_data[i] += m_data[i] - rhs.m_data[i];

       return res;
	}

	template< typename T>
	Matrix4<T>& Matrix4<T>::operator+=( const Matrix4<T>& rhs )
	{
		for( int i = 0; i < 16; ++i )
			m_data[i] += rhs.m_data[i];
		return *this;
	}

	template< typename T>
	Matrix4<T>& Matrix4<T>::operator-=( const Matrix4<T>& rhs )
	{
		for( int i = 0; i < 16; ++i )
			m_data[i] -= rhs.m_data[i];
		return *this;
	}

	template< typename T>
	Matrix4<T> Matrix4<T>::operator*( T val ) const
	{
        Matrix4<T> res;
        for(int i = 0; i < 16; ++i)
            res.m_data[i] = m_data[i] * val;
        return res;
    
	}

	template< typename T>
	bool Matrix4<T>::operator==(const Matrix4<T>& rhs) const
	{
		return equals( rhs );
	}



	template< typename T>
	void Matrix4<T>::setMat3x3( const Matrix3<T>& mat3x3 )
	{
		for( int i = 0; i < 3; ++i )
		{
			m_rows[i][0] = mat3x3[i][0];
			m_rows[i][1] = mat3x3[i][1];
			m_rows[i][2] = mat3x3[i][2];
		}
	}

	template< typename T>
	void Matrix4<T>::setElement( int row, int col, T val )
	{
		m_rows[row][col] = val;
	}



	template< typename T>
	void Matrix4<T>::setElement( int index, T val )
	{
		m_data[index] = val;
	}

	template< typename T>
	void Matrix4<T>::orthographicMatrix2d(T width, T height)
	{
		orthographicMatrix((T) 0.0, (T)width, (T)0.0, (T)height, (T)-1.0, (T)1.0f);
	}


	//////////////////////////////////////////////////////////////////////////
	//\Brief: MatrixX: Operations
	//////////////////////////////////////////////////////////////////////////
	template< typename T>
	MatrixX<T>::MatrixX()
	{
		m_numRows = 0;
		m_numCols = 0;
		m_dataPtr = nullptr;
	}

	template< typename T>
	MatrixX<T>::MatrixX(int numRows, int numCols)
	{
		if (numRows <= 0 || numCols <= 0)
			throw std::exception("Invalid Matrix Dimensions");

		m_numRows = numRows;
		m_numCols = numCols;
		m_dataPtr = new T[m_numRows * m_numCols];
		if (!setIdentity())
			zero();
	}



	template< typename T>
	MatrixX<T>::MatrixX(int numRows, int numCols, const T* dataPtr)
	{
		
		if (numRows <= 0 || numCols <= 0)
			throw std::exception("Invalid Matrix Dimensions");

		m_numRows = numRows;
		m_numCols = numCols;
		m_dataPtr = new T[m_numRows * m_numCols];
		memcpy(m_dataPtr, dataPtr, sizeof(T) * m_numRows *m_numCols);
	}




	template< typename T>
	std::string MatrixX<T>::toString() const
	{
		std::string ret;

		ret += "Rows/Cols [";
		ret += m_numRows;
		ret += ", ";
		ret += m_numCols;
		ret += " ]\n";
		for (int i = 0; i < m_numRows; ++i)
		{
			for (int j = 0; j < m_numCols; ++j)
			{
				char buf[64];
				double elem = (double)getElement(i, j);
				sprintf(buf, "%g ",elem);
				ret += std::string(buf);
			}
			ret += "\n";
		}
		return ret;
	}

	template< typename T>
	MatrixX<T>::MatrixX(int numRows, int numCols, const std::initializer_list<T>& initList)
	{
		if (numRows <= 0 || numCols <= 0)
			throw std::exception("Invalid Matrix Dimensions");
		
		m_numRows = numRows;
		m_numCols = numCols;
		m_dataPtr = new T[m_numRows * m_numCols];
		int count = 0;
		for (auto val : initList)
			m_dataPtr[count++] = val;
	}



	template< typename T>
	MatrixX<T>::MatrixX(int numRows, int numCols, T data)
	{
		m_numRows = numRows;
		m_numCols = numCols;

		m_dataPtr = new T[m_numRows * m_numCols];
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			m_dataPtr[i] = data;
	}


	template< typename T>
	MatrixX<T>::MatrixX(const MatrixX& other)
	{
		m_numRows = other.m_numRows;
		m_numCols = other.m_numCols;

		m_dataPtr = new T[m_numRows * m_numCols];
		memcpy(m_dataPtr, other.m_dataPtr, sizeof(T) * m_numRows * m_numCols);
	}

	template< typename T>
	MatrixX<T>::~MatrixX()
	{
		deAllocate();
	}



	template <typename T>
	template<typename Func>
	void MatrixX<T>::apply(Func fun)
	{
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			m_dataPtr[i] = fun(m_dataPtr[i]);
	}


	template< typename T>
	bool MatrixX<T>::fromString(const std::string& val)
	{
        auto count = Common::CountElements(val);
		if( count < 3 )
			return false;
		auto* data = val.toPointer();
		m_numRows = (int)strtol(data, &data, 10);
		m_numCols = (int)strtol(data, &data, 10);
		allocate( m_numRows, m_numCols );
		for( int i = 0; i < m_numRows; ++i )
		{
			for( int j = 0; j < m_numCols; ++j )
			{
				auto elem = (T) strtod( data, &data );
				m_dataPtr[i][j] = elem;
			}
		}		
		return true;
	}

	template< typename T>
	void MatrixX<T>::swapRow(int row1, int row2)
	{
		int row1Idx = m_numCols * row1;
		int row2Idx = m_numCols * row2;
		for (int i = 0; i < m_numCols; ++i)
			std::swap(m_dataPtr[row1Idx++], m_dataPtr[row2Idx++]);		
	}

	

	template< typename T>
	bool MatrixX<T>::operator==(const MatrixX<T>& rhs) const
	{
		return this->equals(rhs);
	}

	template< typename T>
	bool MatrixX<T>::dimensionEquals(const MatrixX<T>& rhs) const
	{
		return (rhs.m_numRows == m_numRows && rhs.m_numCols == m_numCols);
	}


	template< typename T>
	MatrixX<T> MatrixX<T>::operator+(const MatrixX<T>& rhs) const
	{
		assert(dimensionEquals(rhs));

		MatrixX<T> result(m_numRows, m_numCols);
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			result.m_dataPtr[i] = m_dataPtr[i] + rhs.m_dataPtr[i];
		return result;
	}

	template< typename T>
	MatrixX<T> MatrixX<T>::operator-(const MatrixX<T>& rhs) const
	{
		assert(dimensionEquals(rhs));

		MatrixX<T> result(m_numRows, m_numCols);
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			result.m_dataPtr[i] = m_dataPtr[i] - rhs.m_dataPtr[i];
		return result;
	}


	template< typename T>
	MatrixX<T>& MatrixX<T>::operator+=(const MatrixX<T>& rhs)
	{
		assert(dimensionEquals(rhs));
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			m_dataPtr[i] += rhs.m_dataPtr[i];
		return *this;
	}


	template< typename T>
	MatrixX<T>& MatrixX<T>::operator-=(const MatrixX<T>& rhs)
	{
		assert(dimensionEquals(rhs));
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			m_dataPtr[i] -= rhs.m_dataPtr[i];
		return *this;

	}

	template< typename T>
	MatrixX<T> MatrixX<T>::operator*(T val) const
	{
		MatrixX<T> result(m_numRows, m_numCols);
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			result.m_dataPtr[i] = m_dataPtr[i] * val;
		return result;
	}

	template< typename T>
	MatrixX<T>& MatrixX<T>::operator*=(T val)
	{
		for (int i = 0; i < m_numRows * m_numCols; ++i)
			m_dataPtr[i] = m_dataPtr[i] * val;
		return *this;
	}


	template< typename T>
	VectorX<T> MatrixX<T>::operator*(const VectorX<T>& rhs) const
	{
		assert(rhs.getDimension() == getNumRows());
		VectorX<T> 	result(rhs.getDimension());

		for (int i = 0; i < m_numRows; ++i)
			result[i] = getRow(i).dot(rhs);
		return result;
	}

	template< typename T>
	MatrixX<T>& MatrixX<T>::operator*=(const MatrixX<T>& rhs)
	{
		auto tmp = this->multiply(rhs);
		*this = tmp;
		return *this;
	}

	template< typename T>
	MatrixX<T> MatrixX<T>::operator*(const MatrixX<T>& rhs) const
	{
		return this->multiply(rhs);
	}

	template< typename T>
	T* MatrixX<T>::operator[](int rowIdx)
	{
		assert(rowIdx >= 0);
		assert(rowIdx < m_numRows);
		return &m_dataPtr[rowIdx* m_numCols];
	}

	template< typename T>
	const T* MatrixX<T>::operator[](int rowIdx) const
	{
		assert(rowIdx >= 0);
		assert(rowIdx < m_numRows);
		return &m_dataPtr[rowIdx* m_numCols];
	}

	template< typename T>
	void MatrixX<T>::getColumn(int idx, VectorX<T>& col) const
	{
		assert(col.getDimension() == m_numRows);
		T* startCol = m_dataPtr + idx;
		for (int i = 0; i < m_numRows; ++i)
		{
			col[i] = *startCol;
			startCol += m_numCols;
		}	
	}

	template< typename T>
	void MatrixX<T>::getRow(int idx, VectorX<T>& row) const
	{
		assert(row.getDimension() == m_numCols);
		memcpy(row.m_dataPtr, &m_dataPtr[idx * m_numCols], sizeof(T) * m_numCols);
	}



	template< typename T>
	void MatrixX<T>::setColumn(int col, const VectorX<T>& val)
	{
		assert(col >= 0);
		assert(col < m_numCols);
		assert(m_numRows == val.m_numElements);

		T* startCol = m_dataPtr + col;
		for (int i = 0; i < m_numRows; ++i)
		{
			*startCol = val[i];
			startCol += m_numCols; //update pointer
		}
	}

	template< typename T>
	void MatrixX<T>::setRow(int row, const VectorX<T>& val)
	{
		assert(row >= 0);
		assert(row < m_numRows);
		assert(m_numCols == val.m_numElements);

		memcpy( &m_dataPtr[row * m_numCols], &val.m_dataPtr[0], sizeof(T) * val.m_numElements );
	}

	template< typename T>
	VectorX<T> MatrixX<T>::getRow(int row) const
	{
		return VectorX<T>(&m_dataPtr[row * m_numCols], m_numCols);
	}

	template< typename T>
	VectorX<T> MatrixX<T>::getColumn(int col) const
	{
		VectorX<T> result( m_numRows );
		T* startCol   = m_dataPtr + col;
		for (int i = 0; i < m_numRows; ++i)
		{
			result[i] = *startCol;
			startCol += m_numCols;
		}
		return result;
	}

	template< typename T>
	bool MatrixX<T>::setSubMatrix(int rowStart, int colStart, const MatrixX<T>& val)
	{
		int numRows = val.m_numRows;
		int numCols = val.m_numCols;

		if (((rowStart + numRows) > m_numRows) ||
			((colStart + numCols) > m_numCols))
		{
			return false;
		}

		for (int i = 0; i < numRows; ++i)
			for (int j = 0; j < numCols; ++j)
				setElement(i + rowStart, j + colStart, val.getElement(i, j));
		return true;
	}

	template< typename T>
	MatrixX<T> MatrixX<T>::getSubMatrix(int rowStart, int colStart, int numRows, int numCols) const
	{
		MatrixX<T> result;

		if (((rowStart + numRows) > m_numRows) ||
			((colStart + numCols) > m_numCols))
		{
			return result;
		}
		result.allocate(numRows, numCols);

		for (int i = 0; i < numRows; ++i)
			for (int j = 0; j < numCols; ++j)
				result.setElement(i, j, getElement(i + rowStart, j + colStart));

		return result;
	}


	template< typename T>
	void MatrixX<T>::solveLowerTriangle(const VectorX<T>& b, VectorX<T>& x) const
	{
		x.allocate(b.m_numElements, VectorX<T>::VECX_ALLOC_ZERO );
		int colCount = 1; 	//first solve for Ly = b
		for (int i = 0; i < b.m_numElements; ++i)
		{
			T sum = (T)0.0;
			for (int j = 0; j < colCount; ++j)
				sum += getElement(i, j) * x[j];
			x[i] = (b[i] - sum) / getElement(i, i);
			colCount++;
		}
	}

	template< typename T>
	void MatrixX<T>::solveUpperTriangle(const VectorX<T>& b, VectorX<T>& x) const
	{
		x.allocate(b.m_numElements, VectorX<T>::VECX_ALLOC_ZERO );
		int colCount = 1;
		auto flip = b.m_numElements - 1; //invert row/column iterations
		for (int i = 0; i < b.m_numElements; ++i)
		{
			T sum = (T)0.0;
			int iT = flip - i;
			for (int j = 0; j < colCount; ++j) {
				int jT = flip - j;
				sum += getElement(iT, jT) * x[jT];
			}
			x[iT] = (b[iT] - sum) / getElement(iT, iT);
			colCount++;
		}
	}


	template< typename T>
	bool MatrixX<T>::decomposeLU(MatrixX<T>& L, MatrixX<T>& U) const
	{
		assert(isSquare());
		
		auto dims = m_numRows;
		L.allocate(dims, dims, MATX_RESIZE_IDENTITY );
		U.allocate(dims, dims, MATX_RESIZE_IDENTITY);
		int i, j, k;
		int n = m_numRows;
		T sum = (T) 0.0;
		for (j = 0; j < n; j++)
		{
			for (i = j; i < n; i++)
			{
				sum = (T) 0.0;
				for (k = 0; k < j; k++)
					sum += L.getElement(i, k) * U.getElement(k, j);
				L.m_dataPtr[getElementIdx(i, j)] = getElement(i, j) - sum;
			}

			for (i = j; i < n; i++)
			{
				sum = (T) 0.0;
				for (k = 0; k < j; k++)
					sum += L.getElement(j, k) * U.getElement(k, i);
				if (L.getElement(j, j) < (T) MATRIX_INVERT_EPSILON )
				{
					return false;
				}
				U.m_dataPtr[getElementIdx(j, i)] = (getElement(j, i) - sum) / L.getElement(j, j);
			}
		}

		return true;
	}

	template< typename T>
	bool MatrixX<T>::HouseHolderMatrix(MatrixX<T>& Qn, MatrixX<T>& QA, int curCol /*= 0*/) const
	{
		if (!isSquare())
			return false;

		MatrixX<T> identyMat(m_numRows, m_numRows);
		identyMat.setIdentity();
		auto x = getColumn(curCol);
		auto a = x.length();
		auto reflCol = identyMat.getColumn(curCol) * a;
		auto u = x - reflCol;
		auto v = u.normalized();

		MatrixX<T> vMat(m_numRows, 1, v.toConstPointer());
		MatrixX<T> vMatT = vMat.transposed();

		Qn = identyMat - ((vMat * vMatT) * (T)2.0); //Q1..Q2..QN step of house holder transform
		QA = (Qn) *(*this) ;						   //QN * A will be used as for next step

		return true;
	}

	template< typename T>
	bool MatrixX<T>::decomposeQR(MatrixX<T>& Q, MatrixX<T>& R) const
	{
		if (!isSquare())
			return false;

		int dims = m_numRows;
		int offset = 0;

		MatrixX<float> workMatrix = *this;
		MatrixX<float> QResult(dims, dims); QResult.setIdentity();
		MatrixX<float> Qn;
		MatrixX<float> tmp(dims, dims);

		while (dims >= 2)
		{
			tmp.setIdentity();
			tmp.setSubMatrix(offset, offset, workMatrix.getSubMatrix(offset, offset, dims, dims));
			tmp.HouseHolderMatrix(Qn, workMatrix, offset); 	//generate house holder matrix 
			QResult *= Qn.transposed(); //accumulate result
			dims--;
			offset++;
		}

		Q = QResult;
		R = QResult.transposed() * (*this);
		return true;
	}


	template< typename T>
	bool MatrixX<T>::decomposeSVD( MatrixX<T>& U, MatrixX<T>& V, VectorX<T>& S) const
	{
		int m = m_numRows;
		int n = m_numCols;

		if (m < n)
			return false;

		int flag, i, its, j, jj, k, l, nm;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;

		U = *this;
		S.allocate(m);
		V.allocate(n, n);
		VectorX<double> scratch(n);

		/* Householder reduction to bidiagonal form */
		for (i = 0; i < n; i++)
		{
			/* left-hand reduction */
			l = i + 1;
			scratch[i] = scale * g;
			g = s = scale = 0.0;
			if (i < m)
			{
				for (k = i; k < m; k++)
					scale += std::abs((double)U[k][i]);
				if (scale)
				{
					for (k = i; k < m; k++)
					{
						U[k][i] = (T)((double)U[k][i] / scale);
						s += ((double)U[k][i] * (double)U[k][i]);
					}
					f = (T)U[i][i];
					g = -Math::sign(sqrt(s), f);
					h = f * g - s;
					U[i][i] = (T)(f - g);
					if (i != n - 1)
					{
						for (j = l; j < n; j++)
						{
							for (s = 0.0, k = i; k < m; k++)
								s += ((double)U[k][i] * (double)U[k][j]);
							f = s / h;
							for (k = i; k < m; k++)
								U[k][j] += (T)(f * (double)U[k][i]);
						}
					}
					for (k = i; k < m; k++)
						U[k][i] = (T)((double)U[k][i] * scale);
				}
			}
			S[i] = (T)(scale * g);

			/* right-hand reduction */
			g = s = scale = 0.0;
			if (i < m && i != n - 1)
			{
				for (k = l; k < n; k++)
					scale += std::abs((double)U[i][k]);
				if (scale)
				{
					for (k = l; k < n; k++)
					{
						U[i][k] = (T)((double)U[i][k] / scale);
						s += ((double)U[i][k] * (double)U[i][k]);
					}
					f = (double)U[i][l];
					g = -Math::sign(sqrt(s), f);
					h = f * g - s;
					U[i][l] = (float)(f - g);
					for (k = l; k < n; k++)
						scratch[k] = (double)U[i][k] / h;
					if (i != m - 1)
					{
						for (j = l; j < m; j++)
						{
							for (s = 0.0, k = l; k < n; k++)
								s += ((double)U[j][k] * (double)U[i][k]);
							for (k = l; k < n; k++)
								U[j][k] += (T)(s * scratch[k]);
						}
					}
					for (k = l; k < n; k++)
						U[i][k] = (T)((double)U[i][k] * scale);
				}
			}
			anorm = std::max(anorm, (std::abs((double)S[i]) + std::abs((double)scratch[i])));
		}

		/* accumulate the right-hand transformation */
		for (i = n - 1; i >= 0; i--)
		{
			if (i < n - 1)
			{
				if (g)
				{
					for (j = l; j < n; j++)
						V[j][i] = (T)(((double)U[i][j] / (double)U[i][l]) / g);
					/* double division to avoid underflow */
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < n; k++)
							s += ((double)U[i][k] * (double)V[k][j]);
						for (k = l; k < n; k++)
							V[k][j] += (T)(s * (double)V[k][i]);
					}
				}
				for (j = l; j < n; j++)
					V[i][j] = V[j][i] = 0.0;
			}
			V[i][i] = 1.0;
			g = scratch[i];
			l = i;
		}

		/* accumulate the left-hand transformation */
		for (i = n - 1; i >= 0; i--)
		{
			l = i + 1;
			g = (double)S[i];
			if (i < n - 1)
				for (j = l; j < n; j++)
					U[i][j] = 0.0;
			if (g)
			{
				g = 1.0 / g;
				if (i != n - 1)
				{
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < m; k++)
							s += ((double)U[k][i] * (double)U[k][j]);
						f = (s / (double)U[i][i]) * g;
						for (k = i; k < m; k++)
							U[k][j] += (T)(f * (double)U[k][i]);
					}
				}
				for (j = i; j < m; j++)
					U[j][i] = (T)((double)U[j][i] * g);
			}
			else
			{
				for (j = i; j < m; j++)
					U[j][i] = 0.0;
			}
			++U[i][i];
		}

		/* diagonalize the bidiagonal form */
		for (k = n - 1; k >= 0; k--)
		{                             /* loop over singular values */
			for (its = 0; its < 30; its++)
			{                         /* loop over allowed iterations */
				flag = 1;
				for (l = k; l >= 0; l--)
				{                     /* test for splitting */
					nm = l - 1;
					if (std::abs(scratch[l]) + anorm == anorm)
					{
						flag = 0;
						break;
					}
					if (std::abs((double)S[nm]) + anorm == anorm)
						break;
				}
				if (flag)
				{
				//	c = 0.0;
					s = 1.0;
					for (i = l; i <= k; i++)
					{
						f = s * scratch[i];
						if (std::abs(f) + anorm != anorm)
						{
							g = (double)S[i];
							h = Math::pythag(f, g);
							S[i] = (T)h;
							h = 1.0 / h;
							c = g * h;
							s = (-f * h);
							for (j = 0; j < m; j++)
							{
								y = (double)U[j][nm];
								z = (double)U[j][i];
								U[j][nm] = (T)(y * c + z * s);
								U[j][i] = (T)(z * c - y * s);
							}
						}
					}
				}
				z = (double)S[k];
				if (l == k)
				{                  /* convergence */
					if (z < 0.0)
					{              /* make singular value nonnegative */
						S[k] = (T)(-z);
						for (j = 0; j < n; j++)
							V[j][k] = (-V[j][k]);
					}
					break;
				}
				if (its >= 30)
				{
					//	fprintf(stderr, "No convergence after 30,000! iterations \n");
						//return(0);
					return false;
				}

				/* shift from bottom 2 x 2 minor */
				x = (double)S[l];
				nm = k - 1;
				y = (double)S[nm];
				g = scratch[nm];
				h = scratch[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = Math::pythag(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + Math::sign(g, f))) - h)) / x;

				/* next QR transformation */
				c = s = 1.0;
				for (j = l; j <= nm; j++)
				{
					i = j + 1;
					g = scratch[i];
					y = (double)S[i];
					h = s * g;
					g = c * g;
					z = Math::pythag(f, h);
					scratch[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y = y * c;
					for (jj = 0; jj < n; jj++)
					{
						x = (double)V[jj][j];
						z = (double)V[jj][i];
						V[jj][j] = (T)(x * c + z * s);
						V[jj][i] = (T)(z * c - x * s);
					}
					z = Math::pythag(f, h);
					S[j] = (T)z;
					if (z)
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < m; jj++)
					{
						y = (double)U[jj][j];
						z = (double)U[jj][i];
						U[jj][j] = (T)(y * c + z * s);
						U[jj][i] = (T)(z * c - y * s);
					}
				}
				scratch[l] = 0.0;
				scratch[k] = f;
				S[k] = (T)x;
			}
		}
		return true;
	}



	template< typename T>
	bool MatrixX<T>::decomposeCholesky( MatrixX<T>& L ) const
	{

		if (!isSymmetricPositiveDefinite())
			return false;

		int m = m_numRows;
		int n = m_numCols;

		L.allocate(m, n, MATX_RESIZE_ZERO );

		for (int i = 0; i < n; i++)
		for (int j = 0; j < (i + 1); j++)
		{
			T sum = (T)0.0;
			for (int k = 0; k < j; k++)
				sum += L[i][k] * L[j][k];
			L[i][j] = (i == j) ? sqrt((*this)[i][j] - sum) : 
				((T)1.0 / L[j][j] * ((*this)[i][j] - sum));
		}
		return true;
	}



	template< typename T>
	bool MatrixX<T>::solveLU(const VectorX<T>& b, VectorX<T>& x) const
	{
		MatrixX<T> L, U;
		if (!decomposeLU(L, U))
			return false;
		solveLU(L, U, b, x);
		return true;
	}


	template< typename T>
	bool MatrixX<T>::solveQR(const VectorX<T>& b, VectorX<T>& x) const
	{
		MatrixX<T> Q, R;
		if (!decomposeQR(Q, R))
			return false;
		solveQR(Q, R, b, x);
		return true;
	}

	template< typename T>
	void MatrixX<T>::solveLU(const MatrixX<T>& L, const MatrixX<T>&U, const VectorX<T>& b, VectorX<T>& x)
	{
		VectorX<T> y;
		L.solveLowerTriangle(b, y); //solve for Ly = b
		U.solveUpperTriangle(y, x); //solve for Ux = y		
	}

	template< typename T>
	void MatrixX<T>::solveCholesky(const VectorX<T>& b, VectorX<T>& x) const
	{
		VectorX<T> y;
		solveLowerTriangle(b, y);
		transposed().solveUpperTriangle(y, x);
	}


	template< typename T>
	void MatrixX<T>::solveQR(const MatrixX<T>& Q, const MatrixX<T>& R, const VectorX<T>& b, VectorX<T>& x)
	{
	
		// A = QR
		// Ax = b
		//(QR)x = b
		// y = Qinv * b 
		// Rx = y

		auto y = Q.transposed() * b;
		R.solveUpperTriangle(y, x);		
	}


	template< typename T>
	T MatrixX<T>::trace() const
	{
		return getDiagonal().componentSum();
	}

	template< typename T>
	VectorX<T> MatrixX<T>::getDiagonal() const
	{
		auto min = std::min(m_numRows, m_numCols);
		VectorX<T> result(min);
		for (int i = 0; i < min; ++i)
			result[i] = (*this)[i][i];
		return  result;
	}

	template< typename T>
	void MatrixX<T>::setDiagonal(const VectorX<T>& val)
	{
		auto min = std::min(m_numRows, m_numCols);
		assert(val.getDimension() == min);
		for (int i = 0; i < val.getDimension(); ++i)
			(*this)[i][i] = val[i];
	}

	
	template< typename T>
	int MatrixX<T>::getElementIdx(int row, int col) const
	{
		return (row * m_numCols + col);
	}



	template< typename T>
	void MatrixX<T>::deAllocate()
	{
		m_numRows = 0;
		m_numCols = 0;
		delete[] m_dataPtr;		
		m_dataPtr = nullptr;
	}

	template< typename T>
	void MatrixX<T>::zero()
	{
		for( int i = 0; i < m_numRows * m_numCols; ++i )
			m_dataPtr[i] = (T)0.0;
	}


	template< typename T>
	bool MatrixX<T>::setIdentity()
	{
		if( !isSquare() )
			return false;

		for( int i = 0; i < m_numRows; i++ )
		for( int j = 0; j < m_numCols; j++ )
			(*this)[i][j] = (i == j) ? (T) 1.0 : (T)0.0;
			
		return true;
	}




	template< typename T>
	bool MatrixX<T>::isSquare() const
	{
		return ( m_numRows == m_numCols );
	}

	template< typename T>
	void MatrixX<T>::setElement(  int row, int col, T val )
	{
		m_dataPtr[ row * m_numCols + col] = val;
	}

	template< typename T>
	T MatrixX<T>::getElement( int row, int col ) const
	{
		return m_dataPtr[ row * m_numCols + col];
	}

	template< typename T>
	int MatrixX<T>::getNumRows() const
	{
		return m_numRows;
	}

	template< typename T>
	int MatrixX<T>::getNumColumns() const
	{
		return m_numCols;
	}

	template< typename T>
	MatrixX<T>& MatrixX<T>::operator=( const MatrixX<T>& other )
	{
		allocate( other.m_numRows, other.m_numCols );

		memcpy( m_dataPtr, other.m_dataPtr, sizeof(T) * m_numRows * m_numCols );

		return *this;
	}

	template< typename T>
	void MatrixX<T>::transpose()
	{
		*this = transposed();
	}

	template< typename T>
	MatrixX<T> MatrixX<T>::transposed() const
	{
		MatrixX<T> result(m_numCols, m_numRows);
		for (int i = 0; i < m_numRows; ++i)
		{
			for (int j = 0; j < m_numCols; ++j)
				result.setElement(j, i, getElement(i, j));
		}

		return result;
	}





	template< typename T>
	MatrixX<T> MatrixX<T>::KroneckerProduct(const MatrixX<T>& a, const MatrixX<T>& b)
	{
		MatrixX<T> result(a.m_numRows * b.m_numRows, a.m_numCols * b.m_numCols);
		auto *dataPtr = result.toPointer();
		for (int r1 = 0; r1 < a.m_numRows; ++r1)
		{
			for (int r2 = 0; r2 < b.m_numRows; ++r2)
			{
				//int r3 = r1 + r2;
				for (int c1 = 0; c1 < a.m_numCols; ++c1)
				{
					for (int c2 = 0; c2 < b.m_numCols; ++c2)
					{
						//int c3 = c1 + c2;
						auto val = a[r1][c1] * b[r2][c2];
						*dataPtr++ = a[r1][c1] * b[r2][c2];
					}
				}
			}
		}
		return result;
	}




	template< typename T>
	T* MatrixX<T>::toPointer()
	{
		return &m_dataPtr[0];
	}

	template< typename T>
	const T* MatrixX<T>::toConstPointer() const
	{
		return &m_dataPtr[0];
	}


	template< typename T>
	void MatrixX<T>::setElement( int index, T val )
	{
		m_dataPtr[index] = val;
	}

	template< typename T>
	T MatrixX<T>::getElement( int index ) const
	{
		return m_dataPtr[index];
	}

	template< typename T>
	void MatrixX<T>::allocate( int numRows, int numCols, eMatrixResizeOptions options )
	{
		assert(numRows > 0);
		assert(numCols > 0);
		int newSize = numRows * numCols;
		int oldSize = m_numRows * m_numCols;
		
		//do we need to resize?
		if (oldSize != newSize)
			deAllocate();		
		
		m_numRows = numRows;
		m_numCols = numCols;
		if( !m_dataPtr )
			m_dataPtr = new T[ numRows * numCols ];

		switch (options)
		{
		case MATX_RESIZE_IDENTITY:
			setIdentity();
			break;
		case MATX_RESIZE_ZERO:
			zero();
			break;
		default:
			zero();
			break;
		}
	}

	template< typename T>
	MatrixX<T> MatrixX<T>::pseudoInverse(const MatrixX<T>& U, const MatrixX<T>& V, const VectorX<T>& singVals )
	{
		MatrixX<float> S(U.getNumRows(), U.getNumColumns());
		S.zero();		
		for (int i = 0; i < singVals.getDimension(); ++i)
			S[i][i] = singVals[i] != (T)0.0 ? (T)1.0 / singVals[i] : (T)0.0;
		S = S.transposed();

		return V * S * U.transposed();
	}


	template< typename T>
	bool MatrixX<T>::isDiagonal(T eps /*= (T)FLOAT_EPSILON */) const
	{
		for (int i = 0; i < m_numRows; ++i)
			for (int j = 0; j < m_numCols; ++j)
			{
				if ((i != j) && std::abs((*this)[i][j]) > eps)
					return false;
			}
		return true;
	}




	template< typename T>
	bool MatrixX<T>::isSymmetric(T eps /*= (T)FLOAT_EPSILON*/) const
	{
		if (!isSquare())
			return false;

		for (int i = 0; i < m_numRows; ++i)
			for (int j = 0; j < m_numCols; ++j)
			{
				if (i == j)
					continue;

				auto diff = (*this)[i][j] - (*this)[j][i];
				if (std::abs(diff) >= eps)
					return false;
			}
		return true;
	}

	template< typename T>
	bool MatrixX<T>::isZero(T eps /*= (T) FLOAT_EPSILON */) const
	{
		for (int i = 0; i < m_numRows; ++i)
			for (int j = 0; j < m_numCols; ++j)
			{
				auto val = std::abs((*this)[i][j]);
				if (val >= eps)
					return false;
			}
		return true;
	}

	template< typename T>
	bool MatrixX<T>::isIdentity(T eps /*= (T) FLOAT_EPSILON */) const
	{
		for (int i = 0; i < m_numRows; ++i)
			for (int j = 0; j < m_numCols; ++j)
			{
				auto val = std::abs((*this)[i][j]);
				if (i != j) //off diagonal
				{
					if (val > eps)
						return false;
				}
				else if (std::abs((T)1.0 - val) > eps)  //on diagonal
					return false;
			}
		return true;
	}


	template< typename T>
	bool MatrixX<T>::isSemiPositiveDefinite(T eps /*= (T)FLOAT_EPSILON*/) const
	{
		throw std::exception("Implement me");
		
		return false;
	}



	template< typename T>
	bool MatrixX<T>::isSymmetricPositiveDefinite(T eps /*= (T)FLOAT_EPSILON */) const
	{
		if (!isSymmetric(eps))
			return false;
		return isPositiveDefinite(eps);
	}


	template< typename T>
	bool MatrixX<T>::isPositiveDefinite(T eps /*= (T)FLOAT_EPSILON */) const
	{
		if (!isSquare())
			return false;

		auto m = *this;
		for (int i = 0; i < m_numRows; ++i)
			for (int j = 0; j < m_numCols; j++)
				m[i][j] += (*this)[j][i];

		//all pivots should be non-negative
		for (int i = 0; i < m_numRows; i++)
		{
			for (int j = i; j < m_numCols; j++)
			{
				if (m[j][j] <= eps)
					return false;
			}
			auto d = (T)1.0 / m[i][i];
			for (int j = i + 1; j < m_numCols; j++)
			{
				auto s = d * m[j][i];
				m[j][i] = 0.0f;
				for (int k = i + 1; k < m_numRows; k++)
					m[j][k] -= s * m[i][k];
			}
		}
		return true;
	}



	template< typename T>
	MatrixX<T> MatrixX<T>::multiply( const MatrixX<T>& rhs ) const
	{
		assert( m_numRows == rhs.m_numCols );		

		MatrixX<T> result( m_numRows, rhs.m_numCols );
		T* dstPtr = result.toPointer();
		for (int i = 0; i < m_numRows; ++i)
		{
			auto row = getRow(i);
			for(int j = 0; j < rhs.m_numCols; ++j)
			{
				//multiple each of row 'a' component with each of column 'b component
				//and store it in [row][col] of result matrix
				auto rhsCol = rhs.getColumn(j);
				*dstPtr++ = row.dot(rhsCol);
			}
		}
		return result;
	}

	template< typename T>
	bool MatrixX<T>::equals( const MatrixX<T>& other, T epsilon /*= (T)0.0001 */ ) const
	{
		if( ( m_numRows != other.m_numRows ) || ( m_numCols != other.m_numCols ) )
			return false;

		for( int i = 0; i < ( m_numCols * m_numRows ); ++i )
			if( Abs( m_dataPtr[i] - other.m_dataPtr[i] ) > epsilon )
				return false;

		return true;
	}


	
	template< typename T>	Matrix4<T> Matrix4<T>::IDENTITY = 
					Matrix4<T>(	(T) 1.0, (T) 0.0, (T) 0.0, (T) 0.0,
							(T) 0.0, (T) 1.0, (T) 0.0, (T) 0.0,
							(T) 0.0, (T) 0.0, (T) 1.0, (T) 0.0,
							(T) 0.0, (T) 0.0, (T) 0.0, (T) 1.0  
						);
	template< typename T>	Matrix4<T> Matrix4<T>::ZERO = 
				Matrix4<T>(	(T) 0.0, (T) 0.0, (T) 0.0, (T) 0.0,
							(T) 0.0, (T) 0.0, (T) 0.0, (T) 0.0,
							(T) 0.0, (T) 0.0, (T) 0.0, (T) 0.0,
							(T) 0.0, (T) 0.0, (T) 0.0, (T) 0.0  
		);

	template< typename T>
	Matrix3<T> Matrix3<T>::ZERO = Matrix3<T> (	(T)0.0, (T)0.0, (T)0.0,
												(T)0.0, (T)0.0, (T)0.0,
												(T)0.0, (T)0.0, (T)0.0 );


	 template< typename T>
	 Matrix3<T> Matrix3<T>::IDENTITY = Matrix3<T> ( (T)1.0, (T)0.0, (T)0.0,
												    (T)0.0, (T)1.0, (T)0.0,
												    (T)0.0, (T)0.0, (T)1.0 );

	template< typename T>
	Matrix2<T> Matrix2<T>::ZERO = Matrix2<T> ( (T)0.0, (T)0.0,
									   	       (T)0.0, (T)0.0);
													


	template< typename T>				
	Matrix2<T> Matrix2<T>::IDENTITY  = Matrix2<T> ( (T)1.0, (T)0.0,
												    (T)0.0, (T)1.0 );

	template< typename T>
	Matrix4<T> Matrix4<T>::CreateRotationMatrix(unsigned axis, T radians)
	{
		Matrix4<T> ret;
		
		switch( axis )
		{
		case 0:
			ret.rotationMatrixX( radians );
			break;
		case 1:
			ret.rotationMatrixY( radians );
			break;
		case 2:
			ret.rotationMatrixZ( radians );
			break;
		default:
			assert( false );
		}

		return ret;
	}

  //  namespace ns
   // {
        template <typename T>
        void to_json(IO::JSonObject& obj, const Matrix2<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, Matrix2<T>& v) {
            v.fromJsonObject(obj);
        }

        template <typename T>
        void to_json(IO::JSonObject& obj, const Matrix3<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, Matrix3<T>& v) {
            v.fromJsonObject(obj);
        }

        template <typename T>
        void to_json(IO::JSonObject& obj, const Matrix4<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, Matrix4<T>& v) {
            v.fromJsonObject(obj);
        }

        template <typename T>
        void to_json(IO::JSonObject& obj, const MatrixX<T>& v) {
            obj = v.toJsonObject();
        }
        template <typename T>
        void from_json(const IO::JSonObject& obj, MatrixX<T>& v) {
            v.fromJsonObject(obj);
        }
  //  }

   



