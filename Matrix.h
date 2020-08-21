#ifndef MATRIX_MATH_H
#define MATRIX_MATH_H
#include <sstream>
#include <assert.h>
#include "Vector.h"
#include <Common/StringTools.h>
#include "MathDecl.h"
#include "MathDefs.h" //splitstring

#pragma push_macro("near")
#pragma push_macro("far")

#undef near
#undef far


namespace Math
{
    template<class T>
    inline void sinCos(T angleRadians, T& sinOut, T& cosOut)
    {
        sinOut = sin(angleRadians);
        cosOut = cos(angleRadians);
    }

    template <typename T>
    class Matrix2;
    template <typename T>
    class Matrix3;
    template <typename T>
    class Matrix4;
    
    namespace Mat2
	{
		enum MATRIX2_ELEMS
		{
			m_11,
			m_12,
			m_21,
			m_22
		};
	}

	//////////////////////////////////////////////////////////////////////////
	//Matrix2
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Matrix2
	{
	public:


		Matrix2();
		Matrix2(const Vector2<T>& row1,
		        const Vector2<T>& row2);


		Matrix2(T _11, T _12,
		        T _21, T _22);

		Matrix2(const T* dataPtr);

		Matrix2(const Matrix2<T>& other);

		~Matrix2()
		{
		};

		//////////////////////////////////////////////////////////////////////////
		// \brief: Sets to identity matrix
		//////////////////////////////////////////////////////////////////////////
		void setIdentity();

		//////////////////////////////////////////////////////////////////////////
		// \brief: Clears Matrix
		//////////////////////////////////////////////////////////////////////////
		void zero();


		//////////////////////////////////////////////////////////////////////////
		// \brief: Set's input matrix to inverted of this matrix
		//////////////////////////////////////////////////////////////////////////
		bool invert(Matrix2<T>& invertedMat) const;


		bool isIdentity() const
		{
			return this->equals(Matrix2<T>());
		}


		//////////////////////////////////////////////////////////////////////////
		// \brief: Invert matrix
		//////////////////////////////////////////////////////////////////////////
		Matrix2<T> inverted() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix 
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrix(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the determinant
		//////////////////////////////////////////////////////////////////////////
		T determinant() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the trace
		//////////////////////////////////////////////////////////////////////////
		T trace() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose this matrix
		//////////////////////////////////////////////////////////////////////////
		void transpose();

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose input matrix 
		//////////////////////////////////////////////////////////////////////////
		void transpose(Matrix2<T>& transMat) const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Scale this matrix by input parameter
		//////////////////////////////////////////////////////////////////////////
		void scale(T scaleVal);
        void setScale(const Vector2<T>& scale);
       
		//////////////////////////////////////////////////////////////////////////
		// \brief: Equality test function
		//////////////////////////////////////////////////////////////////////////
		bool equals(const Matrix2<T>& other, T epsilon = (T)0.0001) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply a Vector with this matrix, returns the transformed vector
		//////////////////////////////////////////////////////////////////////////
		Vector2<T> multiply(const Vector2<T>& val) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = A * B
		//////////////////////////////////////////////////////////////////////////
		Matrix2<T> postMultiply(const Matrix2<T>& rhs) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = B * A
		//////////////////////////////////////////////////////////////////////////
		Matrix2<T> preMultiply(const Matrix2<T>& rhs) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns string representation
		//////////////////////////////////////////////////////////////////////////
		std::string toString() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool fromString(const std::string& val);

        bool fromJsonObject(const IO::JSonObject& obj)
        {
            for (int i = 0; i < 4; ++i)
            {
                auto name = "_0" + std::to_string(i);
                m_data[i] = obj.value( name.c_str(), m_data[i]);                
            }
            return true;
        }

        
        IO::JSonObject toJsonObject() const
        {
            IO::JSonObject result;     
            for (int i = 0; i < 4; ++i)
            {
                auto name = "_0" + std::to_string(i);
                result[name] = m_data[i];
            }        
            return result;
        }




		//////////////////////////////////////////////////////////////////////////
		//\brief: Custom function operates on each element
		//////////////////////////////////////////////////////////////////////////
		template <typename Func>
		void apply(Func fun);

		void setElement(int row, int col, T val);
		void setElement(int index, T val);

		T getElement(int row, int col) const;
		T getElement(int index) const;

		T* toPointer();
		const T* toConstPointer() const;
        Matrix3<T> toMatrix3() const;

		T operator [](Mat2::MATRIX2_ELEMS idx) const;
		T& operator [](Mat2::MATRIX2_ELEMS idx);


		const Vector2<T>& operator [](int rowIndex) const;
		Vector2<T>& operator [](int rowIndex);

		Matrix2<T> operator *(const Matrix2<T>& rhs) const;
		Matrix2<T>& operator *=(const Matrix2<T>& rhs);
		Vector2<T> operator *(const Vector2<T>& rhs) const;
		Matrix2<T>& operator *=(T val);
		Matrix2<T> operator *(T val) const;
		Matrix2<T>& operator +=(const Matrix2<T>& rhs);
		Matrix2<T>& operator -=(const Matrix2<T>& rhs);
		Matrix2<T> operator +(const Matrix2<T>& rhs) const;
		Matrix2<T> operator -(const Matrix2<T>& rhs) const;

		Matrix2<T>& operator =(const Matrix2<T>& rhs);
		bool operator ==(const Matrix2<T>& rhs) const;


		static Matrix2<T> IDENTITY;
		static Matrix2<T> ZERO;


#pragma warning ( push )
#pragma warning (disable : 4201 )
	public:
		union
		{
			struct
			{
				Vector2<T> m_rows[2];
			};

			T m_data[4];
		};//union
#pragma warning ( pop )
	};

 
  
   

	namespace Mat3
	{
		enum MATRIX3_ELEMS
		{
			m_11,
			m_12,
			m_13,
			m_21,
			m_22,
			m_23,
			m_31,
			m_32,
			m_33
		};
	}

	//////////////////////////////////////////////////////////////////////////
	//Matrix3
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Matrix3
	{
	public:


		Matrix3();
		Matrix3(const Vector3<T>& row1,
		        const Vector3<T>& row2,
		        const Vector3<T>& row3);

		Matrix3(const Matrix2<T>& rhs)
		{
			m_rows[0] = Vector3<T>(rhs[0], (T)0);
			m_rows[1] = Vector3<T>(rhs[1], (T)0);
			m_rows[2] = Vector3<T>((T)0, (T)0, (T) 1 );
		}
		
		Matrix3(T _11, T _12, T _13,
		        T _21, T _22, T _23,
		        T _31, T _32, T _33);

		Matrix3(const T* dataPtr);

		Matrix3(const Matrix3<T>& other);

		~Matrix3()
		{
		};


		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns a String representation
		//////////////////////////////////////////////////////////////////////////
		std::string toString() const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool fromString(const std::string& val);


		//////////////////////////////////////////////////////////////////////////
		// \brief: Sets to identity matrix
		//////////////////////////////////////////////////////////////////////////
		void setIdentity();

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns if this is an identity matrix
		//////////////////////////////////////////////////////////////////////////
		bool isIdentity() const
		{
			return this->equals(Matrix3<T>());
		}

		//////////////////////////////////////////////////////////////////////////
		// \brief: Clears Matrix
		//////////////////////////////////////////////////////////////////////////
		void zero();


		//////////////////////////////////////////////////////////////////////////
		// \brief: Set's input matrix to inverted of this matrix
		//////////////////////////////////////////////////////////////////////////
		bool invert(Matrix3<T>& invertedMat) const;
				
		//////////////////////////////////////////////////////////////////////////
		// \brief: Invert matrix
		//////////////////////////////////////////////////////////////////////////
		Matrix3<T> inverted() const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix  X
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixX(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix Y
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixY(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix  Z
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixZ(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the determinant
		//////////////////////////////////////////////////////////////////////////
		T determinant() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the trace
		//////////////////////////////////////////////////////////////////////////
		T trace() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose this matrix
		//////////////////////////////////////////////////////////////////////////
		void transpose();

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose input matrix 
		//////////////////////////////////////////////////////////////////////////
		void transpose(Matrix3<T>& transMat) const;


		Matrix3<T> transposed() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Return the adjoint matrix of this matrix
		//////////////////////////////////////////////////////////////////////////
		void adjoint(Matrix3<T>& adjointMat) const;


		bool equals(const Matrix3<T>& other, T epsilon = (T)0.0001) const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Scale this matrix by input parameter
		//////////////////////////////////////////////////////////////////////////
		void scale(T scaleVal);
		void setScale(const Vector3<T>& val);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply a Vector with this matrix, returns the transformed vector
		//////////////////////////////////////////////////////////////////////////
		Vector3<T> multiply(const Vector3<T>& val) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply a Vector with this matrix, returns the transformed vector
		//////////////////////////////////////////////////////////////////////////
		Vector3<T> multiplyTrans(const Vector3<T>& val) const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = A * B
		//////////////////////////////////////////////////////////////////////////
		Matrix3<T> postMultiply(const Matrix3<T>& rhs) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = B * A
		//////////////////////////////////////////////////////////////////////////
		Matrix3<T> preMultiply(const Matrix3<T>& rhs) const;


		void        setElement(int row, int col, T val);
		void        setElement(int index, T val);

		T           getElement(int row, int col) const;
		T           getElement(int index) const;

		T*          toPointer();
		const T*    toConstPointer() const;


        Matrix2<T>  toMatrix2() const;

        Matrix4<T>  toMatrix4() const;

		//////////////////////////////////////////////////////////////////////////
		//\brief: Custom function operates on each element
		//////////////////////////////////////////////////////////////////////////
		template <typename Func>
		void apply(Func fun);


		bool operator ==(const Matrix3<T>& rhs) const;

		T operator [](Mat3::MATRIX3_ELEMS idx) const;
		T& operator [](Mat3::MATRIX3_ELEMS idx);

		const Vector3<T>& operator [](int rowIndex) const;
		Vector3<T>& operator [](int rowIndex);

		Matrix3<T> operator *(const Matrix3<T>& rhs) const;
		Matrix3<T>& operator *=(const Matrix3<T>& rhs);
		Vector3<T> operator *(const Vector3<T>& rhs) const;
		Matrix3<T>& operator *=(T val);
		Matrix3<T> operator *(T val) const;
		Matrix3<T>& operator +=(const Matrix3<T>& rhs);
		Matrix3<T>& operator -=(const Matrix3<T>& rhs);
		Matrix3<T> operator +(const Matrix3<T>& rhs) const;
		Matrix3<T> operator -(const Matrix3<T>& rhs) const;

		Matrix3<T>& operator =(const Matrix3<T>& rhs);


        bool fromJsonObject(const IO::JSonObject& obj)
        {
            for (int i = 0; i < 9; ++i)
            {
                auto name = "_0" + std::to_string(i);
                m_data[i] = obj.value( name.c_str(), m_data[i]);
            }
            return true;
        }


        IO::JSonObject toJsonObject() const
        {
            IO::JSonObject result;
            for (int i = 0; i < 9; ++i)
            {
                auto name = "_0" + std::to_string(i);
                result[name] = m_data[i];
            }
            return result;
        }


		static Matrix3<T> IDENTITY;
		static Matrix3<T> ZERO;


	public:
#pragma warning ( push )
#pragma warning (disable : 4201 )
		union
		{
			struct
			{
				Vector3<T> m_rows[3];
			};

			T m_data[9];
		};//union
#pragma warning ( pop )
	};

  
  

	namespace Mat4
	{
		enum MATRIX4_ELEMS
		{
			m_11,
			m_12,
			m_13,
			m_14,
			m_21,
			m_22,
			m_23,
			m_24,
			m_31,
			m_32,
			m_33,
			m_34,
			m_41,
			m_42,
			m_43,
			m_44
		};
	}

	//////////////////////////////////////////////////////////////////////////
	//Matrix4
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class Matrix4
	{
	public:


		Matrix4();

		Matrix4(const Vector4<T>& row1,
		        const Vector4<T>& row2,
		        const Vector4<T>& row3,
		        const Vector4<T>& row4);		

		Matrix4( const Matrix3<T>& rhs )
			
		{
			m_rows[0] = Vector4<T>(rhs[0], (T)0.0);
			m_rows[1] = Vector4<T>(rhs[1], (T)0.0);
			m_rows[2] = Vector4<T>(rhs[2], (T)0.0);
			m_rows[3] = Vector4<T>((T)0, (T)0, (T)0, (T)1);
		}


		Matrix4(T _11, T _12, T _13, T _14,
		        T _21, T _22, T _23, T _24,
		        T _31, T _32, T _33, T _34,
		        T _41, T _42, T _43, T _44);


		Matrix4(const T* dataPtr);


		Matrix4(const Matrix4<T>& other);


		~Matrix4()
		{
		}


		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns if this is an identity matrix
		//////////////////////////////////////////////////////////////////////////
		bool isIdentity() const
		{
			return this->equals(Matrix4<T>());
		}
		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns a String representation
		//////////////////////////////////////////////////////////////////////////
		std::string toString() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool fromString(const std::string& val);


		//////////////////////////////////////////////////////////////////////////
		// \brief: Sets identity matrix
		//////////////////////////////////////////////////////////////////////////
		void setIdentity();

		void zero();

		void setTranslation(const Vector3<T>& translation);


		void setDirection(const Vector3<T>& direction);
		void setUp(const Vector3<T>& up);
		void setRight(const Vector3<T>& right);

		Vector3<T> getDirection() const;
		Vector3<T> getUp() const;
		Vector3<T> getRight() const;
		Vector3<T> getScale(bool normalize = false) const;
		Vector3<T> getTranslation() const;

		///////////////////////////////////////////////////////////////////////////
		// \brief: Set rotational part of this 4 by 4 mat
		//////////////////////////////////////////////////////////////////////////
		void setMat3x3(const Matrix3<T>& mat3x3);

		void setScale(const Vector3<T>& scale);

		template <class U>
		Matrix4<U> convertTo() const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix  X
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixX(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix Y
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixY(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Set rotation matrix  Z
		//////////////////////////////////////////////////////////////////////////
		void rotationMatrixZ(T angleRadians);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the determinant
		//////////////////////////////////////////////////////////////////////////
		T determinant() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns the trace
		//////////////////////////////////////////////////////////////////////////
		T trace() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose this matrix
		//////////////////////////////////////////////////////////////////////////
		void transpose();

		//////////////////////////////////////////////////////////////////////////
		//\Brief: return a transposed version of this matrix
		//////////////////////////////////////////////////////////////////////////
		Matrix4<T> transposed() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transpose input matrix 
		//////////////////////////////////////////////////////////////////////////
		void transpose(Matrix4<T>& transMat) const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Set input matrix to the invert of this matrix
		//////////////////////////////////////////////////////////////////////////
		bool invert(Matrix4<T>& inverseOut) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Returns an inverted matrix
		//////////////////////////////////////////////////////////////////////////
		Matrix4<T> inverted() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Creates a perspective matrix
		//////////////////////////////////////////////////////////////////////////
		void perspectiveMatrix(T angle, T viewportRatio, T nearPlane, T farPlane);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns a orthographic matrix
		//////////////////////////////////////////////////////////////////////////
		void orthographicMatrix(T left, T right, T bottom, T top, T znear, T zfar);

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Returns a 2d ortho matrix
		//////////////////////////////////////////////////////////////////////////
		void orthographicMatrix2d(T width, T height);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Return the adjoint matrix of this matrix
		//////////////////////////////////////////////////////////////////////////
		void adjoint(Matrix4<T>& adjointMat);

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Test for equality
		//////////////////////////////////////////////////////////////////////////
		bool equals(const Matrix4<T>& other, T epsilon = (T)0.0001) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: 
		//////////////////////////////////////////////////////////////////////////
		Matrix3<T> toMatrix3() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Scale this matrix by input parameter
		//////////////////////////////////////////////////////////////////////////
		void scale(T scaleVal);


        void setScale(const Vector4<T>& val);

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply a Vector with this matrix, returns the transformed vector
		//////////////////////////////////////////////////////////////////////////
		Vector4<T> multiply(const Vector4<T>& val) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply vector3 with this transform
		//////////////////////////////////////////////////////////////////////////
		Vector3<T> multiply(const Vector3<T>& val) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Transform a normal
		//////////////////////////////////////////////////////////////////////////
		Vector3<T> transformNormal(const Vector3<T>& val) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = A * B
		//////////////////////////////////////////////////////////////////////////
		Matrix4<T> postMultiply(const Matrix4<T>& rhs) const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply this matrix with another one, return multiplied matrix(as copy)
		//		   A = B * A
		//////////////////////////////////////////////////////////////////////////
		Matrix4<T> preMultiply(const Matrix4<T>& rhs) const;


		void lookAt(const Vector3<T>& eye, const Vector3<T>& target, const Vector3<T>& up);


        bool fromJsonObject(const IO::JSonObject& obj)
        {
            for (int i = 0; i < 16; ++i)
            {
                auto name = "_0" + std::to_string(i);
                m_data[i] = obj.value(name.c_str(), m_data[i]);
            }
            return true;
        }


        IO::JSonObject toJsonObject() const
        {
            IO::JSonObject result;
            for (int i = 0; i < 16; ++i)
            {
                auto name = "_0" + std::to_string(i);
                result[name] = m_data[i];
            }
            return result;
        }


		void setElement(int row, int col, T val);
		void setElement(int index, T val);

		T getElement(int row, int col) const;
		T getElement(int index) const;

		T* toPointer();
		const T* toConstPointer() const;

		//////////////////////////////////////////////////////////////////////////
		//\brief: Custom function operates on each element
		//////////////////////////////////////////////////////////////////////////
		template <typename Func>
		void apply(Func fun);


		T operator [](Mat4::MATRIX4_ELEMS idx) const;
		T& operator [](Mat4::MATRIX4_ELEMS idx);

		const Vector4<T>& operator [](int rowIndex) const;
		Vector4<T>& operator [](int rowIndex);

		Matrix4<T> operator *(const Matrix4<T>& rhs) const;
		Matrix4<T>& operator *=(const Matrix4<T>& rhs);
		Vector4<T> operator *(const Vector4<T>& rhs) const;
		Vector3<T> operator *(const Vector3<T>& rhs) const;

		Matrix4<T>& operator *=(T val);
		Matrix4<T> operator *(T val) const;
		Matrix4<T>& operator +=(const Matrix4<T>& rhs);
		Matrix4<T>& operator -=(const Matrix4<T>& rhs);
		Matrix4<T> operator +(const Matrix4<T>& rhs) const;
		Matrix4<T> operator -(const Matrix4<T>& rhs) const;

		bool operator ==(const Matrix4<T>& rhs) const;

		Matrix4<T>& operator =(const Matrix4<T>& rhs);


        static Matrix4<T> CreateScaleMatrix(const Vector3<T>& xyz);
        static Matrix4<T> CreateTranslationMatrix(const Vector3<T>& xyz);
		static Matrix4<T> CreateRotationMatrix(unsigned axis, T radians);
		static Matrix4<T> Ortho2DMatrix(T width, T height );
		static Matrix4<T> Perspective(T fovy, T aspect, T nearVal, T farVal);
		static Matrix4<T> Frustum(T left, T right, T bottom, T top, T nearVal, T farVal);

	public:
#pragma warning ( push )
#pragma warning (disable : 4201 )
		union
		{
			struct
			{
				Vector4<T> m_rows[4];
			};

			T m_data[16];
		};//union
#pragma warning ( pop )

		static Matrix4<T> IDENTITY;
		static Matrix4<T> ZERO;
	};

  

	template <typename T>
	static inline T sign(T a, T b)
	{
		return (b >= (T)0.0) ? std::abs(a) : -std::abs(a);
	}


	template <typename T>
	static inline T pythag(T a, T b)
	{
		double at = abs((double)a);
		double bt = abs((double)b);

		if (at > bt)
		{
			auto ct = bt / at;
			return (T)(at * sqrt(1.0 + ct * ct));
		}
		if (bt != (T)0.0)
		{
			auto ct = at / bt;
			return (T)(bt * sqrt(1.0 + ct * ct));
		}
		return (T)0.0;
	}


	//////////////////////////////////////////////////////////////////////////
	//MatrixX
	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	class MatrixX
	{
	public:

		enum eMatrixResizeOptions
		{
			MATX_RESIZE_NONE = 0, //do nothing
			MATX_RESIZE_ZERO = 1, //zero data
			MATX_RESIZE_IDENTITY = 2 //make matrix identity
		};

		MatrixX();
		MatrixX(const MatrixX& other);


		//////////////////////////////////////////////////////////////////////////
		// \brief: Constructor, if matrix is square, generates identity matrix,
		//         else sets' elements to zero
		//////////////////////////////////////////////////////////////////////////
		explicit MatrixX(int numRows, int numCols);

		explicit MatrixX(int numRows, int numCols, T data);

		explicit MatrixX(int numRows, int numCols, const T* dataPtr);

		explicit MatrixX(int numRows, int numCols, const std::initializer_list<T>& initList);


		~MatrixX();


		//////////////////////////////////////////////////////////////////////////
		// \brief: Returns a String representation
		//////////////////////////////////////////////////////////////////////////
		std::string toString() const;

		//////////////////////////////////////////////////////////////////////////
		// \brief: Convert from a string
		//////////////////////////////////////////////////////////////////////////
		bool fromString(const std::string& val);

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Returns a 1 dimensional index for the row and column
		//////////////////////////////////////////////////////////////////////////
		int getElementIdx(int row, int col) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Allocate new data from the heap
		//////////////////////////////////////////////////////////////////////////
		void allocate(int numRows, int numCols, eMatrixResizeOptions = MATX_RESIZE_NONE);

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Returns any memory allocated back to the heap
		//////////////////////////////////////////////////////////////////////////
		void deAllocate();
		//////////////////////////////////////////////////////////////////////////
		//\Brief: Fill with zeros;
		//////////////////////////////////////////////////////////////////////////
		void zero();

		//////////////////////////////////////////////////////////////////////////
		//\brief: Set to identity matrix
		//////////////////////////////////////////////////////////////////////////
		bool setIdentity();

		//////////////////////////////////////////////////////////////////////////
		//\Brief:Return true if this matrix is square
		//////////////////////////////////////////////////////////////////////////
		bool isSquare() const;

		//////////////////////////////////////////////////////////////////////////
		//\brief: Returns if this is a identity matrix
		//////////////////////////////////////////////////////////////////////////
		bool isIdentity(T eps = (T)FLOAT_EPSILON) const;

		//////////////////////////////////////////////////////////////////////////
		//\ Brief:Returns true if all elements are zero, false otherwise
		//////////////////////////////////////////////////////////////////////////
		bool isZero(T eps = (T)FLOAT_EPSILON) const;

		//////////////////////////////////////////////////////////////////////////
		//\ Brief:Returns true if this matrix is symmetric
		//////////////////////////////////////////////////////////////////////////
		bool isSymmetric(T eps = (T)FLOAT_EPSILON) const;

		//////////////////////////////////////////////////////////////////////////
		// \Brief: Is positive definite
		//////////////////////////////////////////////////////////////////////////
		bool isPositiveDefinite(T eps = (T)FLOAT_EPSILON) const;


		//////////////////////////////////////////////////////////////////////////
		// \Brief: Is positive definite
		//////////////////////////////////////////////////////////////////////////
		bool isSemiPositiveDefinite(T eps = (T)FLOAT_EPSILON) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Return true if this matrixd is symmetric and positive definite
		//////////////////////////////////////////////////////////////////////////
		bool isSymmetricPositiveDefinite(T eps = (T)FLOAT_EPSILON) const;

		//////////////////////////////////////////////////////////////////////////
		//\brief: Returns true if this matrix is a diagonal matrix, false otherwise
		//////////////////////////////////////////////////////////////////////////
		bool isDiagonal(T eps = (T)FLOAT_EPSILON) const;


		//////////////////////////////////////////////////////////////////////////
		//\Brief: Return true if this matrix & the rhs have equal dimensions
		//////////////////////////////////////////////////////////////////////////
		bool dimensionEquals(const MatrixX<T>& rhs) const;

		int getNumColumns() const;
		int getNumRows() const;

		T getElement(int row, int col) const;
		void setElement(int row, int col, T val);

		T getElement(int index) const;
		void setElement(int index, T val);

		void swapRow(int r1, int r2);

		void transpose();
		MatrixX<T> transposed() const;

		T* toPointer();
		const T* toConstPointer() const;

		void setRow(int row, const VectorX<T>& val);
		void setColumn(int col, const VectorX<T>& val);

		VectorX<T> getRow(int row) const;
		VectorX<T> getColumn(int col) const;

		void getRow(int idx, VectorX<T>& row) const;
		void getColumn(int idx, VectorX<T>& col) const;

		void setDiagonal(const VectorX<T>& val);
		VectorX<T> getDiagonal() const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: returns the trace of this matrix
		//////////////////////////////////////////////////////////////////////////
		T trace() const;


		//////////////////////////////////////////////////////////////////////////
		//\brief: Custom function operates on each element
		//////////////////////////////////////////////////////////////////////////
		template <typename Func>
		void apply(Func fun);


		MatrixX<T> getSubMatrix(int rowStart, int colStart, int numRows, int numCols) const;
		bool setSubMatrix(int rowStart, int colStart, const MatrixX<T>& val);

		const T* operator[](int rowIdx) const;
		T* operator[](int rowIdx);

		MatrixX<T> operator *(const MatrixX<T>& rhs) const;
		MatrixX<T>& operator *=(const MatrixX<T>& rhs);
		VectorX<T> operator *(const VectorX<T>& rhs) const;

		MatrixX<T>& operator *=(T val);
		MatrixX<T> operator *(T val) const;
		MatrixX<T>& operator +=(const MatrixX<T>& rhs);
		MatrixX<T>& operator -=(const MatrixX<T>& rhs);
		MatrixX<T> operator +(const MatrixX<T>& rhs) const;
		MatrixX<T> operator -(const MatrixX<T>& rhs) const;

		MatrixX<T>& operator =(const MatrixX<T>& other);
		bool operator ==(const MatrixX<T>& rhs) const;


		//////////////////////////////////////////////////////////////////////////
		//\Brief : Solve Ax = b, where A == this
		//////////////////////////////////////////////////////////////////////////
		bool solveQR(const VectorX<T>& b, VectorX<T>& x) const;


		//////////////////////////////////////////////////////////////////////////
		//\Brief : Solve Ax = b, where A == this
		//////////////////////////////////////////////////////////////////////////
		void solveCholesky(const VectorX<T>& b, VectorX<T>& x) const;


		//////////////////////////////////////////////////////////////////////////
		//\Brief : Solve Ax = b, where A == this
		//////////////////////////////////////////////////////////////////////////
		bool solveLU(const VectorX<T>& b, VectorX<T>& x) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief:Given a L & U Matrix solve for Ax = b
		//////////////////////////////////////////////////////////////////////////
		static void solveLU(const MatrixX<T>& L, const MatrixX<T>& U,
		                    const VectorX<T>& b, VectorX<T>& x);
		//////////////////////////////////////////////////////////////////////////
		//\Brief:Given a L & U Matrix solve for Ax = b
		//////////////////////////////////////////////////////////////////////////
		static void solveQR(const MatrixX<T>& Q, const MatrixX<T>& R,
		                    const VectorX<T>& b, VectorX<T>& x);


		//////////////////////////////////////////////////////////////////////////
		// \Brief: Calculates the pseudo inverse: A+ = VS+U*
		//////////////////////////////////////////////////////////////////////////
		static MatrixX<T> pseudoInverse(const MatrixX<T>& U, const MatrixX<T>& V,
		                                const VectorX<T>& S);


		//////////////////////////////////////////////////////////////////////////
		//\Brief: Returns the kronecker product of 2 matrices
		//////////////////////////////////////////////////////////////////////////
		static MatrixX<T> KroneckerProduct(const MatrixX<T>& a, const MatrixX<T>& b);


		//////////////////////////////////////////////////////////////////////////
		//\Brief: Solve Ax = b, using backward substitution
		//////////////////////////////////////////////////////////////////////////
		void solveUpperTriangle(const VectorX<T>& b, VectorX<T>& x) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Solve Ax = b, using forward substitution
		//////////////////////////////////////////////////////////////////////////
		void solveLowerTriangle(const VectorX<T>& b, VectorX<T>& x) const;

		//////////////////////////////////////////////////////////////////////////
		//\brief: Decompose matrix in a upper and lower part, this is 
		//		  based on Crout's Algorithm
		//////////////////////////////////////////////////////////////////////////
		bool decomposeLU(MatrixX<T>& L, MatrixX<T>& U) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Generate a householder reflection matrix 
		//////////////////////////////////////////////////////////////////////////
		bool HouseHolderMatrix(MatrixX<T>& Qn, MatrixX<T>& QA, int curCol = 0) const;


		//////////////////////////////////////////////////////////////////////////
		// \Brief: Decompose this matrix using cholesky factorization
		//		   Only the lower triangular matrix is set since:
		//		   Acholesky = L * Lt
		//////////////////////////////////////////////////////////////////////////
		bool decomposeCholesky(MatrixX<T>& L) const;

		//////////////////////////////////////////////////////////////////////////
		//\Brief: Perform QR decomposition on this matrix
		//////////////////////////////////////////////////////////////////////////
		bool decomposeQR(MatrixX<T>& Q, MatrixX<T>& R) const;


		//////////////////////////////////////////////////////////////////////////
		//\Brief: Decompose this matrix in left( U ) right ( V ) and singular (S) values
		//		  N.B. 'V' should be transposed after calling this function
		//////////////////////////////////////////////////////////////////////////
		bool decomposeSVD(MatrixX<T>& U, MatrixX<T>& V, VectorX<T>& S) const;


		//////////////////////////////////////////////////////////////////////////
		// \brief: Multiply NxM matrix by IxJ matrix, return result which is
		//		   a NxJ matrix
		//////////////////////////////////////////////////////////////////////////
		MatrixX<T> multiply(const MatrixX<T>& rhs) const;

		//////////////////////////////////////////////////////////////////////////
		// \Brief: Compare this matrix with another one
		//////////////////////////////////////////////////////////////////////////
		bool equals(const MatrixX<T>& other, T eps = (T)FLOAT_EPSILON) const;


        bool fromJsonObject(const IO::JSonObject& obj)
        {
            int rows = obj.value("Rows", 0 );
            int cols = obj.value("Cols", 0 );
            if (!rows || !cols)
                return false;

            allocate(rows, cols);

            for (int row = 0; row < getNumRows(); ++row)
            for (int col = 0; col < getNumColumns(); ++col)
            {
                    std::string name = "_" + std::to_string(row) + "_" + std::to_string(col);                                 
                    T tmp =  obj.value( name.c_str(), (T) 0.0 );
                    setElement( row, col, tmp );
            }
            return true;
        }


        IO::JSonObject toJsonObject() const
        {
            IO::JSonObject result;
            result["Rows"] = m_numRows;
            result["Cols"] = m_numCols;
            for (int row = 0; row < getNumRows(); ++row)
            for (int col = 0; col < getNumColumns(); ++col)
            {
                    std::string name = "_" + std::to_string(row) + "_" + std::to_string(col);
                    result[name.c_str()] = getElement(row, col);
            }
            return result;
        }


	public:

		T* m_dataPtr;
		int m_numCols;
		int m_numRows;
	};


	typedef Matrix2<float> Matrix2f;
	typedef Matrix2<double> Matrix2d;
	typedef Matrix3<float> Matrix3f;
	typedef Matrix3<double> Matrix3d;
	typedef Matrix4<float> Matrix4f;
	typedef Matrix4<double> Matrix4d;


#include "Matrix.inl"
};


#pragma pop_macro ("near")
#pragma pop_macro ("far")

#endif
