#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
/*qubit中的複數的虛數部份以及實數部分都介於-1~1
 * 因此使用long long配合定點數可以有效的解決運算中的浮點誤差問題
 * 以2^62為縮放基準
 */

class Qubit_Simulation{

	public:

	static constexpr const double Double_Epsilon=1e-9;
	static constexpr const long long int Fixed_Point=(1LL<<62);
	static constexpr const int Fixed_shift=62;
	static constexpr const double Root_half=0.70710678118;
	static constexpr const double Pi=3.141592658979;
	static constexpr const unsigned long long int Fixed_Pi=14488038916154245120ULL;

	static long long int FixedPoint(const double input){

		return static_cast<long long>(input*Fixed_Point);

	}

	static double FixedPointToDouble(const long long int input){

		double result=1.0*input/Fixed_Point;

		return result;

	}

	struct FixedComplex{

		public:

		long long int Real=0,Imaginary=0;

			FixedComplex():Real(0),Imaginary(0){};
			explicit FixedComplex(double A):Real(FixedPoint(A)),Imaginary(0){};
			explicit FixedComplex(double A,double B):Real(FixedPoint(A)),Imaginary(FixedPoint(B)){};
			FixedComplex(long long int A):Real(A),Imaginary(0){};
			FixedComplex(long long int A,long long int B):Real(A),Imaginary(B){};

			friend FixedComplex operator+(FixedComplex& c,const double d)noexcept{
				return FixedComplex(c.Real+FixedPoint(d),c.Imaginary);
			}
			friend FixedComplex operator+(const double d,FixedComplex& c)noexcept{
				return FixedComplex(c.Real+FixedPoint(d),c.Imaginary);
			}
			friend FixedComplex operator+(const FixedComplex& c,const long long int d)noexcept{
				return FixedComplex(c.Real+d,c.Imaginary);
			}
			friend FixedComplex operator+(const long long int d,const FixedComplex& c)noexcept{
				return FixedComplex(c.Real+d,c.Imaginary);
			}
			FixedComplex operator+(const FixedComplex& other)const& noexcept{
				return FixedComplex(Real+other.Real,Imaginary+other.Imaginary);
			}


			friend FixedComplex operator-(const FixedComplex& c,const double d)noexcept{
				return FixedComplex(c.Real-FixedPoint(d),c.Imaginary);
			}
			friend FixedComplex operator-(const double d,const FixedComplex& c)noexcept{
				return FixedComplex(FixedPoint(d)-c.Real,c.Imaginary);
			}
			friend FixedComplex operator-(const FixedComplex& c,const long long int d)noexcept{
				return FixedComplex(c.Real-d,c.Imaginary);
			}
			friend FixedComplex operator-(const long long int d,const FixedComplex& c)noexcept{
				return FixedComplex(d-c.Real,c.Imaginary);
			}
			FixedComplex operator-(const FixedComplex& other)const& noexcept{
				return FixedComplex(Real-other.Real,Imaginary-other.Imaginary);
			}


			friend FixedComplex operator*(const FixedComplex& c,const double d)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*FixedPoint(d))>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*FixedPoint(d))>>Fixed_shift;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator*(const double d,const FixedComplex& c)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*FixedPoint(d))>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*FixedPoint(d))>>Fixed_shift;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator*(const FixedComplex& c,const long long d)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*d)>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*d)>>Fixed_shift;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator*(const long long d,const FixedComplex& c)noexcept{
				__int128_t Real_result=((__int128_t)c.Real*d)>>Fixed_shift,
						   Imaginary_result=((__int128_t)c.Imaginary*d)>>Fixed_shift;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			FixedComplex operator*(const FixedComplex& other)const& noexcept{
				__int128_t Real_result=((__int128_t)Real*other.Real-(__int128_t)Imaginary*other.Imaginary)>>Fixed_shift,
						   Imaginary_result=((__int128_t)Real*other.Imaginary+(__int128_t)Imaginary*other.Real)>>Fixed_shift;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}


			friend FixedComplex operator/(const FixedComplex& c,const double d)noexcept{
				__int128_t Real_result=(((__int128_t)c.Real<<Fixed_shift)/FixedPoint(d)),
						   Imaginary_result=(((__int128_t)c.Imaginary<<Fixed_shift)/FixedPoint(d));
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator/(const double& d,const FixedComplex& c)noexcept{
				__int128_t denominator=((__int128_t)c.Real*c.Real+c.Imaginary*c.Imaginary);
				__int128_t Real_result=FixedPoint(d)>c.Real?
									   ((__int128_t)FixedPoint(d)<<Fixed_shift)/denominator*c.Real:
									   ((__int128_t)c.Real<<Fixed_shift)/denominator*FixedPoint(d),

						   Imaginary_result=FixedPoint(d)>c.Imaginary?
								            -((__int128_t)FixedPoint(d)<<Fixed_shift)/denominator*c.Imaginary:
											-((__int128_t)c.Imaginary<<Fixed_shift)/denominator*FixedPoint(d);

				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator/(const FixedComplex& c,const long long d)noexcept{
				__int128_t Real_result=(((__int128_t)c.Real<<Fixed_shift)/d),
						   Imaginary_result=(((__int128_t)c.Imaginary<<Fixed_shift)/d);
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			friend FixedComplex operator/(const long long& d,const FixedComplex& c)noexcept{
				__int128_t denominator=((__int128_t)c.Real*c.Real+c.Imaginary*c.Imaginary);
				__int128_t Real_result=d>c.Real?
									   ((__int128_t)d<<Fixed_shift)/denominator*c.Real:
									   ((__int128_t)c.Real<<Fixed_shift)/denominator*d,

						   Imaginary_result=d>c.Imaginary?
								            -((__int128_t)d<<Fixed_shift)/denominator*c.Imaginary:
											-((__int128_t)c.Imaginary<<Fixed_shift)/denominator*d;
				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}
			FixedComplex operator/(const FixedComplex& other)const& noexcept{
				__int128_t denominator=((__int128_t)other.Real*other.Real+(__int128_t)other.Imaginary*other.Imaginary);
				__int128_t Real_result=(Real>other.Real?
									   (((__int128_t)Real<<Fixed_shift)/denominator*other.Real):
									   (((__int128_t)other.Real<<Fixed_shift)/denominator*Real))
									  +(Imaginary>other.Imaginary?
									   (((__int128_t)Imaginary<<Fixed_shift)/denominator*other.Imaginary):
									   (((__int128_t)other.Imaginary<<Fixed_shift)/denominator*Imaginary)),
						   Imaginary_result=(Imaginary>other.Real?
								   	   	    (((__int128_t)Imaginary<<Fixed_shift)/denominator*other.Real):
										    (((__int128_t)other.Real<<Fixed_shift)/denominator*Imaginary))
										   -(Real>other.Imaginary?
											(((__int128_t)Real<<Fixed_shift)/denominator*other.Imaginary):
											(((__int128_t)other.Imaginary<<Fixed_shift)/denominator*Real));

				return FixedComplex(static_cast<long long>(Real_result),static_cast<long long>(Imaginary_result));
			}

	};

	struct Complex{

		public:

		double Real=0.0,Imaginary=0.0;

		Complex():Real(0),Imaginary(0){};
		explicit Complex(double A):Real(A),Imaginary(0){};
		explicit Complex(double A,double B):Real(A),Imaginary(B){};

		friend Complex operator+(Complex& c,const double d)noexcept{
			return Complex(c.Real+d,c.Imaginary);
		}
		friend Complex operator+(const double d,Complex& c)noexcept{
			return Complex(c.Real+d,c.Imaginary);
		}
		Complex operator+(const Complex& other)const& noexcept{
			return Complex(Real+other.Real,Imaginary+other.Imaginary);
		}


		friend Complex operator-(const Complex& c,const double d)noexcept{
			return Complex(c.Real-d,c.Imaginary);
		}
		friend Complex operator-(const double d,const Complex& c)noexcept{
			return Complex(d-c.Real,c.Imaginary);
		}
		Complex operator-(const Complex& other)const& noexcept{
			return Complex(Real-other.Real,Imaginary-other.Imaginary);
		}


		friend Complex operator*(const Complex& c,const double d)noexcept{
			return Complex(c.Real*d,c.Imaginary*d);
		}
		friend Complex operator*(const double d,const Complex& c)noexcept{
			return Complex(c.Real*d,c.Imaginary*d);
		}
		Complex operator*(const Complex& other)const& noexcept{
			return Complex(Real*other.Real-Imaginary*other.Imaginary,Real*other.Imaginary+Imaginary*other.Real);
		}


		friend Complex operator/(const Complex& c,const double d)noexcept{
			if(Qubit_Simulation::QMath::AbsoluteValue(d))
				return Complex();
			return Complex(c.Real/d,c.Imaginary/d);
		}
		friend Complex operator/(const double& d,const Complex& c)noexcept{
			if(Qubit_Simulation::QMath::IsZero(c))
				return Complex();
			return Complex((d*c.Real)/(c.Real*c.Real+c.Imaginary*c.Imaginary),
						   (-d*c.Imaginary)/(c.Real*c.Real+c.Imaginary*c.Imaginary));
		}
		Complex operator/(const Complex& other)const& noexcept{
			if(Qubit_Simulation::QMath::IsZero(other))
				return Complex();
			return Complex((Real*other.Real+Imaginary*other.Imaginary)/(other.Real*other.Real+other.Imaginary*other.Imaginary),
						   (Imaginary*other.Real-Real*other.Imaginary)/(other.Real*other.Real+other.Imaginary*other.Imaginary));
		}

		inline void InputComplex(){

			std::cin>>Real>>Imaginary;

		}
		static inline void OutputComplex(const Complex& C){

			std::cout<<C.Real<<','<<C.Imaginary;

		}

	};

	struct Matrix{

	public:

		std::vector<std::vector<Complex>> data={{}};

		Matrix():data({{}}){};
		Matrix(const std::vector<double> &A){

			std::vector<std::vector<Complex>> Result(1,std::vector<Complex>(A.size()));

			for(int y=0;y<int(A.size());y++)
				Result[0][y]=Complex(A[y]);

			data=Result;

		};
		Matrix(const std::vector<Complex> &A):data({A}){};
		Matrix(const std::vector<std::vector<double>> &A){

			if(A.size()==0)
				data={{}};

			std::vector<std::vector<Complex>> Result(A.size(),std::vector<Complex>(A[0].size()));

			for(int x=0;x<int(A.size());x++)
				for(int y=0;y<int(A[x].size());y++)
					Result[x][y]=Complex(A[x][y]);

			data=Result;

		};
		Matrix(const std::vector<std::vector<Complex>> &A):data(A){};

		std::vector<Complex>& operator[](const int x){return data[x];};
		const std::vector<Complex>& operator[](const int x)const{return data[x];};

		bool IsMatrix() const{

			if(data.size()==0)
				return 0;

			int Column=data[0].size();
			for(int x=0;x<int(data.size());x++)
				if(int(data[x].size())!=Column)
					return 0;

			return 1;

		}
		bool IsSquareMatrix() const{

			return (IsMatrix()&&data.size()==data[0].size());

		}

		Matrix operator+(const Matrix& other)const& noexcept{

			if(!IsMatrix()||!other.IsMatrix())
				return {};

			if((data.size()!=other.data.size())||(data.size()==0))
				return {};

			if(data[0].size()!=other.data[0].size())
					return {};

			std::vector<std::vector<Complex>> Result(data.size(),std::vector<Complex>(data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=data[x][y]+other.data[x][y];

			return Matrix(Result);

		}
		friend Matrix operator+(const Matrix& A,const std::vector<std::vector<double>>& B){
			return A+Matrix(B);
		}
		friend Matrix operator+(const std::vector<std::vector<double>>& B,const Matrix& A){
			return Matrix(B)+A;
		}

		Matrix operator-(const Matrix& other)const& noexcept{

			if(!IsMatrix()||!other.IsMatrix())
				return {};

			if((data.size()!=other.data.size())||(data.size()==0))
				return {};

			if(data[0].size()!=other.data[0].size())
					return {};

			std::vector<std::vector<Complex>> Result(data.size(),std::vector<Complex>(data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=data[x][y]-other.data[x][y];

			return Matrix(Result);

		}
		friend Matrix operator-(const Matrix& A,const std::vector<std::vector<double>>& B){
			return A-Matrix(B);
		}
		friend Matrix operator-(const std::vector<std::vector<double>>& B,const Matrix& A){
			return Matrix(B)-A;
		}


		friend Matrix operator*(const Matrix& A,const double Scalar)noexcept{

			if(!A.IsMatrix())
				return {};

			std::vector<std::vector<Complex>> Result(A.data.size(),std::vector<Complex>(A.data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=A.data[x][y]*Scalar;

			return Matrix(Result);

		}
		friend Matrix operator*(const double Scalar,const Matrix& A)noexcept{

			return A*Scalar;

		}
		friend Matrix operator*(const Matrix& A,const Complex& Scalar)noexcept{

			if(!A.IsMatrix())
				return {};

			std::vector<std::vector<Complex>> Result(A.data.size(),std::vector<Complex>(A.data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=A.data[x][y]*Scalar;

			return Matrix(Result);

		}
		friend Matrix operator*(const Complex& Scalar,const Matrix& A)noexcept{

			return A*Scalar;

		}
		Matrix operator*(const Matrix& other)const& noexcept{

			if(!IsMatrix()||!other.IsMatrix())
				return {};

			if(data[0].size()!=other.data.size())
				return {};

			std::vector<std::vector<Complex>> Result(data.size(),std::vector<Complex>(other.data[0].size()));
			for(int Leftx=0;Leftx<int(data.size());Leftx++)
				for(int Righty=0;Righty<int(other.data[0].size());Righty++){

					Complex sum=Complex(0);
					for(int count=0;count<int(data[0].size());count++)
						sum=sum+(data[Leftx][count]*other[count][Righty]);

					Result[Leftx][Righty]=sum;

				}

			return Matrix(Result);

		}
		friend Matrix operator*(const Matrix& A,const std::vector<std::vector<double>> B)noexcept{

			return A*Matrix(B);

		}
		friend Matrix operator*(const std::vector<std::vector<double>> B,const Matrix& A)noexcept{

			return Matrix(B)*A;

		}

		friend Matrix operator/(const Matrix& A,const double Scalar)noexcept{

			if(!A.IsMatrix())
				return {};

			std::vector<std::vector<Complex>> Result(A.data.size(),std::vector<Complex>(A.data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=A.data[x][y]/Scalar;

			return Matrix(Result);

		}
		friend Matrix operator/(const Matrix& A,const Complex& Scalar)noexcept{

			if(!A.IsMatrix())
				return {};

			std::vector<std::vector<Complex>> Result(A.data.size(),std::vector<Complex>(A.data[0].size()));
			for(int x=0;x<int(Result.size());x++)
				for(int y=0;y<int(Result[x].size());y++)
					Result[x][y]=A.data[x][y]/Scalar;

			return Matrix(Result);

		}

		inline void InputMatrix(){

			int Row,Column;
			std::cout<<"The number of the Row and the Column"<<std::endl;
			std::cin>>Row>>Column;

			data=std::vector<std::vector<Complex>>(Row,std::vector<Complex>(Column,Complex()));
			for(int x=0;x<Row;x++)
				for(int y=0;y<Column;y++)
					data[x][y].InputComplex();

		}
		static inline void OutputMatrix(const Matrix& A){

			if(!A.IsMatrix())
				return ;

			for(int x=0;x<int(A.data.size());x++){
				for(int y=0;y<int(A[x].size());y++)
					Complex::OutputComplex(A[x][y]),std::cout<<' '<<' ';
				std::cout<<std::endl;
			}

			std::cout<<std::endl;

		}

	};


	class QMath{

	public:

		static inline bool IsZero(const Complex &C){

			return ((AbsoluteValue(C.Real)<Double_Epsilon)&&(AbsoluteValue(C.Imaginary)<Double_Epsilon));

		}

		static inline double AbsoluteValue(const double &D){

			return (D>0?D:-D);

		}

		static inline long long int AbsoluteValue(const FixedComplex &C){

			__int128_t Result=((__int128_t)C.Real*C.Real+(__int128_t)C.Imaginary*C.Imaginary);

			return NewtonSqrt(Result);

		}

		static inline double AbsoluteValue(const Complex &C){
			return sqrt((C.Real*C.Real)+(C.Imaginary*C.Imaginary));
		}

		static inline long long int NewtonSqrt(const __int128_t input){

			if(input==0)
				return 0;

			__int128_t Result=input;
			__int128_t x=(Result+1)/2;

			while(x<Result){

				Result=x;
				x=(Result+input/Result)/2;

			}

			return static_cast<long long int>(Result);

		}

		static inline Complex Conjugate(const Complex &C){

			return Complex(C.Real,-C.Imaginary);

		}


		static inline bool IsSameMatrix(const Matrix &Matrix1,const Matrix &Matrix2){


			if(Matrix1.data.size()!=Matrix2.data.size()||!Matrix1.IsMatrix()||!Matrix2.IsMatrix())
				return 0;

			if(Matrix1[0].size()!=Matrix2[0].size())
				return 0;

			for(int x=0;x<int(Matrix1.data.size());x++)
				for(int y=0;y<int(Matrix1[x].size());y++)
					if(!IsZero(Matrix1[x][y]-Matrix2[x][y]))
						return 0;

			return 1;

		}

		static inline bool IsDiagonalMatrix(const Matrix &A){

			if(!A.IsSquareMatrix())
				return 0;

			for(int x=0;x<int(A.data.size());x++)
				for(int y=0;y<int(A[0].size());y++)
					if(x!=y&&!IsZero(A[x][y]))
						return 0;

			return 1;

		}

		static inline double NormOfVector(const std::vector<Complex> &Vector){

			double Sum=0;
			for(int count=0;count<int(Vector.size());count++)
				Sum+=(Vector[count].Real*Vector[count].Real)+(Vector[count].Imaginary*Vector[count].Imaginary);

			return sqrt(Sum);

		}

		static inline std::vector<Complex> ScalarMultiplicationVector(const Complex &Scalar,const std::vector<Complex> &Vector){

			std::vector<Complex> Result(Vector.size());
			for(int count=0;count<int(Vector.size());count++)
				Result[count]=Vector[count]*Scalar;

			return Result;

		}

		static inline Matrix GenerateReflectionMatrix(const std::vector<Complex> &NormalVector){

			return GenerateUnitMatrix(NormalVector.size())-((2/DotProduct(NormalVector,NormalVector))*OuterProduct(NormalVector,NormalVector));

		}

		static inline Complex DotProduct(const std::vector<Complex> &Vector1,const std::vector<Complex> &Vector2){

			if(Vector1.size()!=Vector2.size())
				return Complex();

			Complex Sum=Complex();
			for(int count=0;count<int(Vector1.size());count++)
				Sum=Sum+Conjugate(Vector1[count])*Vector2[count];

			return Sum;

		}

		static inline Matrix OuterProduct(const std::vector<Complex> &Vector1,const std::vector<Complex> &Vector2){

			if(Vector1.size()!=Vector2.size())
				return {};

			return Matrix(Vector1)*TransposeMatrix(Matrix(Vector2));

		}

		static Matrix GenerateUnitMatrix(const int &order){

			Matrix UnitMatrix=Matrix(std::vector<std::vector<Complex>>(order,std::vector<Complex>(order,Complex())));
			for(int count=0;count<order;count++)
				UnitMatrix[count][count]=Complex(1);

			return UnitMatrix;

		}

		static inline void AdditionOneRowToAnotherRow(const Complex &Scalar,const int &AddedRow,const int &Row,Matrix &A){

			if(!A.IsMatrix())
				return ;

			if(AddedRow>=int(A.data.size())||AddedRow<0||Row>=int(A.data.size())||Row<0||Row==AddedRow)
				return ;

			for(int y=0;y<int(A[0].size());y++)
				A[Row][y]=A[Row][y]+A[AddedRow][y]*Scalar;

		}

		static inline void MultipleOfARow(const Complex &Scalar,const int &Row,Matrix &A){

			if(!A.IsMatrix())
				return;

			if(Row>=int(A.data.size())||Row<0||(Scalar.Real==0&&Scalar.Imaginary==0))
				return;

			for(int y=0;y<int(A[Row].size());y++)
				A[Row][y]=A[Row][y]*Scalar;

		}

		static inline void InterchangeOfTwoRow(const std::pair<int,int> &InterchangeRow,Matrix &A){

			if(!A.IsMatrix())
				return;

			if(InterchangeRow.first>=int(A.data.size())||InterchangeRow.first<0||InterchangeRow.second>=int(A.data.size())||InterchangeRow.second<0)
				return;

			std::swap(A[InterchangeRow.first],A[InterchangeRow.second]);

		}

		static inline Matrix HermitianTransposeMatrix(Matrix A){

			if(!A.IsSquareMatrix())
				return Matrix();

			for(int x=0;x<int(A.data.size());x++){

				A[x][x]=Conjugate(A[x][x]);

				for(int y=0;y<x;y++)
					std::swap(A[x][y],A[y][x]),
					A[x][y]=Conjugate(A[x][y]),
					A[y][x]=Conjugate(A[y][x]);

			}

			return A;

		}

		static inline Matrix TransposeMatrix(Matrix A){

			if(!A.IsSquareMatrix())
				return Matrix();

			for(int x=0;x<int(A.data.size());x++)
				for(int y=0;y<x;y++)
					std::swap(A[x][y],A[y][x]);

			return A;

		}

		static Complex DeterminantsOfMatrix(const Matrix &A){

			if(!A.IsSquareMatrix())
				return Complex();

			if(A.data.size()==1)
				return A[0][0];

			if(A.data.size()==2)
				return A[0][0]*A[1][1]-A[0][1]*A[1][0];

			Complex Sum=Complex();

			for(int y=0;y<int(A[0].size());y++){

				Matrix ReducedMatrix=Matrix(std::vector<std::vector<Complex>>(A.data.size()-1));
				Complex now=A[0][y];

				for(int x=1;x<int(A.data.size());x++)
					for(int Column=0;Column<int(A[0].size());Column++)
						if(y!=Column)
							ReducedMatrix[x-1].push_back(A[x][Column]);

				if(y%2==1)
					now=now-1;

				Sum=Sum+now*DeterminantsOfMatrix(ReducedMatrix);

			}

			return Sum;

		}

		static Matrix GenerateInverseMatrix(const Matrix &A){

			if(!A.IsSquareMatrix())
				return {};

			if(IsZero(DeterminantsOfMatrix(A))||A[0].size()==0)
				return {};

			Matrix Now=A,Inverse=GenerateUnitMatrix(A.data.size());
			for(int y=0;y<int(Now.data.size());y++){

				if(IsZero(Now[y][y])){

					for(int x=0;x<int(Now.data.size());x++)
						if(!IsZero(Now[x][y])){
							InterchangeOfTwoRow({x,y},Now);
							InterchangeOfTwoRow({x,y},Inverse);
							break;
						}

				}

				MultipleOfARow(1/Now[y][y],y,Inverse);
				MultipleOfARow(1/Now[y][y],y,Now);

				for(int x=0;x<int(Now.data.size());x++)
					if(x!=y&&!IsZero(Now[x][y])){
						AdditionOneRowToAnotherRow(Complex()-Now[x][y],y,x,Inverse);
						AdditionOneRowToAnotherRow(Complex()-Now[x][y],y,x,Now);
					}

			}

			return Inverse;

		}

		static Matrix GenerateQuantumLogicGate(const Matrix &A){

		}

		static Matrix GenerateDiagonalMatrix(const Matrix &A){

			if(!A.IsSquareMatrix()||A.data.size()==0)
				return {};

			if(A.data.size()==1)
				return A;

			Matrix Now=A,Front=Matrix();
			std::pair<Matrix,Matrix> QR;
			while(!IsDiagonalMatrix(Now)&&!IsSameMatrix(Front,Now)){

				Front=Now;
				QR=GenerateQRFactorization(Now);
				Now=QR.second*QR.first;

			}

			for(int x=0;x<int(Now.data.size());x++)
				for(int y=0;y<int(Now[x].size());y++)
					if(x!=y)
						Now[x][y]=Complex();

			return Now;

		}

		static std::pair<Matrix,Matrix> GenerateQRFactorization(const Matrix &A){

			if(!A.IsSquareMatrix()||A.data.size()==0||A.data.size()==1)
				return {{},{}};


			Matrix Q=GenerateUnitMatrix(A.data.size()),R=A;
			for(int y=0;y<int(A.data.size()-1);y++){

				std::vector<Complex> Set(A.data.size()-y);
				for(int x=y;x<int(A.data.size());x++)
					Set[x-y]=R[x][y];


				Set[0]=Set[0]-NormOfVector(Set);
				Set=ScalarMultiplicationVector(Complex(1/NormOfVector(Set)),Set);
				Matrix Caculate=GenerateReflectionMatrix(Set);
				Caculate=GenerateBlockEmbeddingForHouseholder(A.data.size(),Caculate);

				Q=Q*Caculate;
				R=Caculate*R;

			}

			return {Q,R};

		}

		static Matrix GenerateBlockEmbeddingForHouseholder(const int &order,const Matrix &A){

			if(!A.IsSquareMatrix()||int(A.data.size())>order)
				return {};

			if(order==int(A.data.size()))
				return A;

			Matrix Result=Matrix(std::vector<std::vector<Complex>>(order,std::vector<Complex>(order,Complex())));
			for(int count=0;count<order;count++)
				Result[count][count]=Complex(1);

			for(int x=order-A.data.size();x<order;x++)
				for(int y=order-A.data.size();y<order;y++)
					Result[x][y]=A[x-order+A.data.size()][y-order+A.data.size()];

			return Result;

		}





	};

	private:

	int Qubit_Amount;
	std::vector<std::pair<bool,int>> Qubit_Set_Observation;
	std::unordered_map<int,std::pair<std::vector<int>,std::vector<FixedComplex>>> Entangled_Qubit_Set_Pointer;
	std::unordered_map<int,std::pair<std::vector<int>,std::vector<FixedComplex>>*> Entangled_Qubit_Set;

	void CombineEntangledQubitSet(int situation1,int situation2){

		if(Entangled_Qubit_Set.find(situation1)==Entangled_Qubit_Set.end()||Entangled_Qubit_Set.find(situation2)==Entangled_Qubit_Set.end())
			return ;

		std::vector<int> First_Qubit_Name_Set=Entangled_Qubit_Set[situation1]->first,
						 Second_Qubit_Name_Set=Entangled_Qubit_Set[situation2]->first;

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			First_Qubit_Name_Set.push_back(Second_Qubit_Name_Set[count]);

		Entangled_Qubit_Set[situation1]->first=First_Qubit_Name_Set;

		std::vector<FixedComplex> First_Qubit_Set=Entangled_Qubit_Set[situation1]->second,
							 Second_Qubit_Set=Entangled_Qubit_Set[situation2]->second,
							 Result(First_Qubit_Set.size()*Second_Qubit_Set.size(),FixedComplex());

		for(unsigned long long x=0;x<First_Qubit_Set.size();x++)
			for(unsigned long long y=0;y<Second_Qubit_Set.size();y++)
				if((First_Qubit_Set[x].Real!=0||First_Qubit_Set[x].Imaginary!=0)&&(Second_Qubit_Set[y].Real!=0||Second_Qubit_Set[y].Imaginary!=0))
					Result[x*Second_Qubit_Set.size()+y]=First_Qubit_Set[x]*Second_Qubit_Set[y];

		Entangled_Qubit_Set[situation1]->second=Result;

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			if(Second_Qubit_Name_Set[count]!=situation1)
				Entangled_Qubit_Set_Pointer.erase(Second_Qubit_Name_Set[count]);

		for(unsigned long long count=0;count<Second_Qubit_Name_Set.size();count++)
			Entangled_Qubit_Set[Second_Qubit_Name_Set[count]]=Entangled_Qubit_Set[situation1];

	}

	constexpr inline unsigned long long BitRemova(unsigned long long num,const unsigned long long situation){

		return (num&((1ULL<<situation)-1ULL))|((num>>(situation+1ULL))<<situation);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num,unsigned long long situation1,unsigned long long situation2){

		if(situation1<situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;
		}

		num=(num&((1ULL<<situation1)-1ULL))|((num>>(situation1+1ULL))<<situation1);

		return (num&((1ULL<<situation2)-1ULL))|((num>>(situation2+1ULL))<<situation2);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num,unsigned long long situation1,unsigned long long situation2,unsigned long long situation3){

		if(situation1<situation2){
			situation1=situation1^situation2;
			situation2=situation1^situation2;
			situation1=situation1^situation2;
		}
		if(situation1<situation3){
			situation1=situation1^situation3;
			situation3=situation1^situation3;
			situation1=situation1^situation3;
		}

		num=(num&((1ULL<<situation1)-1ULL))|((num>>(situation1+1ULL))<<situation1);

		if(situation2<situation3){
			situation2=situation2^situation3;
			situation3=situation2^situation3;
			situation2=situation2^situation3;
		}
		num=(num&((1ULL<<situation2)-1ULL))|((num>>(situation2+1ULL))<<situation2);

		return (num&((1ULL<<situation3)-1ULL))|((num>>(situation3+1ULL))<<situation3);

	}

	constexpr inline unsigned long long BitAdd(unsigned long long num,const int situation,bool add){

		return ((num>>situation)<<(situation+1ULL))|(add<<situation)|(num&((1ULL<<situation)-1ULL));

	}
	inline unsigned long long BitAdd(unsigned long long num,std::vector<std::pair<int,bool>> Control){

		std::sort(Control.begin(),Control.end());

		for(int count=0;count<int(Control.size());count++)
			num=BitAdd(num,Control[count].first,Control[count].second);

		return num;

	}

	int GetSituation(const int situation){

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=-1;
		for(unsigned long long count=0;count<Qubit_Set.size();count++)
			if(Qubit_Set[count]==situation)
				Qubit_Situation=count;

		if(Qubit_Situation==-1)
			return -1;

		return Qubit_Set.size()-Qubit_Situation-1;

	}

	bool IsQubitUnoperable(const int situation){

		return (Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end()||IsObservered(situation));

	}

	public:

	Qubit_Simulation():Qubit_Amount(0){Entangled_Qubit_Set.clear();Qubit_Set_Observation.clear();Entangled_Qubit_Set_Pointer.clear();};
	Qubit_Simulation(int num){

		Qubit_Amount=num;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

		for(int count=0;count<num;count++)
			Entangled_Qubit_Set_Pointer[count]={{count},{FixedComplex(Fixed_Point),FixedComplex()}},
			Qubit_Set_Observation.push_back({0,-1}),
			Entangled_Qubit_Set[count]=&Entangled_Qubit_Set_Pointer[count];

	};

	void ResetQubitSet(){

		Qubit_Amount=0;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

	}

	void OuputEntangledQubitSet(const int situation){

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return;

		if(Entangled_Qubit_Set[situation]==nullptr)
			return;

		std::vector<FixedComplex> Output_Qubit=Entangled_Qubit_Set[situation]->second;
		for(unsigned long long count=0;count<Output_Qubit.size();count++)
			std::cout<<FixedPointToDouble(Output_Qubit[count].Real)<<' '<<FixedPointToDouble(Output_Qubit[count].Imaginary)<<std::endl;
			//std::cout<<Output_Qubit[count].Real<<' '<<Output_Qubit[count].Imaginary<<std::endl;

	}

	bool IsObservered(const int situation){

		return Qubit_Set_Observation[situation].first;

	}

	void CutQubitSet(const int situation,long long int Probability0,long long int Probability1){

		if(Probability0==0)
			Qubit_Set_Observation[situation]={1,1},Probability1=Fixed_Point;
		else if(Probability0==1)
			Qubit_Set_Observation[situation]={1,0},Probability0=Fixed_Point;
		else{

			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::uniform_int_distribution<long long int> dist(0,Fixed_Point);
			long long int RandNum=dist(gen);

			if(RandNum<=Probability0)
				Qubit_Set_Observation[situation]={1,0};
			else
				Qubit_Set_Observation[situation]={1,1};

		}

		long long int State=Qubit_Set_Observation[situation].second;
		long long int Probability=(State==0)?Probability0:Probability1;

		std::vector<FixedComplex> Qubit_Set_Value=Entangled_Qubit_Set[situation]->second,
							 New_Qubit_Set_Value((Qubit_Set_Value.size()>>1)-(Qubit_Set_Value.size()==2),FixedComplex());

		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=GetSituation(situation);

		for(unsigned long long int count=0;count<(Qubit_Set_Value.size()>>1)-(Qubit_Set_Value.size()==2);count++)
			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)].Imaginary!=0)
				New_Qubit_Set_Value[count]=Qubit_Set_Value[BitAdd(count,Qubit_Situation,State)]/QMath::NewtonSqrt((__int128_t)Probability<<Fixed_shift);

		Entangled_Qubit_Set[situation]->second=New_Qubit_Set_Value;

		Qubit_Set.erase(std::remove(Qubit_Set.begin(),Qubit_Set.end(),situation),Qubit_Set.end());
		Entangled_Qubit_Set[situation]->first=Qubit_Set;

		Entangled_Qubit_Set[situation]=nullptr;

	}

	int ObserverQubit(const int situation){

		if(IsObservered(situation))
			return Qubit_Set_Observation[situation].second;

		if(Entangled_Qubit_Set.find(situation)==Entangled_Qubit_Set.end())
			return -1;

		std::pair<std::vector<int>,std::vector<FixedComplex>>* Qubit_Set_Pointer=Entangled_Qubit_Set[situation];
		std::vector<int> Qubit_Set=Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation=GetSituation(situation);

		std::vector<FixedComplex> Qubit_Set_Value=Entangled_Qubit_Set[situation]->second;
		long long int Probability0=0,Probability1=0;

		for(unsigned long long int count=0;count<Qubit_Set_Value.size()>>1;count++){

			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)].Imaginary!=0){
				__int128_t Caculate=static_cast<__int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Situation,0ULL)]));
				Caculate*=Caculate;
				Caculate>>=Fixed_shift;
				Probability0+=Caculate;
			}

			if(Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)].Imaginary!=0){
				__int128_t Caculate=static_cast<__int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Situation,1ULL)]));
				Caculate*=Caculate;
				Caculate>>=Fixed_shift;
				Probability1+=Caculate;
			}

		}

		CutQubitSet(situation,Probability0,Probability1);
		Qubit_Set=Qubit_Set_Pointer->first;
		Qubit_Set_Value=Qubit_Set_Pointer->second;

		for(long long int setting=Qubit_Set.size()-1;setting>=0;setting--){

			Qubit_Set=Qubit_Set_Pointer->first;
			Qubit_Set_Value=Qubit_Set_Pointer->second;

			Probability0=0,Probability1=0;

			for(unsigned long long int count=0;count<Qubit_Set_Value.size()>>1;count++){

				if(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)].Imaginary!=0){
					__int128_t Caculate=static_cast<__int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,0ULL)]));
					Caculate*=Caculate;
					Caculate>>=Fixed_shift;
					Probability0+=Caculate;
				}

				if(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)].Real!=0||Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)].Imaginary!=0){
					__int128_t Caculate=static_cast<__int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count,Qubit_Set.size()-setting-1,1ULL)]));
					Caculate*=Caculate;
					Caculate>>=Fixed_shift;
					Probability1+=Caculate;
				}

			}

			if(Probability0==0||Probability1==0)
				CutQubitSet(Qubit_Set[setting],Probability0,Probability1);

		}

		return Qubit_Set_Observation[situation].second;

	}

	int GenerateQubit(){

		Qubit_Amount++;

		Entangled_Qubit_Set_Pointer[Qubit_Amount-1]={{Qubit_Amount-1},{FixedComplex(Fixed_Point),FixedComplex()}};
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Entangled_Qubit_Set.size()];
		Qubit_Set_Observation.push_back({0,-1});

		return Qubit_Amount-1;

	}

	int GeneratePositiveBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{FixedComplex(Root_half),FixedComplex(),FixedComplex(),FixedComplex(Root_half)}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GenerateNegativeBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{FixedComplex(Root_half),FixedComplex(),FixedComplex(),FixedComplex(-Root_half)}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GeneratePositiveOrthogonalBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{FixedComplex(),FixedComplex(Root_half),FixedComplex(Root_half),FixedComplex()}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}

	int GenerateSingletBellState(){

		Qubit_Amount+=2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount-2]={{Qubit_Amount-2,Qubit_Amount-1},{FixedComplex(),FixedComplex(Root_half),FixedComplex(-Root_half),FixedComplex()}};
		Entangled_Qubit_Set[Qubit_Amount-2]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];
		Entangled_Qubit_Set[Qubit_Amount-1]=&Entangled_Qubit_Set_Pointer[Qubit_Amount-2];

		Qubit_Set_Observation.push_back({0,-1});
		Qubit_Set_Observation.push_back({0,-1});
		return Qubit_Amount-2;

	}


	//單量子位元邏輯閘


	void CostomGate(const int qubit,std::vector<std::vector<Complex>> Matrix){





	}

	void HadamardGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++){

			FixedComplex First=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)],
					 	 Second=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)];

			 Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]=(First*(Root_half)+(Second*(Root_half)));
			 Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=(First*(Root_half))-(Second*(Root_half));

		}

	}

	void XGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			std::swap(Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)],
					  Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]);

	}

	void YGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++){

			std::swap(Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)],
					  Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]);

			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]*FixedComplex(0LL,-Fixed_Point);
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(0LL,Fixed_Point);

		}

	}

	void ZGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*(-Fixed_Point);

	}

	void SGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(0ULL,Fixed_Point);

	}

	void SdgGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(0ULL,-Fixed_Point);

	}

	void TGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(Root_half,Root_half);


	}

	void TdgGate(const int qubit){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++)
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(Root_half,-Root_half);

	}

	void RxGate(const int qubit,const double angle){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		double Cos_Angle=cos(angle*Pi/180/2),
			   Sin_Angle=sin(angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++){

			FixedComplex First=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)],
						Second=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)];

			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]=First*FixedPoint(Cos_Angle)+Second*FixedComplex(0,FixedPoint(-Sin_Angle));
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=First*FixedComplex(0,FixedPoint(-Sin_Angle))+Second*FixedPoint(Cos_Angle);

		}

	}

	void RyGate(const int qubit,const double angle){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		double Cos_Angle=cos(angle*Pi/180/2),
			   Sin_Angle=sin(angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++){

			FixedComplex First=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)],
						 Second=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)];

			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]=First*FixedPoint(Cos_Angle)+Second*FixedPoint(-Sin_Angle);
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=First*FixedPoint(-Sin_Angle)+Second*FixedPoint(Cos_Angle);

		}

	}

	void RzGate(const int qubit,const double angle){

		if(IsQubitUnoperable(qubit))
			return ;

		int Qubit_Situation=GetSituation(qubit);

		double Cos_Angle=cos(angle*Pi/180/2),
			   Sin_Angle=sin(angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit]->second.size()>>1;count++){

			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,0ULL)]*FixedComplex(FixedPoint(Cos_Angle),FixedPoint(-Sin_Angle));
			Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]=Entangled_Qubit_Set[qubit]->second[BitAdd(count,Qubit_Situation,1ULL)]*FixedComplex(FixedPoint(Cos_Angle),FixedPoint(Sin_Angle));

		}

	}

	void UGate(const int qubit, const double theta, const double phi, const double lambda){

		RzGate(qubit,theta);
		RyGate(qubit,phi);
		RzGate(qubit,lambda);

	}


	//雙量子位元邏輯閘

	void CNOTGate(const int controlQubit,const int targetQubit){

		if(IsQubitUnoperable(controlQubit)||
		   IsQubitUnoperable(targetQubit)||
		   controlQubit==targetQubit)
			return ;

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[targetQubit]->first)
			CombineEntangledQubitSet(controlQubit,targetQubit);

		int Qubit_Situation1=GetSituation(controlQubit),
			Qubit_Situation2=GetSituation(targetQubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[controlQubit]->second.size()>>2;count++)
			std::swap(Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1}})],
					  Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})]);

	}

	void CZGate(const int controlQubit,const int targetQubit){

		if(IsQubitUnoperable(controlQubit)||
		   IsQubitUnoperable(targetQubit)||
		   controlQubit==targetQubit)
			return ;

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[targetQubit]->first)
			CombineEntangledQubitSet(controlQubit,targetQubit);

		int Qubit_Situation1=GetSituation(controlQubit),
			Qubit_Situation2=GetSituation(targetQubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[controlQubit]->second.size()>>2;count++)
			Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1}})]=
			Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1}})]*FixedComplex(-Fixed_Point);

	}

	void SqrtSWAPGate(const int qubit1,const int qubit2){

		if(IsQubitUnoperable(qubit1)||
		   IsQubitUnoperable(qubit2)||
		   qubit1==qubit2)
			return ;

		if(Entangled_Qubit_Set[qubit1]->first!=Entangled_Qubit_Set[qubit2]->first)
			CombineEntangledQubitSet(qubit1,qubit2);

		int Qubit_Situation1=GetSituation(qubit1),
			Qubit_Situation2=GetSituation(qubit2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit1]->second.size()>>2;count++){

			FixedComplex First=Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})],
						 Second=Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})];

			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})]=First*FixedComplex(Fixed_Point>>1,Fixed_Point>>1)+Second*FixedComplex(Fixed_Point>>1,-(Fixed_Point>>1));
			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})]=First*FixedComplex(Fixed_Point>>1,-(Fixed_Point>>1))+Second*FixedComplex(Fixed_Point>>1,Fixed_Point>>1);

		}
	}

 	void SWAPGate(const int qubit1,const int qubit2){

		if(IsQubitUnoperable(qubit1)||
		   IsQubitUnoperable(qubit2)||
		   qubit1==qubit2)
			return ;

		if(Entangled_Qubit_Set[qubit1]->first!=Entangled_Qubit_Set[qubit2]->first)
			CombineEntangledQubitSet(qubit1,qubit2);

		int Qubit_Situation1=GetSituation(qubit1),
			Qubit_Situation2=GetSituation(qubit2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit1]->second.size()>>2;count++)
			std::swap(Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})],
					  Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})]);

	}

	void ISWAPGate(const int qubit1,const int qubit2){

		if(IsQubitUnoperable(qubit1)||
		   IsQubitUnoperable(qubit2)||
		   qubit1==qubit2)
			return ;

		if(Entangled_Qubit_Set[qubit1]->first!=Entangled_Qubit_Set[qubit2]->first)
			CombineEntangledQubitSet(qubit1,qubit2);

		int Qubit_Situation1=GetSituation(qubit1),
			Qubit_Situation2=GetSituation(qubit2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[qubit1]->second.size()>>2;count++){
			std::swap(Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})],
					  Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})]);
			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})]=
			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0}})]*FixedComplex(0,Fixed_Point);
			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})]=
			Entangled_Qubit_Set[qubit1]->second[BitAdd(count,{{Qubit_Situation1,0},{Qubit_Situation2,1}})]*FixedComplex(0,Fixed_Point);

		}
	}


	//三量子量子位元邏輯閘

	void CSWAPGate(const int controlQubit,const int qubit1,const int qubit2){

		if(IsQubitUnoperable(controlQubit)||
		   IsQubitUnoperable(qubit1)||
		   IsQubitUnoperable(qubit2)||
		   controlQubit==qubit1||
		   controlQubit==qubit2||
		   qubit1==qubit2)
			return ;

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[qubit1]->first)
			CombineEntangledQubitSet(controlQubit,qubit1);

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[qubit2]->first)
			CombineEntangledQubitSet(controlQubit,qubit2);

		int Qubit_Situation1=GetSituation(controlQubit),
			Qubit_Situation2=GetSituation(qubit1),
			Qubit_Situation3=GetSituation(qubit2);


		for(unsigned long long count=0;count<Entangled_Qubit_Set[controlQubit]->second.size()>>3;count++)
			std::swap(Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,0}})],
					  Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,0},{Qubit_Situation3,1}})]);

	}

	void CCNOTGate(const int controlQubit,const int controlQubit2,const int targetQubit){

		if(IsQubitUnoperable(controlQubit)||
		   IsQubitUnoperable(controlQubit2)||
		   IsQubitUnoperable(targetQubit)||
		   controlQubit==controlQubit2||
		   controlQubit==targetQubit||
		   controlQubit2==targetQubit)
			return ;

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[controlQubit2]->first)
			CombineEntangledQubitSet(controlQubit,controlQubit2);

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[targetQubit]->first)
			CombineEntangledQubitSet(controlQubit,targetQubit);

		int Qubit_Situation1=GetSituation(controlQubit),
			Qubit_Situation2=GetSituation(controlQubit2),
			Qubit_Situation3=GetSituation(targetQubit);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[controlQubit]->second.size()>>3;count++)
			std::swap(Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,0}})],
					  Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,1}})]);

	}

	void DeutschGate(const int controlQubit,const int controlQubit2,const int targetQubit,const double angle){

		if(IsQubitUnoperable(controlQubit)||
		   IsQubitUnoperable(controlQubit2)||
		   IsQubitUnoperable(targetQubit)||
		   controlQubit==controlQubit2||
		   controlQubit==targetQubit||
		   controlQubit2==targetQubit)
			return ;

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[controlQubit2]->first)
			CombineEntangledQubitSet(controlQubit,controlQubit2);

		if(Entangled_Qubit_Set[controlQubit]->first!=Entangled_Qubit_Set[targetQubit]->first)
			CombineEntangledQubitSet(controlQubit,targetQubit);

		int Qubit_Situation1=GetSituation(controlQubit),
			Qubit_Situation2=GetSituation(controlQubit2),
			Qubit_Situation3=GetSituation(targetQubit);

		double Cos_Angle=cos(angle*Pi/180/2),
			   Sin_Angle=sin(angle*Pi/180/2);

		for(unsigned long long count=0;count<Entangled_Qubit_Set[controlQubit]->second.size()>>3;count++){

			FixedComplex First=Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,0}})],
						 Second=Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,1}})];

			Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,0}})]=First*FixedComplex(0,Cos_Angle)+Second*Sin_Angle;
			Entangled_Qubit_Set[controlQubit]->second[BitAdd(count,{{Qubit_Situation1,1},{Qubit_Situation2,1},{Qubit_Situation3,1}})]=First*Sin_Angle+Second*FixedComplex(0,-Cos_Angle);
		}
	}


};

int main() {





	Qubit_Simulation a=Qubit_Simulation(3);

	Qubit_Simulation::Matrix A;
	A.InputMatrix();
	//Qubit_Simulation::Complex::OutputComplex(Qubit_Simulation::QMath::DeterminantsOfMatrix(A));
	Qubit_Simulation::Matrix::OutputMatrix(Qubit_Simulation::QMath::GenerateQRFactorization(A).first);



	/*
	a.HadamardGate(0);
	a.CNOTGate(0,1);
	a.OuputEntangledQubitSet(0);
	a.HadamardGate(1);

	a.OuputEntangledQubitSet(0);
	std::cout<<std::endl;
	std::cout<<a.ObserverQubit(0)<<std::endl;
	std::cout<<a.ObserverQubit(1)<<std::endl;


*/
/*
	a.PositiveGenerateBellState();

	a.Hadamard(1);
	a.CNOT(0,1);

	//a.CCNOT(0,1,2);
	//a.CSWAP(0,1,2);
	a.OuputEntangledQubitSet(0);*/


}
