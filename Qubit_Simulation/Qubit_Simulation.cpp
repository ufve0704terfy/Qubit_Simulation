#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
#include <cmath>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;

static constexpr const double Double_Epsilon = 1e-9;
static constexpr const long long int Fixed_Point = (1LL << 62);
static constexpr const int Fixed_shift = 62;
static constexpr const double Root_half = 0.70710678118;
static constexpr const double Pi = 3.141592658979;
static constexpr const unsigned long long int Fixed_Pi = 14488038916154245120ULL;
static constexpr const long long int Tolerent_Round = 1LL << 10;

/*qubit中的複數的虛數部份以及實數部分都介於-1~1
 * 因此使用long long配合定點數可以有效的解決運算中的浮點誤差問題
 * 以2^62為縮放基準
 */
struct EntangledGroupSnapshot {
	std::vector<int>                       qubit_ids;
	std::vector<std::pair<double, double>>  amplitudes;  // {real, imag}
};

enum class GateType {
	Identity,
	// 單 qubit
	H, X, Y, Z, S, Sdg, T, Tdg,
	Rx, Ry, Rz,
	// 雙 qubit
	CNOT, CZ, SWAP, iSWAP, SqrtSWAP,
	// 三 qubit
	CCNOT, CSWAP, Deutsch,
	// 測量
	Measure, Parity,
	// 控制流
	Conditional, Label, Goto
};

// ══════════════════════════════════════════════════════════
// 單一指令
// ══════════════════════════════════════════════════════════
struct CircuitInstruction {
	GateType                        type = GateType::Identity;
	std::vector<int>                qubits = {};
	double                          param = 0.0;
	int                             expected_val = 0;    // Conditional 用
	std::string                     label = "";   // Label / Goto 用
	std::vector<CircuitInstruction> sub_insts = {};   // Conditional 用
};

class FixedComplex {

public:

	static long long int FixedPoint(const double input) {

		return static_cast<long long>(input * Fixed_Point);

	}
	static double FixedPointToDouble(const long long int input) {

		double result = 1.0 * input / Fixed_Point;

		return result;

	}

	long long int Real = 0, Imaginary = 0;

	FixedComplex() :Real(0), Imaginary(0) {};
	explicit FixedComplex(double A) :Real(FixedPoint(A)), Imaginary(0) {};
	explicit FixedComplex(double A, double B) :Real(FixedPoint(A)), Imaginary(FixedPoint(B)) {};
	FixedComplex(long long int A) :Real(A), Imaginary(0) {};
	FixedComplex(long long int A, long long int B) :Real(A), Imaginary(B) {};

	friend FixedComplex operator+(FixedComplex& c, const double d)noexcept {
		return FixedComplex(c.Real + FixedPoint(d), c.Imaginary);
	}
	friend FixedComplex operator+(const double d, FixedComplex& c)noexcept {
		return FixedComplex(c.Real + FixedPoint(d), c.Imaginary);
	}
	friend FixedComplex operator+(const FixedComplex& c, const long long int d)noexcept {
		return FixedComplex(c.Real + d, c.Imaginary);
	}
	friend FixedComplex operator+(const long long int d, const FixedComplex& c)noexcept {
		return FixedComplex(c.Real + d, c.Imaginary);
	}
	FixedComplex operator+(const FixedComplex& other)const& noexcept {
		return FixedComplex(Real + other.Real, Imaginary + other.Imaginary);
	}


	friend FixedComplex operator-(const FixedComplex& c, const double d)noexcept {
		return FixedComplex(c.Real - FixedPoint(d), c.Imaginary);
	}
	friend FixedComplex operator-(const double d, const FixedComplex& c)noexcept {
		return FixedComplex(FixedPoint(d) - c.Real, c.Imaginary);
	}
	friend FixedComplex operator-(const FixedComplex& c, const long long int d)noexcept {
		return FixedComplex(c.Real - d, c.Imaginary);
	}
	friend FixedComplex operator-(const long long int d, const FixedComplex& c)noexcept {
		return FixedComplex(d - c.Real, c.Imaginary);
	}
	FixedComplex operator-(const FixedComplex& other)const& noexcept {
		return FixedComplex(Real - other.Real, Imaginary - other.Imaginary);
	}


	friend FixedComplex operator*(const FixedComplex& c, const double d)noexcept {
		int128_t Real_result = ((int128_t)c.Real * FixedPoint(d)) >> Fixed_shift,
			Imaginary_result = ((int128_t)c.Imaginary * FixedPoint(d)) >> Fixed_shift;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator*(const double d, const FixedComplex& c)noexcept {
		int128_t Real_result = ((int128_t)c.Real * FixedPoint(d)) >> Fixed_shift,
			Imaginary_result = ((int128_t)c.Imaginary * FixedPoint(d)) >> Fixed_shift;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator*(const FixedComplex& c, const long long d)noexcept {
		int128_t Real_result = ((int128_t)c.Real * d) >> Fixed_shift,
			Imaginary_result = ((int128_t)c.Imaginary * d) >> Fixed_shift;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator*(const long long d, const FixedComplex& c)noexcept {
		int128_t Real_result = ((int128_t)c.Real * d) >> Fixed_shift,
			Imaginary_result = ((int128_t)c.Imaginary * d) >> Fixed_shift;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	FixedComplex operator*(const FixedComplex& other)const& noexcept {
		int128_t Real_result = ((int128_t)Real * other.Real - (int128_t)Imaginary * other.Imaginary) >> Fixed_shift,
			Imaginary_result = ((int128_t)Real * other.Imaginary + (int128_t)Imaginary * other.Real) >> Fixed_shift;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}


	friend FixedComplex operator/(const FixedComplex& c, const double d)noexcept {
		int128_t Real_result = (((int128_t)c.Real << Fixed_shift) / FixedPoint(d)),
			Imaginary_result = (((int128_t)c.Imaginary << Fixed_shift) / FixedPoint(d));
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator/(const double& d, const FixedComplex& c)noexcept {
		int128_t denominator = ((int128_t)c.Real * c.Real + c.Imaginary * c.Imaginary);
		int128_t Real_result = FixedPoint(d) > c.Real ?
			((int128_t)FixedPoint(d) << Fixed_shift) / denominator * c.Real :
			((int128_t)c.Real << Fixed_shift) / denominator * FixedPoint(d),

			Imaginary_result = FixedPoint(d) > c.Imaginary ?
			-((int128_t)FixedPoint(d) << Fixed_shift) / denominator * c.Imaginary :
			-((int128_t)c.Imaginary << Fixed_shift) / denominator * FixedPoint(d);

		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator/(const FixedComplex& c, const long long d)noexcept {
		int128_t Real_result = (((int128_t)c.Real << Fixed_shift) / d),
			Imaginary_result = (((int128_t)c.Imaginary << Fixed_shift) / d);
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	friend FixedComplex operator/(const long long& d, const FixedComplex& c)noexcept {
		int128_t denominator = ((int128_t)c.Real * c.Real + c.Imaginary * c.Imaginary);
		int128_t Real_result = d > c.Real ?
			((int128_t)d << Fixed_shift) / denominator * c.Real :
			((int128_t)c.Real << Fixed_shift) / denominator * d,

			Imaginary_result = d > c.Imaginary ?
			-((int128_t)d << Fixed_shift) / denominator * c.Imaginary :
			-((int128_t)c.Imaginary << Fixed_shift) / denominator * d;
		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}
	FixedComplex operator/(const FixedComplex& other)const& noexcept {
		int128_t denominator = ((int128_t)other.Real * other.Real + (int128_t)other.Imaginary * other.Imaginary);
		int128_t Real_result = (Real > other.Real ?
			(((int128_t)Real << Fixed_shift) / denominator * other.Real) :
			(((int128_t)other.Real << Fixed_shift) / denominator * Real))
			+ (Imaginary > other.Imaginary ?
				(((int128_t)Imaginary << Fixed_shift) / denominator * other.Imaginary) :
				(((int128_t)other.Imaginary << Fixed_shift) / denominator * Imaginary)),
			Imaginary_result = (Imaginary > other.Real ?
				(((int128_t)Imaginary << Fixed_shift) / denominator * other.Real) :
				(((int128_t)other.Real << Fixed_shift) / denominator * Imaginary))
			- (Real > other.Imaginary ?
				(((int128_t)Real << Fixed_shift) / denominator * other.Imaginary) :
				(((int128_t)other.Imaginary << Fixed_shift) / denominator * Real));

		return FixedComplex(static_cast<long long>(Real_result), static_cast<long long>(Imaginary_result));
	}

	static inline long long int AbsoluteValue(const FixedComplex& C) {

		int128_t Result = ((int128_t)C.Real * C.Real + (int128_t)C.Imaginary * C.Imaginary);

		return NewtonSqrt(Result);

	}

	static inline long long int NewtonSqrt(const int128_t input) {

		if (input == 0)
			return 0;

		int128_t Result = input;
		int128_t x = (Result + 1) / 2;

		while (x < Result) {

			Result = x;
			x = (Result + input / Result) / 2;

		}

		return static_cast<long long int>(Result);

	}

};

class Complex {

public:

	double Real = 0.0, Imaginary = 0.0;

	Complex() :Real(0), Imaginary(0) {};
	explicit Complex(double A) :Real(A), Imaginary(0) {};
	explicit Complex(double A, double B) :Real(A), Imaginary(B) {};

	friend Complex operator+(Complex& c, const double d)noexcept {
		return Complex(c.Real + d, c.Imaginary);
	}
	friend Complex operator+(const double d, Complex& c)noexcept {
		return Complex(c.Real + d, c.Imaginary);
	}
	Complex operator+(const Complex& other)const& noexcept {
		return Complex(Real + other.Real, Imaginary + other.Imaginary);
	}


	friend Complex operator-(const Complex& c, const double d)noexcept {
		return Complex(c.Real - d, c.Imaginary);
	}
	friend Complex operator-(const double d, const Complex& c)noexcept {
		return Complex(d - c.Real, c.Imaginary);
	}
	Complex operator-(const Complex& other)const& noexcept {
		return Complex(Real - other.Real, Imaginary - other.Imaginary);
	}


	friend Complex operator*(const Complex& c, const double d)noexcept {
		return Complex(c.Real * d, c.Imaginary * d);
	}
	friend Complex operator*(const double d, const Complex& c)noexcept {
		return Complex(c.Real * d, c.Imaginary * d);
	}
	Complex operator*(const Complex& other)const& noexcept {
		return Complex(Real * other.Real - Imaginary * other.Imaginary, Real * other.Imaginary + Imaginary * other.Real);
	}


	friend Complex operator/(const Complex& c, const double d)noexcept {
		if (std::abs(d) == 0.0)
			return Complex();
		return Complex(c.Real / d, c.Imaginary / d);
	}
	friend Complex operator/(const double& d, const Complex& c)noexcept {
		if (IsZero(c))
			return Complex();
		return Complex((d * c.Real) / (c.Real * c.Real + c.Imaginary * c.Imaginary),
			(-d * c.Imaginary) / (c.Real * c.Real + c.Imaginary * c.Imaginary));
	}
	Complex operator/(const Complex& other)const& noexcept {
		if (IsZero(other))
			return Complex();
		return Complex((Real * other.Real + Imaginary * other.Imaginary) / (other.Real * other.Real + other.Imaginary * other.Imaginary),
			(Imaginary * other.Real - Real * other.Imaginary) / (other.Real * other.Real + other.Imaginary * other.Imaginary));
	}

	inline void InputComplex() {

		std::cin >> Real >> Imaginary;

	}
	static inline void OutputComplex(const Complex& C) {

		std::cout << C.Real << ',' << C.Imaginary;

	}

	static inline bool IsZero(const Complex& C) {

		return ((std::abs(C.Real) < Double_Epsilon) && (std::abs(C.Imaginary) < Double_Epsilon));

	}

	static inline double AbsoluteValue(const Complex& C) {
		return sqrt((C.Real * C.Real) + (C.Imaginary * C.Imaginary));
	}

};

class Matrix {

public:

	std::vector<std::vector<Complex>> data = { {} };

	Matrix() :data({ {} }) {};
	Matrix(const std::vector<double>& A) {

		std::vector<std::vector<Complex>> Result(1, std::vector<Complex>(A.size()));

		for (int y = 0; y<int(A.size()); y++)
			Result[0][y] = Complex(A[y]);

		data = Result;

	};
	Matrix(const std::vector<Complex>& A) :data({ A }) {};
	Matrix(const std::vector<std::vector<double>>& A) {

		if (A.size() == 0)
			data = { {} };

		std::vector<std::vector<Complex>> Result(A.size(), std::vector<Complex>(A[0].size()));

		for (int x = 0; x<int(A.size()); x++)
			for (int y = 0; y<int(A[x].size()); y++)
				Result[x][y] = Complex(A[x][y]);

		data = Result;

	};
	Matrix(const std::vector<std::vector<Complex>>& A) :data(A) {};

	std::vector<Complex>& operator[](const int x) { return data[x]; };
	const std::vector<Complex>& operator[](const int x)const { return data[x]; };

	bool IsMatrix() const {

		if (data.size() == 0)
			return 0;

		int Column = data[0].size();
		for (int x = 0; x<int(data.size()); x++)
			if (int(data[x].size()) != Column)
				return 0;

		return 1;

	}
	bool IsSquareMatrix() const {

		return (IsMatrix() && data.size() == data[0].size());

	}

	Matrix operator+(const Matrix& other)const& noexcept {

		if (!IsMatrix() || !other.IsMatrix())
			return {};

		if ((data.size() != other.data.size()) || (data.size() == 0))
			return {};

		if (data[0].size() != other.data[0].size())
			return {};

		std::vector<std::vector<Complex>> Result(data.size(), std::vector<Complex>(data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = data[x][y] + other.data[x][y];

		return Matrix(Result);

	}
	friend Matrix operator+(const Matrix& A, const std::vector<std::vector<double>>& B) {
		return A + Matrix(B);
	}
	friend Matrix operator+(const std::vector<std::vector<double>>& B, const Matrix& A) {
		return Matrix(B) + A;
	}

	Matrix operator-(const Matrix& other)const& noexcept {

		if (!IsMatrix() || !other.IsMatrix())
			return {};

		if ((data.size() != other.data.size()) || (data.size() == 0))
			return {};

		if (data[0].size() != other.data[0].size())
			return {};

		std::vector<std::vector<Complex>> Result(data.size(), std::vector<Complex>(data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = data[x][y] - other.data[x][y];

		return Matrix(Result);

	}
	friend Matrix operator-(const Matrix& A, const std::vector<std::vector<double>>& B) {
		return A - Matrix(B);
	}
	friend Matrix operator-(const std::vector<std::vector<double>>& B, const Matrix& A) {
		return Matrix(B) - A;
	}


	friend Matrix operator*(const Matrix& A, const double Scalar)noexcept {

		if (!A.IsMatrix())
			return {};

		std::vector<std::vector<Complex>> Result(A.data.size(), std::vector<Complex>(A.data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = A.data[x][y] * Scalar;

		return Matrix(Result);

	}
	friend Matrix operator*(const double Scalar, const Matrix& A)noexcept {

		return A * Scalar;

	}
	friend Matrix operator*(const Matrix& A, const Complex& Scalar)noexcept {

		if (!A.IsMatrix())
			return {};

		std::vector<std::vector<Complex>> Result(A.data.size(), std::vector<Complex>(A.data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = A.data[x][y] * Scalar;

		return Matrix(Result);

	}
	friend Matrix operator*(const Complex& Scalar, const Matrix& A)noexcept {

		return A * Scalar;

	}
	Matrix operator*(const Matrix& other)const& noexcept {

		if (!IsMatrix() || !other.IsMatrix())
			return {};

		if (data[0].size() != other.data.size())
			return {};

		std::vector<std::vector<Complex>> Result(data.size(), std::vector<Complex>(other.data[0].size()));
		for (int Leftx = 0; Leftx<int(data.size()); Leftx++)
			for (int Righty = 0; Righty<int(other.data[0].size()); Righty++) {

				Complex sum = Complex(0);
				for (int count = 0; count<int(data[0].size()); count++)
					sum = sum + (data[Leftx][count] * other[count][Righty]);

				Result[Leftx][Righty] = sum;

			}

		return Matrix(Result);

	}
	friend Matrix operator*(const Matrix& A, const std::vector<std::vector<double>> B)noexcept {

		return A * Matrix(B);

	}
	friend Matrix operator*(const std::vector<std::vector<double>> B, const Matrix& A)noexcept {

		return Matrix(B) * A;

	}

	friend Matrix operator/(const Matrix& A, const double Scalar)noexcept {

		if (!A.IsMatrix())
			return {};

		std::vector<std::vector<Complex>> Result(A.data.size(), std::vector<Complex>(A.data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = A.data[x][y] / Scalar;

		return Matrix(Result);

	}
	friend Matrix operator/(const Matrix& A, const Complex& Scalar)noexcept {

		if (!A.IsMatrix())
			return {};

		std::vector<std::vector<Complex>> Result(A.data.size(), std::vector<Complex>(A.data[0].size()));
		for (int x = 0; x<int(Result.size()); x++)
			for (int y = 0; y<int(Result[x].size()); y++)
				Result[x][y] = A.data[x][y] / Scalar;

		return Matrix(Result);

	}

	inline void InputMatrix() {

		int Row, Column;
		std::cout << "The number of the Row and the Column" << std::endl;
		std::cin >> Row >> Column;

		data = std::vector<std::vector<Complex>>(Row, std::vector<Complex>(Column, Complex()));
		for (int x = 0; x < Row; x++)
			for (int y = 0; y < Column; y++)
				data[x][y].InputComplex();

	}
	static inline void OutputMatrix(const Matrix& A) {

		if (!A.IsMatrix())
			return;

		for (int x = 0; x<int(A.data.size()); x++) {
			for (int y = 0; y<int(A[x].size()); y++)
				Complex::OutputComplex(A[x][y]), std::cout << ' ' << ' ';
			std::cout << std::endl;
		}

		std::cout << std::endl;

	}

};

class QMath {

public:

	static inline bool IsZero(const Complex& C) {

		return ((AbsoluteValue(C.Real) < Double_Epsilon) && (AbsoluteValue(C.Imaginary) < Double_Epsilon));

	}

	static inline double AbsoluteValue(const double& D) {

		return (D > 0 ? D : -D);

	}

	static inline long long int AbsoluteValue(const FixedComplex& C) {

		int128_t Result = ((int128_t)C.Real * C.Real + (int128_t)C.Imaginary * C.Imaginary);

		return NewtonSqrt(Result);

	}

	static inline double AbsoluteValue(const Complex& C) {
		return sqrt((C.Real * C.Real) + (C.Imaginary * C.Imaginary));
	}

	static inline long long int NewtonSqrt(const int128_t input) {

		if (input == 0)
			return 0;

		int128_t Result = input;
		int128_t x = (Result + 1) / 2;

		while (x < Result) {

			Result = x;
			x = (Result + input / Result) / 2;

		}

		return static_cast<long long int>(Result);

	}

	static inline Complex Conjugate(const Complex& C) {

		return Complex(C.Real, -C.Imaginary);

	}

	static inline bool IsSameMatrix(const Matrix& Matrix1, const Matrix& Matrix2) {


		if (Matrix1.data.size() != Matrix2.data.size() || !Matrix1.IsMatrix() || !Matrix2.IsMatrix())
			return 0;

		if (Matrix1[0].size() != Matrix2[0].size())
			return 0;

		for (int x = 0; x<int(Matrix1.data.size()); x++)
			for (int y = 0; y<int(Matrix1[x].size()); y++)
				if (!IsZero(Matrix1[x][y] - Matrix2[x][y]))
					return 0;

		return 1;

	}

	static inline bool IsDiagonalMatrix(const Matrix& A) {

		if (!A.IsSquareMatrix())
			return 0;

		for (int x = 0; x<int(A.data.size()); x++)
			for (int y = 0; y<int(A[0].size()); y++)
				if (x != y && !IsZero(A[x][y]))
					return 0;

		return 1;

	}

	static inline double NormOfVector(const std::vector<Complex>& Vector) {

		double Sum = 0;
		for (int count = 0; count<int(Vector.size()); count++)
			Sum += (Vector[count].Real * Vector[count].Real) + (Vector[count].Imaginary * Vector[count].Imaginary);

		return sqrt(Sum);

	}

	static inline std::vector<Complex> ScalarMultiplicationVector(const Complex& Scalar, const std::vector<Complex>& Vector) {

		std::vector<Complex> Result(Vector.size());
		for (int count = 0; count<int(Vector.size()); count++)
			Result[count] = Vector[count] * Scalar;

		return Result;

	}

	static inline Matrix GenerateReflectionMatrix(const std::vector<Complex>& NormalVector) {

		return GenerateUnitMatrix(NormalVector.size()) - ((2 / DotProduct(NormalVector, NormalVector)) * OuterProduct(NormalVector, NormalVector));

	}

	static inline Complex DotProduct(const std::vector<Complex>& Vector1, const std::vector<Complex>& Vector2) {

		if (Vector1.size() != Vector2.size())
			return Complex();

		Complex Sum = Complex();
		for (int count = 0; count<int(Vector1.size()); count++)
			Sum = Sum + Conjugate(Vector1[count]) * Vector2[count];

		return Sum;

	}

	static inline Matrix OuterProduct(const std::vector<Complex>& Vector1, const std::vector<Complex>& Vector2) {

		if (Vector1.size() != Vector2.size())
			return {};

		return Matrix(Vector1) * TransposeMatrix(Matrix(Vector2));

	}

	static Matrix GenerateUnitMatrix(const int& order) {

		Matrix UnitMatrix = Matrix(std::vector<std::vector<Complex>>(order, std::vector<Complex>(order, Complex())));
		for (int count = 0; count < order; count++)
			UnitMatrix[count][count] = Complex(1);

		return UnitMatrix;

	}

	static inline void AdditionOneRowToAnotherRow(const Complex& Scalar, const int& AddedRow, const int& Row, Matrix& A) {

		if (!A.IsMatrix())
			return;

		if (AddedRow >= int(A.data.size()) || AddedRow < 0 || Row >= int(A.data.size()) || Row < 0 || Row == AddedRow)
			return;

		for (int y = 0; y<int(A[0].size()); y++)
			A[Row][y] = A[Row][y] + A[AddedRow][y] * Scalar;

	}

	static inline void MultipleOfARow(const Complex& Scalar, const int& Row, Matrix& A) {

		if (!A.IsMatrix())
			return;

		if (Row >= int(A.data.size()) || Row < 0 || (Scalar.Real == 0 && Scalar.Imaginary == 0))
			return;

		for (int y = 0; y<int(A[Row].size()); y++)
			A[Row][y] = A[Row][y] * Scalar;

	}

	static inline void InterchangeOfTwoRow(const std::pair<int, int>& InterchangeRow, Matrix& A) {

		if (!A.IsMatrix())
			return;

		if (InterchangeRow.first >= int(A.data.size()) || InterchangeRow.first < 0 || InterchangeRow.second >= int(A.data.size()) || InterchangeRow.second < 0)
			return;

		std::swap(A[InterchangeRow.first], A[InterchangeRow.second]);

	}

	static inline Matrix HermitianTransposeMatrix(Matrix A) {

		if (!A.IsSquareMatrix())
			return Matrix();

		for (int x = 0; x<int(A.data.size()); x++) {

			A[x][x] = Conjugate(A[x][x]);

			for (int y = 0; y < x; y++)
				std::swap(A[x][y], A[y][x]),
				A[x][y] = Conjugate(A[x][y]),
				A[y][x] = Conjugate(A[y][x]);

		}

		return A;

	}

	static inline Matrix TransposeMatrix(Matrix A) {

		if (!A.IsSquareMatrix())
			return Matrix();

		for (int x = 0; x<int(A.data.size()); x++)
			for (int y = 0; y < x; y++)
				std::swap(A[x][y], A[y][x]);

		return A;

	}

	static Complex DeterminantsOfMatrix(const Matrix& A) {

		if (!A.IsSquareMatrix())
			return Complex();

		if (A.data.size() == 1)
			return A[0][0];

		if (A.data.size() == 2)
			return A[0][0] * A[1][1] - A[0][1] * A[1][0];

		Complex Sum = Complex();

		for (int y = 0; y<int(A[0].size()); y++) {

			Matrix ReducedMatrix = Matrix(std::vector<std::vector<Complex>>(A.data.size() - 1));
			Complex now = A[0][y];

			for (int x = 1; x<int(A.data.size()); x++)
				for (int Column = 0; Column<int(A[0].size()); Column++)
					if (y != Column)
						ReducedMatrix[x - 1].push_back(A[x][Column]);

			if (y % 2 == 1)
				now = now - 1;

			Sum = Sum + now * DeterminantsOfMatrix(ReducedMatrix);

		}

		return Sum;

	}

	static Matrix GenerateInverseMatrix(const Matrix& A) {

		if (!A.IsSquareMatrix())
			return {};

		if (IsZero(DeterminantsOfMatrix(A)) || A[0].size() == 0)
			return {};

		Matrix Now = A, Inverse = GenerateUnitMatrix(A.data.size());
		for (int y = 0; y<int(Now.data.size()); y++) {

			if (IsZero(Now[y][y])) {

				for (int x = 0; x<int(Now.data.size()); x++)
					if (!IsZero(Now[x][y])) {
						InterchangeOfTwoRow({ x,y }, Now);
						InterchangeOfTwoRow({ x,y }, Inverse);
						break;
					}

			}

			MultipleOfARow(1 / Now[y][y], y, Inverse);
			MultipleOfARow(1 / Now[y][y], y, Now);

			for (int x = 0; x<int(Now.data.size()); x++)
				if (x != y && !IsZero(Now[x][y])) {
					AdditionOneRowToAnotherRow(Complex() - Now[x][y], y, x, Inverse);
					AdditionOneRowToAnotherRow(Complex() - Now[x][y], y, x, Now);
				}

		}

		return Inverse;

	}

	/*static Matrix GenerateQuantumLogicGate(const Matrix& A) {

	}*/

	static Matrix GenerateDiagonalMatrix(const Matrix& A) {

		if (!A.IsSquareMatrix() || A.data.size() == 0)
			return {};

		if (A.data.size() == 1)
			return A;

		Matrix Now = A, Front = Matrix();
		std::pair<Matrix, Matrix> QR;
		while (!IsDiagonalMatrix(Now) && !IsSameMatrix(Front, Now)) {

			Front = Now;
			QR = GenerateQRFactorization(Now);
			Now = QR.second * QR.first;

		}

		for (int x = 0; x<int(Now.data.size()); x++)
			for (int y = 0; y<int(Now[x].size()); y++)
				if (x != y)
					Now[x][y] = Complex();

		return Now;

	}

	static std::pair<Matrix, Matrix> GenerateQRFactorization(const Matrix& A) {

		if (!A.IsSquareMatrix() || A.data.size() == 0 || A.data.size() == 1)
			return { {},{} };


		Matrix Q = GenerateUnitMatrix(A.data.size()), R = A;
		for (int y = 0; y<int(A.data.size() - 1); y++) {

			std::vector<Complex> Set(A.data.size() - y);
			for (int x = y; x<int(A.data.size()); x++)
				Set[x - y] = R[x][y];


			Set[0] = Set[0] - NormOfVector(Set);
			Set = ScalarMultiplicationVector(Complex(1 / NormOfVector(Set)), Set);
			Matrix Caculate = GenerateReflectionMatrix(Set);
			Caculate = GenerateBlockEmbeddingForHouseholder(A.data.size(), Caculate);

			Q = Q * Caculate;
			R = Caculate * R;

		}

		return { Q,R };

	}

	static Matrix GenerateBlockEmbeddingForHouseholder(const int& order, const Matrix& A) {

		if (!A.IsSquareMatrix() || int(A.data.size()) > order)
			return {};

		if (order == int(A.data.size()))
			return A;

		Matrix Result = Matrix(std::vector<std::vector<Complex>>(order, std::vector<Complex>(order, Complex())));
		for (int count = 0; count < order; count++)
			Result[count][count] = Complex(1);

		for (int x = order - A.data.size(); x < order; x++)
			for (int y = order - A.data.size(); y < order; y++)
				Result[x][y] = A[x - order + A.data.size()][y - order + A.data.size()];

		return Result;

	}





};

class Qubit_Simulation {

public:

	static long long int FixedPoint(const double input) {

		return static_cast<long long>(input * Fixed_Point);

	}
	static double FixedPointToDouble(const long long int input) {

		double result = 1.0 * input / Fixed_Point;

		return result;

	}

private:

	int Qubit_Amount;
	std::vector<std::pair<bool, int>> Qubit_Set_Observation;
	std::unordered_map<int, std::pair<std::vector<int>, std::vector<FixedComplex>>> Entangled_Qubit_Set_Pointer;
	std::unordered_map<int, std::pair<std::vector<int>, std::vector<FixedComplex>>*> Entangled_Qubit_Set;

	void CombineEntangledQubitSet(int situation1, int situation2) {

		if (Entangled_Qubit_Set.find(situation1) == Entangled_Qubit_Set.end() || Entangled_Qubit_Set.find(situation2) == Entangled_Qubit_Set.end())
			return;

		std::vector<int> First_Qubit_Name_Set = Entangled_Qubit_Set[situation1]->first,
			Second_Qubit_Name_Set = Entangled_Qubit_Set[situation2]->first;

		for (unsigned long long count = 0; count < Second_Qubit_Name_Set.size(); count++)
			First_Qubit_Name_Set.push_back(Second_Qubit_Name_Set[count]);

		Entangled_Qubit_Set[situation1]->first = First_Qubit_Name_Set;

		std::vector<FixedComplex> First_Qubit_Set = Entangled_Qubit_Set[situation1]->second,
			Second_Qubit_Set = Entangled_Qubit_Set[situation2]->second,
			Result(First_Qubit_Set.size() * Second_Qubit_Set.size(), FixedComplex());

		for (unsigned long long x = 0; x < First_Qubit_Set.size(); x++)
			for (unsigned long long y = 0; y < Second_Qubit_Set.size(); y++)
				if ((First_Qubit_Set[x].Real != 0 || First_Qubit_Set[x].Imaginary != 0) && (Second_Qubit_Set[y].Real != 0 || Second_Qubit_Set[y].Imaginary != 0))
					Result[x * Second_Qubit_Set.size() + y] = First_Qubit_Set[x] * Second_Qubit_Set[y];

		Entangled_Qubit_Set[situation1]->second = Result;

		for (unsigned long long count = 0; count < Second_Qubit_Name_Set.size(); count++)
			if (Second_Qubit_Name_Set[count] != situation1)
				Entangled_Qubit_Set_Pointer.erase(Second_Qubit_Name_Set[count]);

		for (unsigned long long count = 0; count < Second_Qubit_Name_Set.size(); count++)
			Entangled_Qubit_Set[Second_Qubit_Name_Set[count]] = Entangled_Qubit_Set[situation1];

	}

	constexpr inline unsigned long long BitRemova(unsigned long long num, const unsigned long long situation) {

		return (num & ((1ULL << situation) - 1ULL)) | ((num >> (situation + 1ULL)) << situation);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num, unsigned long long situation1, unsigned long long situation2) {

		if (situation1 < situation2) {
			situation1 = situation1 ^ situation2;
			situation2 = situation1 ^ situation2;
			situation1 = situation1 ^ situation2;
		}

		num = (num & ((1ULL << situation1) - 1ULL)) | ((num >> (situation1 + 1ULL)) << situation1);

		return (num & ((1ULL << situation2) - 1ULL)) | ((num >> (situation2 + 1ULL)) << situation2);

	}
	constexpr inline unsigned long long BitRemova(unsigned long long num, unsigned long long situation1, unsigned long long situation2, unsigned long long situation3) {

		if (situation1 < situation2) {
			situation1 = situation1 ^ situation2;
			situation2 = situation1 ^ situation2;
			situation1 = situation1 ^ situation2;
		}
		if (situation1 < situation3) {
			situation1 = situation1 ^ situation3;
			situation3 = situation1 ^ situation3;
			situation1 = situation1 ^ situation3;
		}

		num = (num & ((1ULL << situation1) - 1ULL)) | ((num >> (situation1 + 1ULL)) << situation1);

		if (situation2 < situation3) {
			situation2 = situation2 ^ situation3;
			situation3 = situation2 ^ situation3;
			situation2 = situation2 ^ situation3;
		}
		num = (num & ((1ULL << situation2) - 1ULL)) | ((num >> (situation2 + 1ULL)) << situation2);

		return (num & ((1ULL << situation3) - 1ULL)) | ((num >> (situation3 + 1ULL)) << situation3);

	}

	constexpr inline unsigned long long BitAdd(unsigned long long num, const int situation, bool add) {

		return ((num >> situation) << (situation + 1ULL)) | (add << situation) | (num & ((1ULL << situation) - 1ULL));

	}
	inline unsigned long long BitAdd(unsigned long long num, std::vector<std::pair<int, bool>> Control) {

		std::sort(Control.begin(), Control.end());

		for (int count = 0; count<int(Control.size()); count++)
			num = BitAdd(num, Control[count].first, Control[count].second);

		return num;

	}

	int GetSituation(const int situation) {

		std::vector<int> Qubit_Set = Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation = -1;
		for (unsigned long long count = 0; count < Qubit_Set.size(); count++)
			if (Qubit_Set[count] == situation)
				Qubit_Situation = count;

		if (Qubit_Situation == -1)
			return -1;

		return Qubit_Set.size() - Qubit_Situation - 1;

	}

	bool IsQubitUnoperable(const int situation) {

		return (Entangled_Qubit_Set.find(situation) == Entangled_Qubit_Set.end() || IsObservered(situation));

	}

	bool IsSeparable(const int qubit) {

		if (IsQubitUnoperable(qubit)) return false;

		// 只有一個 qubit 的 set 本來就是獨立的
		if (Entangled_Qubit_Set[qubit]->first.size() == 1) return true;

		int q_pos = GetSituation(qubit);
		std::vector<FixedComplex>& state = Entangled_Qubit_Set[qubit]->second;
		size_t half = state.size() >> 1;

		// 找第一個非零的參考振幅對
		long long ref0_R = 0, ref0_I = 0;
		long long ref1_R = 0, ref1_I = 0;
		bool found_ref = false;

		for (size_t count = 0; count < half; count++) {

			unsigned long long idx0 = BitAdd(count, q_pos, 0ULL);
			unsigned long long idx1 = BitAdd(count, q_pos, 1ULL);

			FixedComplex amp0 = state[idx0];
			FixedComplex amp1 = state[idx1];

			bool amp0_zero = (amp0.Real == 0 && amp0.Imaginary == 0);
			bool amp1_zero = (amp1.Real == 0 && amp1.Imaginary == 0);

			if (!found_ref) {
				if (!amp0_zero || !amp1_zero) {
					ref0_R = amp0.Real; ref0_I = amp0.Imaginary;
					ref1_R = amp1.Real; ref1_I = amp1.Imaginary;
					found_ref = true;
				}
				continue;
			}

			// 驗證交叉乘積：amp0 * ref1 == amp1 * ref0
			// (amp0_R + i·amp0_I)(ref1_R - i·ref1_I)
			// == (amp1_R + i·amp1_I)(ref0_R - i·ref0_I)
			int128_t lhs_R = (int128_t)amp0.Real * ref1_R + (int128_t)amp0.Imaginary * ref1_I;
			int128_t lhs_I = (int128_t)amp0.Imaginary * ref1_R - (int128_t)amp0.Real * ref1_I;

			int128_t rhs_R = (int128_t)amp1.Real * ref0_R + (int128_t)amp1.Imaginary * ref0_I;
			int128_t rhs_I = (int128_t)amp1.Imaginary * ref0_R - (int128_t)amp1.Real * ref0_I;

			if (QMath::NewtonSqrt((lhs_R - rhs_R) * (lhs_R - rhs_R) + (lhs_I - rhs_I) * (lhs_I - rhs_I)) > Tolerent_Round)
				return false;
		}

		return true;
	}

	void ExtractSeparatedQubit(const int qubit) {

		if (!IsSeparable(qubit)) return;
		if (Entangled_Qubit_Set[qubit]->first.size() == 1) return;

		int q_pos = GetSituation(qubit);
		std::vector<FixedComplex>& state = Entangled_Qubit_Set[qubit]->second;
		size_t half = state.size() >> 1;


		int128_t P0 = 0, P1 = 0;
		for (size_t count = 0; count < half; count++) {
			FixedComplex amp0 = state[BitAdd(count, q_pos, 0ULL)];
			FixedComplex amp1 = state[BitAdd(count, q_pos, 1ULL)];

			int128_t abs0 = QMath::AbsoluteValue(amp0);
			int128_t abs1 = QMath::AbsoluteValue(amp1);
			P0 += (abs0 * abs0) >> Fixed_shift;
			P1 += (abs1 * abs1) >> Fixed_shift;
		}

		long long abs_a = QMath::NewtonSqrt(P0 << Fixed_shift);
		long long abs_b = QMath::NewtonSqrt(P1 << Fixed_shift);

		FixedComplex qubit_amp0 = FixedComplex(abs_a);
		FixedComplex qubit_amp1 = FixedComplex();

		for (size_t count = 0; count < half; count++) {
			FixedComplex amp0 = state[BitAdd(count, q_pos, 0ULL)];
			FixedComplex amp1 = state[BitAdd(count, q_pos, 1ULL)];

			bool amp0_zero = (amp0.Real == 0 && amp0.Imaginary == 0);

			if (!amp0_zero) {

				qubit_amp1 = (amp1 / amp0) * abs_a;
				break;
			}

			if (amp0_zero && (amp1.Real != 0 || amp1.Imaginary != 0)) {
				qubit_amp0 = FixedComplex(0LL);
				qubit_amp1 = FixedComplex(abs_b);
				break;
			}
		}

		std::vector<FixedComplex> new_state(half, FixedComplex());
		bool use_a = (qubit_amp0.Real != 0 || qubit_amp0.Imaginary != 0);

		for (size_t count = 0; count < half; count++) {
			FixedComplex amp = use_a
				? state[BitAdd(count, q_pos, 0ULL)]
				: state[BitAdd(count, q_pos, 1ULL)];

			FixedComplex divisor = use_a ? qubit_amp0 : qubit_amp1;
			new_state[count] = amp / divisor;
		}

		std::vector<int>& qubit_names = Entangled_Qubit_Set[qubit]->first;
		qubit_names.erase(
			std::remove(qubit_names.begin(), qubit_names.end(), qubit),
			qubit_names.end()
		);

		Entangled_Qubit_Set[qubit]->second = new_state;

		Entangled_Qubit_Set_Pointer[qubit] = {
			{ qubit },
			{ qubit_amp0, qubit_amp1 }
		};
		Entangled_Qubit_Set[qubit] = &Entangled_Qubit_Set_Pointer[qubit];
	}

public:

	Qubit_Simulation() :Qubit_Amount(0) { Entangled_Qubit_Set.clear(); Qubit_Set_Observation.clear(); Entangled_Qubit_Set_Pointer.clear(); };
	Qubit_Simulation(int num) {

		Qubit_Amount = num;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

		for (int count = 0; count < num; count++)
			Entangled_Qubit_Set_Pointer[count] = { {count},{FixedComplex(Fixed_Point),FixedComplex()} },
			Qubit_Set_Observation.push_back({ 0,-1 }),
			Entangled_Qubit_Set[count] = &Entangled_Qubit_Set_Pointer[count];

	};

	void ResetQubitSet() {

		Qubit_Amount = 0;
		Entangled_Qubit_Set.clear();
		Qubit_Set_Observation.clear();
		Entangled_Qubit_Set_Pointer.clear();

	}

	void OuputEntangledQubitSet(const int situation) {

		if (Entangled_Qubit_Set.find(situation) == Entangled_Qubit_Set.end())
			return;

		if (Entangled_Qubit_Set[situation] == nullptr)
			return;

		std::vector<FixedComplex> Output_Qubit = Entangled_Qubit_Set[situation]->second;
		for (unsigned long long count = 0; count < Output_Qubit.size(); count++)
			std::cout << FixedPointToDouble(Output_Qubit[count].Real) << ' ' << FixedPointToDouble(Output_Qubit[count].Imaginary) << std::endl;
		//std::cout<<Output_Qubit[count].Real<<' '<<Output_Qubit[count].Imaginary<<std::endl;

	}

	bool IsObservered(const int situation) {

		return Qubit_Set_Observation[situation].first;

	}

	void CutQubitSet(const int situation, long long int Probability0, long long int Probability1) {

		if (Probability0 == 0)
			Qubit_Set_Observation[situation] = { 1,1 }, Probability1 = Fixed_Point;
		else if (Probability0 == 1)
			Qubit_Set_Observation[situation] = { 1,0 }, Probability0 = Fixed_Point;
		else {

			std::random_device rd;
			std::mt19937_64 gen(rd());
			std::uniform_int_distribution<long long int> dist(0, Fixed_Point);
			long long int RandNum = dist(gen);

			if (RandNum <= Probability0)
				Qubit_Set_Observation[situation] = { 1,0 };
			else
				Qubit_Set_Observation[situation] = { 1,1 };

		}

		long long int State = Qubit_Set_Observation[situation].second;
		long long int Probability = (State == 0) ? Probability0 : Probability1;

		std::vector<FixedComplex> Qubit_Set_Value = Entangled_Qubit_Set[situation]->second,
			New_Qubit_Set_Value((Qubit_Set_Value.size() >> 1) - (Qubit_Set_Value.size() == 2), FixedComplex());

		std::vector<int> Qubit_Set = Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation = GetSituation(situation);

		for (unsigned long long int count = 0; count < (Qubit_Set_Value.size() >> 1) - (Qubit_Set_Value.size() == 2); count++)
			if (Qubit_Set_Value[BitAdd(count, Qubit_Situation, State)].Real != 0 || Qubit_Set_Value[BitAdd(count, Qubit_Situation, State)].Imaginary != 0)
				New_Qubit_Set_Value[count] = Qubit_Set_Value[BitAdd(count, Qubit_Situation, State)] / QMath::NewtonSqrt((int128_t)Probability << Fixed_shift);

		Entangled_Qubit_Set[situation]->second = New_Qubit_Set_Value;

		Qubit_Set.erase(std::remove(Qubit_Set.begin(), Qubit_Set.end(), situation), Qubit_Set.end());
		Entangled_Qubit_Set[situation]->first = Qubit_Set;

		Entangled_Qubit_Set[situation] = nullptr;

	}

	int ObserverQubit(const int situation) {

		if (IsObservered(situation))
			return Qubit_Set_Observation[situation].second;

		if (Entangled_Qubit_Set.find(situation) == Entangled_Qubit_Set.end())
			return -1;

		std::pair<std::vector<int>, std::vector<FixedComplex>>* Qubit_Set_Pointer = Entangled_Qubit_Set[situation];
		std::vector<int> Qubit_Set = Entangled_Qubit_Set[situation]->first;

		int Qubit_Situation = GetSituation(situation);

		std::vector<FixedComplex> Qubit_Set_Value = Entangled_Qubit_Set[situation]->second;
		long long int Probability0 = 0, Probability1 = 0;

		for (unsigned long long int count = 0; count < Qubit_Set_Value.size() >> 1; count++) {

			if (Qubit_Set_Value[BitAdd(count, Qubit_Situation, 0ULL)].Real != 0 || Qubit_Set_Value[BitAdd(count, Qubit_Situation, 0ULL)].Imaginary != 0) {
				int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Situation, 0ULL)]));
				Caculate *= Caculate;
				Caculate >>= Fixed_shift;
				Probability0 += static_cast<long long>(Caculate);
			}

			if (Qubit_Set_Value[BitAdd(count, Qubit_Situation, 1ULL)].Real != 0 || Qubit_Set_Value[BitAdd(count, Qubit_Situation, 1ULL)].Imaginary != 0) {
				int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Situation, 1ULL)]));
				Caculate *= Caculate;
				Caculate >>= Fixed_shift;
				Probability1 += static_cast<long long>(Caculate);
			}

		}

		CutQubitSet(situation, Probability0, Probability1);
		Qubit_Set = Qubit_Set_Pointer->first;
		Qubit_Set_Value = Qubit_Set_Pointer->second;

		for (long long int setting = Qubit_Set.size() - 1; setting >= 0; setting--) {

			Qubit_Set = Qubit_Set_Pointer->first;
			Qubit_Set_Value = Qubit_Set_Pointer->second;

			Probability0 = 0, Probability1 = 0;

			for (unsigned long long int count = 0; count < Qubit_Set_Value.size() >> 1; count++) {

				if (Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 0ULL)].Real != 0 || Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 0ULL)].Imaginary != 0) {
					int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 0ULL)]));
					Caculate *= Caculate;
					Caculate >>= Fixed_shift;
					Probability0 += static_cast<long long>(Caculate);
				}

				if (Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 1ULL)].Real != 0 || Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 1ULL)].Imaginary != 0) {
					int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 1ULL)]));
					Caculate *= Caculate;
					Caculate >>= Fixed_shift;
					Probability1 += static_cast<long long>(Caculate);
				}

			}

			if (Probability0 == 0 || Probability1 == 0)
				CutQubitSet(Qubit_Set[setting], Probability0, Probability1);

		}

		return Qubit_Set_Observation[situation].second;

	}

	int MeasureParity(std::vector<int> qubits) {

		if (qubits.empty()) return -1;
		for (int q : qubits)
			if (IsQubitUnoperable(q)) return -1;

		// 合併所有 qubit 到同一個 set
		for (size_t i = 1; i < qubits.size(); i++)
			if (Entangled_Qubit_Set[qubits[0]]->first != Entangled_Qubit_Set[qubits[i]]->first)
				CombineEntangledQubitSet(qubits[0], qubits[i]);

		auto& state = Entangled_Qubit_Set[qubits[0]]->second;

		std::vector<int> positions;
		for (int q : qubits)
			positions.push_back(GetSituation(q));

		// 計算偶/奇宇稱的機率
		long long P_even = 0, P_odd = 0;
		for (size_t idx = 0; idx < state.size(); idx++) {

			int parity = 0;
			for (int pos : positions)
				parity ^= ((idx >> pos) & 1);

			int128_t abs_val = QMath::AbsoluteValue(state[idx]);
			long long prob = static_cast<long long>((abs_val * abs_val) >> Fixed_shift);
			(parity == 0 ? P_even : P_odd) += prob;

		}

		// 隨機決定宇稱
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_int_distribution<long long> dist(0, Fixed_Point);
		int result = (dist(gen) <= P_even) ? 0 : 1;
		long long norm = (result == 0) ? P_even : P_odd;
		long long norm_factor = QMath::NewtonSqrt((int128_t)norm << Fixed_shift);

		// 坍縮並歸一化（不坍縮個別 qubit，只投影到宇稱子空間）
		for (size_t idx = 0; idx < state.size(); idx++) {

			int parity = 0;
			for (int pos : positions)
				parity ^= ((idx >> pos) & 1);

			if (parity != result) state[idx] = FixedComplex();
			else                  state[idx] = state[idx] / norm_factor;

		}

		return result;
	}

	int GenerateQubit() {

		Qubit_Amount++;

		Entangled_Qubit_Set_Pointer[Qubit_Amount - 1] = { {Qubit_Amount - 1},{FixedComplex(Fixed_Point),FixedComplex()} };
		Entangled_Qubit_Set[Qubit_Amount - 1] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 1];
		Qubit_Set_Observation.push_back({ 0,-1 });

		return Qubit_Amount - 1;

	}

	int GeneratePositiveBellState() {

		Qubit_Amount += 2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount - 2] = { {Qubit_Amount - 2,Qubit_Amount - 1},{FixedComplex(Root_half),FixedComplex(),FixedComplex(),FixedComplex(Root_half)} };
		Entangled_Qubit_Set[Qubit_Amount - 2] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];
		Entangled_Qubit_Set[Qubit_Amount - 1] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];

		Qubit_Set_Observation.push_back({ 0,-1 });
		Qubit_Set_Observation.push_back({ 0,-1 });
		return Qubit_Amount - 2;

	}

	int GenerateNegativeBellState() {

		Qubit_Amount += 2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount - 2] = { {Qubit_Amount - 2,Qubit_Amount - 1},{FixedComplex(Root_half),FixedComplex(),FixedComplex(),FixedComplex(-Root_half)} };
		Entangled_Qubit_Set[Qubit_Amount - 2] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];
		Entangled_Qubit_Set[Qubit_Amount - 1] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];

		Qubit_Set_Observation.push_back({ 0,-1 });
		Qubit_Set_Observation.push_back({ 0,-1 });
		return Qubit_Amount - 2;

	}

	int GeneratePositiveOrthogonalBellState() {

		Qubit_Amount += 2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount - 2] = { {Qubit_Amount - 2,Qubit_Amount - 1},{FixedComplex(),FixedComplex(Root_half),FixedComplex(Root_half),FixedComplex()} };
		Entangled_Qubit_Set[Qubit_Amount - 2] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];
		Entangled_Qubit_Set[Qubit_Amount - 1] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];

		Qubit_Set_Observation.push_back({ 0,-1 });
		Qubit_Set_Observation.push_back({ 0,-1 });
		return Qubit_Amount - 2;

	}

	int GenerateSingletBellState() {

		Qubit_Amount += 2;
		Entangled_Qubit_Set_Pointer[Qubit_Amount - 2] = { {Qubit_Amount - 2,Qubit_Amount - 1},{FixedComplex(),FixedComplex(Root_half),FixedComplex(-Root_half),FixedComplex()} };
		Entangled_Qubit_Set[Qubit_Amount - 2] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];
		Entangled_Qubit_Set[Qubit_Amount - 1] = &Entangled_Qubit_Set_Pointer[Qubit_Amount - 2];

		Qubit_Set_Observation.push_back({ 0,-1 });
		Qubit_Set_Observation.push_back({ 0,-1 });
		return Qubit_Amount - 2;

	}

	int GetQubitAmount() const { return Qubit_Amount; }

	int  GetMeasurementResult(int qubit) {
		if (!IsObservered(qubit)) return -1;
			return Qubit_Set_Observation[qubit].second;
	}

	std::pair<double, double> GetProbabilities(int qubit) {

		if (IsObservered(qubit)) return {0.0, 0.0};

		if (Entangled_Qubit_Set.find(qubit) == Entangled_Qubit_Set.end()) return {0.0, 0.0};

		auto& state = Entangled_Qubit_Set[qubit]->second;
		int q_pos   = GetSituation(qubit);
		size_t half = state.size() >> 1;
		double P0 = 0.0, P1 = 0.0;

		for (size_t c = 0; c < half; c++) {

			auto a0 = state[BitAdd(c, q_pos, 0ULL)];
			auto a1 = state[BitAdd(c, q_pos, 1ULL)];

			double r0 = FixedPointToDouble(a0.Real), i0 = FixedPointToDouble(a0.Imaginary);
			double r1 = FixedPointToDouble(a1.Real), i1 = FixedPointToDouble(a1.Imaginary);

			P0 += r0*r0 + i0*i0;
			P1 += r1*r1 + i1*i1;

		}

		return {P0, P1};

	}

	std::vector<EntangledGroupSnapshot> GetEntangledGroups() {

		std::vector<EntangledGroupSnapshot> result;

		std::vector<bool> visited(Qubit_Amount, false);

		for (int q = 0; q < Qubit_Amount; q++) {

			if (visited[q] || IsObservered(q)) continue;

			if (Entangled_Qubit_Set.find(q) == Entangled_Qubit_Set.end()) continue;

			auto* ptr = Entangled_Qubit_Set[q];

		    if (!ptr) continue;
		         // 只處理 primary（first qubit id == q）
		    if (ptr->first.empty() || ptr->first[0] != q) continue;

			EntangledGroupSnapshot info;

			info.qubit_ids = ptr->first;

			for (int id : info.qubit_ids) visited[id] = true;
			for (const auto& amp : ptr->second)
				info.amplitudes.push_back({
					FixedPointToDouble(amp.Real),
					FixedPointToDouble(amp.Imaginary)
				});
			result.push_back(info);
		}

		return result;

	}


	void ApplySingleQubitGate(const int qubit, const GateType type, const double param = 0.0) {

		if (IsQubitUnoperable(qubit)) return;
		int q_pos = GetSituation(qubit);
		auto& state = Entangled_Qubit_Set[qubit]->second;
		size_t half = state.size() >> 1;

		double half_rad, c, s;
		if (type == GateType::Rx || type == GateType::Ry || type == GateType::Rz)
			half_rad = param * Pi / 180.0 / 2.0, c = cos(half_rad), s = sin(half_rad);

		for (size_t i = 0; i < half; i++) {

			auto& s0 = state[BitAdd(i, q_pos, 0ULL)];
			auto& s1 = state[BitAdd(i, q_pos, 1ULL)];
			FixedComplex tmp;

			switch (type) {

				// ── Permutation ─────────────────────────────
			case GateType::X:
				std::swap(s0, s1);
				break;

				// ── PhaseOnly（不動 s0）────────────────────
			case GateType::Z:
				s1.Real = -s1.Real; s1.Imaginary = -s1.Imaginary;
				break;
			case GateType::S:
				tmp = s1;
				s1.Real = -tmp.Imaginary; s1.Imaginary = tmp.Real;
				break;
			case GateType::Sdg:
				tmp = s1;
				s1.Real = tmp.Imaginary; s1.Imaginary = -tmp.Real;
				break;
			case GateType::T:
				s1 = s1 * FixedComplex(FixedPoint(Root_half), FixedPoint(Root_half));
				break;
			case GateType::Tdg:
				s1 = s1 * FixedComplex(FixedPoint(Root_half), FixedPoint(-Root_half));
				break;

				// ── PhaseOnly（s0 和 s1 都動）──────────────
			case GateType::Rz:
				s0 = s0 * FixedComplex(FixedPoint(c), FixedPoint(-s));
				s1 = s1 * FixedComplex(FixedPoint(c), FixedPoint(s));
				break;

				// ── FullUnitary off-diagonal（Y）───────────
			case GateType::Y:
				std::swap(s0, s1);
				tmp = s0;
				s0.Real = -tmp.Imaginary; s0.Imaginary = tmp.Real;   // × (-i)
				tmp = s1;
				s1.Real = tmp.Imaginary; s1.Imaginary = -tmp.Real;  // × i  (swap後再乘)
				// 校正：Y = [[0,-i],[i,0]]
				// swap後 s0=原s1, s1=原s0
				// s0 × (-i), s1 × i
				break;

				// ── FullUnitary（H, Rx, Ry）─────────────────
			case GateType::H: {
				FixedComplex a = s0, b = s1;
				s0 = a * FixedPoint(Root_half) + b * FixedPoint(Root_half);
				s1 = a * FixedPoint(Root_half) + b * FixedPoint(-Root_half);
				break;
			}
			case GateType::Rx: {
				FixedComplex a = s0, b = s1;
				s0 = a * FixedPoint(c) + b * FixedComplex(0LL, FixedPoint(-s));
				s1 = a * FixedComplex(0LL, FixedPoint(-s)) + b * FixedPoint(c);
				break;
			}
			case GateType::Ry: {
				FixedComplex a = s0, b = s1;
				s0 = a * FixedPoint(c) + b * FixedPoint(-s);
				s1 = a * FixedPoint(s) + b * FixedPoint(c);
				break;
			}
			default: break;
			}
		}
	}

	void ApplyTwoQubitGate(const int q1, const int q2, const GateType type, const double param = 0.0) {

		if (IsQubitUnoperable(q1) || IsQubitUnoperable(q2) || q1 == q2) return;
		if (Entangled_Qubit_Set[q1]->first != Entangled_Qubit_Set[q2]->first)
			CombineEntangledQubitSet(q1, q2);

		int qs1 = GetSituation(q1), qs2 = GetSituation(q2);
		auto& state = Entangled_Qubit_Set[q1]->second;
		size_t quarter = state.size() >> 2;

		for (size_t i = 0; i < quarter; i++) {

			auto& s00 = state[BitAdd(i, { {qs1,0},{qs2,0} })];
			auto& s01 = state[BitAdd(i, { {qs1,0},{qs2,1} })];
			auto& s10 = state[BitAdd(i, { {qs1,1},{qs2,0} })];
			auto& s11 = state[BitAdd(i, { {qs1,1},{qs2,1} })];

			switch (type) {

			case GateType::CNOT:
				std::swap(s10, s11);
				break;

			case GateType::CZ:
				s11.Real = -s11.Real; s11.Imaginary = -s11.Imaginary;
				break;

			case GateType::SWAP:
				std::swap(s01, s10);
				break;

			case GateType::iSWAP: {
				std::swap(s01, s10);
				FixedComplex tmp = s10;
				s10.Real = -tmp.Imaginary; s10.Imaginary = tmp.Real;  // × i
				tmp = s01;
				s01.Real = -tmp.Imaginary; s01.Imaginary = tmp.Real;  // × i
				break;
			}
			case GateType::SqrtSWAP: {
				FixedComplex a = s01, b = s10;
				s01 = a * FixedComplex(Fixed_Point >> 1, Fixed_Point >> 1)
					+ b * FixedComplex(Fixed_Point >> 1, -(Fixed_Point >> 1));
				s10 = a * FixedComplex(Fixed_Point >> 1, -(Fixed_Point >> 1))
					+ b * FixedComplex(Fixed_Point >> 1, Fixed_Point >> 1);
				break;
			}
			default: break;
			}
		}

		std::vector<int> qubits_in_set = Entangled_Qubit_Set[q1]->first;

		for (int count = 0; count < qubits_in_set.size(); count++)
			ExtractSeparatedQubit(qubits_in_set[count]);

	}

	void ApplyThreeQubitGate(const int q1, const int q2, const int q3,
		const GateType type, const double param = 0.0) {

		if (IsQubitUnoperable(q1) || IsQubitUnoperable(q2) || IsQubitUnoperable(q3) ||
			q1 == q2 || q1 == q3 || q2 == q3) return;
		if (Entangled_Qubit_Set[q1]->first != Entangled_Qubit_Set[q2]->first)
			CombineEntangledQubitSet(q1, q2);
		if (Entangled_Qubit_Set[q1]->first != Entangled_Qubit_Set[q3]->first)
			CombineEntangledQubitSet(q1, q3);

		int qs1 = GetSituation(q1), qs2 = GetSituation(q2), qs3 = GetSituation(q3);
		auto& state = Entangled_Qubit_Set[q1]->second;
		size_t eighth = state.size() >> 3;

		for (size_t i = 0; i < eighth; i++) {

			auto& s110 = state[BitAdd(i, { {qs1,1},{qs2,1},{qs3,0} })];
			auto& s111 = state[BitAdd(i, { {qs1,1},{qs2,1},{qs3,1} })];
			auto& s101 = state[BitAdd(i, { {qs1,1},{qs2,0},{qs3,1} })];

			switch (type) {

			case GateType::CCNOT:
				std::swap(s110, s111);
				break;

			case GateType::CSWAP:
				std::swap(s110, s101);
				break;

			case GateType::Deutsch: {
				double half_rad = param * Pi / 180.0 / 2.0;
				double c = cos(half_rad), s = sin(half_rad);
				FixedComplex a = s110, b = s111;
				s110 = a * FixedComplex(0LL, FixedPoint(c))
					+ b * FixedComplex(FixedPoint(s), 0LL);
				s111 = a * FixedComplex(FixedPoint(s), 0LL)
					+ b * FixedComplex(0LL, FixedPoint(-c));
				break;
			}
			default: break;
			}
		}

		std::vector<int> qubits_in_set = Entangled_Qubit_Set[q1]->first;

		for (int count = 0; count < qubits_in_set.size(); count++)
			ExtractSeparatedQubit(qubits_in_set[count]);

	}

	//單量子位元邏輯閘

	// 單 qubit
	void HadamardGate(int q) { ApplySingleQubitGate(q, GateType::H); }
	void XGate(int q) { ApplySingleQubitGate(q, GateType::X); }
	void YGate(int q) { ApplySingleQubitGate(q, GateType::Y); }
	void ZGate(int q) { ApplySingleQubitGate(q, GateType::Z); }
	void SGate(int q) { ApplySingleQubitGate(q, GateType::S); }
	void SdgGate(int q) { ApplySingleQubitGate(q, GateType::Sdg); }
	void TGate(int q) { ApplySingleQubitGate(q, GateType::T); }
	void TdgGate(int q) { ApplySingleQubitGate(q, GateType::Tdg); }
	void RxGate(int q, double a) { ApplySingleQubitGate(q, GateType::Rx, a); }
	void RyGate(int q, double a) { ApplySingleQubitGate(q, GateType::Ry, a); }
	void RzGate(int q, double a) { ApplySingleQubitGate(q, GateType::Rz, a); }
	// 雙 qubit
	void CNOTGate(int c, int t) { ApplyTwoQubitGate(c, t, GateType::CNOT); }
	void CZGate(int c, int t) { ApplyTwoQubitGate(c, t, GateType::CZ); }
	void SWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::SWAP); }
	void ISWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::iSWAP); }
	void SqrtSWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::SqrtSWAP); }
	// 三 qubit
	void CCNOTGate(int c1, int c2, int t) { ApplyThreeQubitGate(c1, c2, t, GateType::CCNOT); }
	void CSWAPGate(int c, int a, int b) { ApplyThreeQubitGate(c, a, b, GateType::CSWAP); }
	void DeutschGate(int c1, int c2, int t, double angle)
	{
		ApplyThreeQubitGate(c1, c2, t, GateType::Deutsch, angle);
	}




public:







};



// ══════════════════════════════════════════════════════════════
// 快照資料結構
// ══════════════════════════════════════════════════════════════

struct StepSnapshot {
	int                                    step_index = 0;
	std::string                            step_label = "";
	std::vector<EntangledGroupSnapshot>    groups = {};
	std::map<int, std::pair<double, double>>probabilities = {};  // qubit → {P0,P1}
	std::map<int, int>                     measurements = {};  // qubit → 結果（-1=未測）
};

struct CircuitResult {
	int                        num_qubits = 0;
	std::vector<StepSnapshot>  snapshots = {};
};

// ══════════════════════════════════════════════════════════════
// CircuitExpression
//
// 格式：
//   <qubit數>, {op, qubit [, qubit...] [, param]}, ...
//
// 範例：
//   "2, {H,0}, {CNOT,0,1}, {M,0}, {M,1}"
//   "3, {H,0}, {Rx,0,90}, {Deutsch,0,1,2,45}, {Parity,0,1,2}"
//   "3, {Label,start}, {H,0}, {If,0,1,{X,1},{Goto,start}}"
// ══════════════════════════════════════════════════════════════
class CircuitExpression {

public:

	int                             num_qubits = 0;
	std::vector<CircuitInstruction> instructions = {};

	// ──────────────────────────────────────────────────────────
	// Parse：字串 → CircuitExpression
	// ──────────────────────────────────────────────────────────
	static CircuitExpression Parse(const std::string& expr) {

		CircuitExpression result;

		std::string s;
		for (char ch : expr)
			if (!std::isspace(static_cast<unsigned char>(ch)))
				s += ch;

		size_t first_brace = s.find('{');
		if (first_brace == std::string::npos)
			throw std::runtime_error("CircuitExpression: 找不到任何 {}");

		result.num_qubits = std::stoi(s.substr(0, first_brace - 1));

		for (const auto& block : SplitTopLevelBlocks(s.substr(first_brace)))
			result.instructions.push_back(ParseBlock(block));

		return result;
	}

	// ──────────────────────────────────────────────────────────
	// Execute：執行電路（不記錄快照）
	// ──────────────────────────────────────────────────────────
	void Execute(Qubit_Simulation& sim) const {

		sim.ResetQubitSet();

		for (int count = 0; count < num_qubits; count++)
			sim.GenerateQubit();

		std::unordered_map<std::string, int> label_map = BuildLabelMap();

		int ip = 0;
		while (ip < static_cast<int>(instructions.size())) {
			int jump = ExecuteOne(instructions[ip], sim, label_map);
			ip = (jump >= 0) ? jump : ip + 1;
		}
	}

	// ──────────────────────────────────────────────────────────
	// ExecuteWithCapture：執行電路並在每步記錄快照
	// ──────────────────────────────────────────────────────────
	CircuitResult ExecuteWithCapture(Qubit_Simulation& sim) const {

		sim.ResetQubitSet();

		for (int count = 0; count < num_qubits; count++)
			sim.GenerateQubit();

		CircuitResult result;
		result.num_qubits = num_qubits;

		std::unordered_map<std::string, int> label_map = BuildLabelMap();

		int ip = 0;
		while (ip < static_cast<int>(instructions.size())) {

			const auto& inst = instructions[ip];
			int jump = ExecuteOne(inst, sim, label_map);

			result.snapshots.push_back(TakeSnapshot(inst, ip, sim));

			ip = (jump >= 0) ? jump : ip + 1;
		}

		return result;
	}

	// ──────────────────────────────────────────────────────────
	// GateTypeName：enum → 字串
	// ──────────────────────────────────────────────────────────
	static std::string GateTypeName(GateType type) {
		switch (type) {
		case GateType::Identity:    return "Identity";
		case GateType::H:           return "H";
		case GateType::X:           return "X";
		case GateType::Y:           return "Y";
		case GateType::Z:           return "Z";
		case GateType::S:           return "S";
		case GateType::Sdg:         return "Sdg";
		case GateType::T:           return "T";
		case GateType::Tdg:         return "Tdg";
		case GateType::Rx:          return "Rx";
		case GateType::Ry:          return "Ry";
		case GateType::Rz:          return "Rz";
		case GateType::CNOT:        return "CNOT";
		case GateType::CZ:          return "CZ";
		case GateType::SWAP:        return "SWAP";
		case GateType::iSWAP:       return "iSWAP";
		case GateType::SqrtSWAP:    return "SqrtSWAP";
		case GateType::CCNOT:       return "CCNOT";
		case GateType::CSWAP:       return "CSWAP";
		case GateType::Deutsch:     return "Deutsch";
		case GateType::Measure:     return "M";
		case GateType::Parity:      return "Parity";
		case GateType::Conditional: return "If";
		case GateType::Label:       return "Label";
		case GateType::Goto:        return "Goto";
		default:                    return "?";
		}
	}

	// ──────────────────────────────────────────────────────────
	// InstructionToString：指令 → 可讀字串
	// ──────────────────────────────────────────────────────────
	static std::string InstructionToString(const CircuitInstruction& inst) {
		std::ostringstream ss;
		ss << "{" << GateTypeName(inst.type);
		if (!inst.label.empty())
			ss << "," << inst.label;
		for (int q : inst.qubits)
			ss << ",Q" << q;
		if (inst.param != 0.0)
			ss << "," << std::fixed << std::setprecision(1) << inst.param;
		if (inst.type == GateType::Conditional)
			ss << ",exp=" << inst.expected_val;
		for (const auto& sub : inst.sub_insts)
			ss << " " << InstructionToString(sub);
		ss << "}";
		return ss.str();
	}

private:

	// ══════════════════════════════════════════════════════════
	// 解析
	// ══════════════════════════════════════════════════════════

	static GateType ParseGateName(const std::string& name) {
		static const std::unordered_map<std::string, GateType> table = {
			{"H",        GateType::H},        {"X",       GateType::X},
			{"Y",        GateType::Y},        {"Z",       GateType::Z},
			{"S",        GateType::S},        {"Sdg",     GateType::Sdg},
			{"T",        GateType::T},        {"Tdg",     GateType::Tdg},
			{"Rx",       GateType::Rx},       {"Ry",      GateType::Ry},
			{"Rz",       GateType::Rz},
			{"CNOT",     GateType::CNOT},     {"CZ",      GateType::CZ},
			{"SWAP",     GateType::SWAP},     {"iSWAP",   GateType::iSWAP},
			{"SqrtSWAP", GateType::SqrtSWAP},
			{"CCNOT",    GateType::CCNOT},    {"CSWAP",   GateType::CSWAP},
			{"Deutsch",  GateType::Deutsch},
			{"M",        GateType::Measure},  {"Measure", GateType::Measure},
			{"Parity",   GateType::Parity},
			{"If",       GateType::Conditional},
			{"Label",    GateType::Label},
			{"Goto",     GateType::Goto},
		};
		auto it = table.find(name);
		return (it != table.end()) ? it->second : GateType::Identity;
	}

	// {qubit 數, param 數}，-1 = 可變 qubit 數（Parity）
	static std::pair<int, int> GateSignature(GateType type) {
		switch (type) {
		case GateType::H:        case GateType::X:
		case GateType::Y:        case GateType::Z:
		case GateType::S:        case GateType::Sdg:
		case GateType::T:        case GateType::Tdg:
		case GateType::Measure:  return { 1, 0 };
		case GateType::Rx:       case GateType::Ry:
		case GateType::Rz:       return { 1, 1 };
		case GateType::CNOT:     case GateType::CZ:
		case GateType::SWAP:     case GateType::iSWAP:
		case GateType::SqrtSWAP: return { 2, 0 };
		case GateType::CCNOT:    case GateType::CSWAP:
			return { 3, 0 };
		case GateType::Deutsch:  return { 3, 1 };
		case GateType::Parity:   return { -1, 0 };
		default:                 return { 0, 0 };
		}
	}

	// 頂層 {} 切分（處理巢狀）
	static std::vector<std::string> SplitTopLevelBlocks(const std::string& s) {
		std::vector<std::string> result;
		int depth = 0;
		std::string cur;
		for (char ch : s) {
			if (ch == '{') {
				if (depth > 0) cur += ch;
				depth++;
			}
			else if (ch == '}') {
				depth--;
				if (depth > 0) cur += ch;
				else {
					if (!cur.empty()) result.push_back(cur);
					cur.clear();
				}
			}
			else {
				if (depth > 0) cur += ch;
			}
		}
		return result;
	}

	// token 切分（不含 {} 的普通逗號分隔）
	static std::vector<std::string> SplitTokens(const std::string& s) {
		std::vector<std::string> tokens;
		std::string cur;
		for (size_t i = 0; i <= s.size(); i++) {
			if (i == s.size() || s[i] == ',') {
				if (!cur.empty()) tokens.push_back(cur);
				cur.clear();
			}
			else {
				cur += s[i];
			}
		}
		return tokens;
	}

	static CircuitInstruction ParseBlock(const std::string& block) {

		size_t first_comma = block.find(',');
		std::string gate_name = (first_comma == std::string::npos)
			? block : block.substr(0, first_comma);

		CircuitInstruction inst;
		inst.type = ParseGateName(gate_name);

		if (first_comma == std::string::npos)
			return inst;

		std::string rest = block.substr(first_comma + 1);

		// Label / Goto
		if (inst.type == GateType::Label || inst.type == GateType::Goto) {
			inst.label = rest;
			return inst;
		}

		// Conditional：If,<qubit>,<expected>,{sub},...
		if (inst.type == GateType::Conditional) {
			size_t first_brace = rest.find('{');
			std::string header = (first_brace == std::string::npos)
				? rest : rest.substr(0, first_brace - 1);
			std::vector<std::string> hdr = SplitTokens(header);
			if (hdr.size() >= 2) {
				inst.qubits.push_back(std::stoi(hdr[0]));
				inst.expected_val = std::stoi(hdr[1]);
			}
			if (first_brace != std::string::npos)
				for (const auto& sub : SplitTopLevelBlocks(rest.substr(first_brace)))
					inst.sub_insts.push_back(ParseBlock(sub));
			return inst;
		}

		// 普通閘
		std::pair<int, int> sig = GateSignature(inst.type);
		int nq = sig.first;
		int np = sig.second;

		if (nq == -1) {
			// Parity：全部 token 都是 qubit
			for (const auto& t : SplitTokens(rest))
				inst.qubits.push_back(std::stoi(t));
		}
		else {
			std::vector<std::string> tokens = SplitTokens(rest);
			for (int i = 0; i < nq && i < static_cast<int>(tokens.size()); i++)
				inst.qubits.push_back(std::stoi(tokens[i]));
			if (np > 0 && nq < static_cast<int>(tokens.size()))
				inst.param = std::stod(tokens[nq]);
		}

		return inst;
	}

	// ══════════════════════════════════════════════════════════
	// 執行
	// ══════════════════════════════════════════════════════════

	std::unordered_map<std::string, int> BuildLabelMap() const {
		std::unordered_map<std::string, int> m;
		for (int i = 0; i < static_cast<int>(instructions.size()); i++)
			if (instructions[i].type == GateType::Label)
				m[instructions[i].label] = i;
		return m;
	}

	// 回傳 -1 = 正常繼續，>= 0 = 跳躍到該 index
	static int ExecuteOne(
		const CircuitInstruction& inst,
		Qubit_Simulation& sim,
		const std::unordered_map<std::string, int>& label_map) {

		switch (inst.type) {

		case GateType::H:    case GateType::X:   case GateType::Y:
		case GateType::Z:    case GateType::S:   case GateType::Sdg:
		case GateType::T:    case GateType::Tdg:
			sim.ApplySingleQubitGate(inst.qubits[0], inst.type);
			return -1;

		case GateType::Rx:   case GateType::Ry:  case GateType::Rz:
			sim.ApplySingleQubitGate(inst.qubits[0], inst.type, inst.param);
			return -1;

		case GateType::CNOT:     case GateType::CZ:
		case GateType::SWAP:     case GateType::iSWAP:
		case GateType::SqrtSWAP:
			sim.ApplyTwoQubitGate(inst.qubits[0], inst.qubits[1], inst.type);
			return -1;

		case GateType::CCNOT:    case GateType::CSWAP:
			sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1],
				inst.qubits[2], inst.type);
			return -1;

		case GateType::Deutsch:
			sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1],
				inst.qubits[2], inst.type, inst.param);
			return -1;

		case GateType::Measure:
			sim.ObserverQubit(inst.qubits[0]);
			return -1;

		case GateType::Parity:
			sim.MeasureParity(inst.qubits);
			return -1;

		case GateType::Label:
			return -1;

		case GateType::Goto: {
			auto it = label_map.find(inst.label);
			return (it != label_map.end()) ? it->second : -1;
		}

		case GateType::Conditional: {
			int meas = sim.ObserverQubit(inst.qubits[0]);
			if (meas == inst.expected_val)
				for (const auto& sub : inst.sub_insts) {
					int jump = ExecuteOne(sub, sim, label_map);
					if (jump >= 0) return jump;
				}
			return -1;
		}

		default:
			return -1;
		}
	}

	// ══════════════════════════════════════════════════════════
	// 快照
	// ══════════════════════════════════════════════════════════

	static StepSnapshot TakeSnapshot(const CircuitInstruction& inst,
		int                        step,
		Qubit_Simulation& sim) {
		StepSnapshot snap;
		snap.step_index = step;
		snap.step_label = InstructionToString(inst);
		snap.groups = sim.GetEntangledGroups();

		int nq = sim.GetQubitAmount();
		for (int q = 0; q < nq; q++) {
			int meas = sim.GetMeasurementResult(q);
			snap.measurements[q] = meas;
			if (meas == -1)
				snap.probabilities[q] = sim.GetProbabilities(q);
		}

		return snap;
	}

};

class CircuitPrinter {

public:

	// ── 主入口：電路圖 + 每步快照 ──────────────────────────────
	static void Print(const CircuitExpression& circuit,
		const CircuitResult& result,
		std::ostream& out = std::cout) {
		PrintHeader(circuit, out);
		PrintCircuitDiagram(circuit, out);
		PrintDivider('═', 50, out);
		for (const auto& snap : result.snapshots)
			PrintSnapshot(snap, result.num_qubits, out);
		PrintDivider('═', 50, out);
	}

	// ── 只印電路圖 ─────────────────────────────────────────────
	static void PrintDiagram(const CircuitExpression& circuit,
		std::ostream& out = std::cout) {
		PrintHeader(circuit, out);
		PrintCircuitDiagram(circuit, out);
	}

	// ── 只印全部快照 ───────────────────────────────────────────
	static void PrintAllSnapshots(const CircuitResult& result,
		std::ostream& out = std::cout) {
		for (const auto& snap : result.snapshots)
			PrintSnapshot(snap, result.num_qubits, out);
	}

	// ── 印單一快照 ─────────────────────────────────────────────
	static void PrintSnapshot(const StepSnapshot& snap,
		int                 num_qubits,
		std::ostream& out = std::cout) {

		PrintDivider('─', 50, out);
		out << " Step " << std::setw(2) << snap.step_index + 1
			<< " │ " << snap.step_label << "\n";
		PrintDivider('─', 50, out);

		// 測量結果行（已坍縮的 qubit）
		bool any_measured = false;
		for (int q = 0; q < num_qubits; q++) {
			auto it = snap.measurements.find(q);
			if (it != snap.measurements.end() && it->second != -1) {
				if (!any_measured) {
					out << " 測量結果\n";
					any_measured = true;
				}
				out << "   Q[" << q << "] = " << it->second << "\n";
			}
		}

		// 各糾纏組的狀態向量
		if (!snap.groups.empty()) {
			out << " 量子態\n";
			for (const auto& group : snap.groups)
				PrintGroup(group, snap.probabilities, out);
		}
	}

private:

	// ══════════════════════════════════════════════════════════
	// 電路圖
	// ══════════════════════════════════════════════════════════

	static void PrintCircuitDiagram(const CircuitExpression& circuit,
		std::ostream& out) {

		int nq = circuit.num_qubits;

		// ── 把指令分配到時間槽 ─────────────────────────────
		struct TimeSlot {
			std::vector<bool>              occupied;
			std::vector<CircuitInstruction>inst_at;
			TimeSlot(int n)
				: occupied(n, false), inst_at(n) {
			}
		};

		std::vector<TimeSlot> slots;

		for (const auto& inst : circuit.instructions) {

			// 收集這條指令涉及的所有 qubit
			std::vector<int> targets = CollectQubits(inst, nq);

			// Label / Goto / Conditional 獨占一欄
			bool standalone = (inst.type == GateType::Label ||
				inst.type == GateType::Goto ||
				inst.type == GateType::Conditional);

			int slot_idx = -1;
			if (!standalone) {
				for (int s = static_cast<int>(slots.size()) - 1; s >= 0; s--) {
					bool conflict = false;
					for (int q : targets)
						if (slots[s].occupied[q]) { conflict = true; break; }
					if (!conflict) { slot_idx = s; break; }
				}
			}

			if (slot_idx == -1) {
				slots.push_back(TimeSlot(nq));
				slot_idx = static_cast<int>(slots.size()) - 1;
			}

			for (int q : targets) {
				slots[slot_idx].occupied[q] = true;
				slots[slot_idx].inst_at[q] = inst;
			}
		}

		// ── 計算欄寬 ──────────────────────────────────────
		const int COL_W = 8;

		// ── 表頭 ──────────────────────────────────────────
		out << "\n";
		out << "       ";
		for (int s = 0; s < static_cast<int>(slots.size()); s++) {
			std::string hdr = "T" + std::to_string(s);
			out << std::left << std::setw(COL_W) << hdr;
		}
		out << "\n";

		// ── 每個 qubit 一行 ───────────────────────────────
		for (int q = 0; q < nq; q++) {
			std::ostringstream row;
			row << " Q[" << q << "]  ";
			for (const auto& slot : slots) {
				if (slot.occupied[q])
					row << std::left << std::setw(COL_W)
					<< GateSymbol(slot.inst_at[q], q);
				else
					row << std::left << std::setw(COL_W) << "──────";
			}
			out << row.str() << "\n";
		}
		out << "\n";
	}

	// 收集一條指令（含 sub_insts）涉及的所有合法 qubit
	static std::vector<int> CollectQubits(const CircuitInstruction& inst,
		int nq) {
		std::vector<int> result;
		for (int q : inst.qubits)
			if (q >= 0 && q < nq) result.push_back(q);
		for (const auto& sub : inst.sub_insts)
			for (int q : sub.qubits)
				if (q >= 0 && q < nq) result.push_back(q);
		// 去重
		std::sort(result.begin(), result.end());
		result.erase(std::unique(result.begin(), result.end()), result.end());
		return result;
	}

	// ── 閘符號 ────────────────────────────────────────────────
	static std::string GateSymbol(const CircuitInstruction& inst, int q) {
		switch (inst.type) {
		case GateType::H:           return "[H]";
		case GateType::X:           return "[X]";
		case GateType::Y:           return "[Y]";
		case GateType::Z:           return "[Z]";
		case GateType::S:           return "[S]";
		case GateType::Sdg:         return "[Sdg]";
		case GateType::T:           return "[T]";
		case GateType::Tdg:         return "[Tdg]";
		case GateType::Rx:          return "[Rx]";
		case GateType::Ry:          return "[Ry]";
		case GateType::Rz:          return "[Rz]";
		case GateType::Measure:     return "[M]";
		case GateType::Parity:      return "[Par]";
		case GateType::Label:       return "[LBL:" + inst.label + "]";
		case GateType::Goto:        return "[->:" + inst.label + "]";
		case GateType::Conditional: return "[If]";
		case GateType::SWAP:        return "[x]";
		case GateType::iSWAP:       return "[ix]";
		case GateType::SqrtSWAP:    return "[Vx]";
		case GateType::CZ:
			return "[CZ]";
		case GateType::CNOT:
			return (inst.qubits.size() > 0 && q == inst.qubits[0])
				? "[ctrl]" : "[NOT]";
		case GateType::CCNOT:
			return (inst.qubits.size() > 1 && q == inst.qubits[2])
				? "[NOT]" : "[ctrl]";
		case GateType::CSWAP:
			return (inst.qubits.size() > 0 && q == inst.qubits[0])
				? "[ctrl]" : "[x]";
		case GateType::Deutsch:
			return (inst.qubits.size() > 1 && q == inst.qubits[2])
				? "[D]" : "[ctrl]";
		default:
			return "[?]";
		}
	}

	// ══════════════════════════════════════════════════════════
	// 糾纏組輸出
	// ══════════════════════════════════════════════════════════

	static void PrintGroup(
		const EntangledGroupSnapshot& group,
		const std::map<int, std::pair<double, double>>& probabilities,
		std::ostream& out) {

		int n = static_cast<int>(group.qubit_ids.size());

		// ── 組標題 ────────────────────────────────────────
		if (n == 1) {
			out << "   Q[" << group.qubit_ids[0] << "] (獨立)\n";
		}
		else {
			out << "   {";
			for (int i = 0; i < n; i++) {
				if (i) out << ", ";
				out << "Q[" << group.qubit_ids[i] << "]";
			}
			out << "} (糾纏)\n";
		}

		// ── 狀態向量 ──────────────────────────────────────
		out << "   " << std::left << std::setw(6) << "State"
			<< std::right << std::setw(20) << "振幅"
			<< std::right << std::setw(10) << "機率"
			<< "\n";
		out << "   " << std::string(36, '-') << "\n";

		for (int idx = 0; idx < static_cast<int>(group.amplitudes.size()); idx++) {
			double re = group.amplitudes[idx].first;
			double im = group.amplitudes[idx].second;
			double prob = re * re + im * im;
			if (prob < 1e-9) continue;

			out << "   |" << StateBinary(idx, n) << ">  "
				<< std::right << std::setw(18) << FormatComplex(re, im)
				<< std::right << std::setw(8) << FormatPercent(prob)
				<< "\n";
		}

		// ── 各 qubit 的 P(0) / P(1)（獨立 qubit 才印）──
		if (n == 1) {
			int q = group.qubit_ids[0];
			auto it = probabilities.find(q);
			if (it != probabilities.end())
				out << "   P(0)=" << FormatPercent(it->second.first)
				<< "  P(1)=" << FormatPercent(it->second.second) << "\n";
		}

		// 糾纏組：逐一列出每個 qubit 的 reduced probability
		if (n > 1) {
			out << "   各 qubit 邊際機率\n";
			for (int qi = 0; qi < n; qi++) {
				int q = group.qubit_ids[qi];
				auto it = probabilities.find(q);
				if (it != probabilities.end())
					out << "     Q[" << q << "]"
					<< "  P(0)=" << FormatPercent(it->second.first)
					<< "  P(1)=" << FormatPercent(it->second.second)
					<< "\n";
			}
		}

		out << "\n";
	}

	// ══════════════════════════════════════════════════════════
	// 格式化工具
	// ══════════════════════════════════════════════════════════

	static std::string StateBinary(int idx, int n) {
		std::string s(n, '0');
		for (int i = 0; i < n; i++)
			s[n - 1 - i] = ((idx >> i) & 1) ? '1' : '0';
		return s;
	}

	static std::string FormatComplex(double re, double im) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(4) << re;
		if (im >= 0) ss << "+" << std::fixed << std::setprecision(4) << im << "i";
		else         ss << std::fixed << std::setprecision(4) << im << "i";
		return ss.str();
	}

	static std::string FormatPercent(double p) {
		std::ostringstream ss;
		ss << std::fixed << std::setprecision(1) << p * 100.0 << "%";
		return ss.str();
	}

	static void PrintDivider(char ch, int width, std::ostream& out) {
		out << " " << std::string(width, ch) << "\n";
	}

	static void PrintHeader(const CircuitExpression& circuit, std::ostream& out) {
		out << "\n";
		PrintDivider('=', 50, out);
		out << " Quantum Circuit Simulation\n";
		out << " Qubits: " << circuit.num_qubits
			<< "  Steps: " << circuit.instructions.size() << "\n";
		PrintDivider('=', 50, out);
	}

};


int main() {

	Qubit_Simulation sim = Qubit_Simulation();

	auto circuit = CircuitExpression::Parse("3, {H,0}, {CNOT,0,1}");
	//circuit.Execute(sim);
	CircuitResult result = circuit.ExecuteWithCapture(sim);
	CircuitPrinter::Print(circuit, result);
	//sim.OuputEntangledQubitSet(0);

}



/*
	Qubit_Simulation a = Qubit_Simulation(3);


	Qubit_Simulation::Matrix A;


	A.InputMatrix();
	//Qubit_Simulation::Complex::OutputComplex(Qubit_Simulation::QMath::DeterminantsOfMatrix(A));
	Qubit_Simulation::Matrix::OutputMatrix(Qubit_Simulation::QMath::GenerateQRFactorization(A).first);



a.HadamardGate(0);
a.CNOTGate(0, 1);
a.OuputEntangledQubitSet(0);
//a.HadamardGate(1);
//a.OuputEntangledQubitSet(0);

a.CNOTGate(1, 2);
a.OuputEntangledQubitSet(0);

a.CNOTGate(2, 1);
a.OuputEntangledQubitSet(0);
std::cout << std::endl;


std::cout << a.ObserverQubit(0) << std::endl;
std::cout << a.ObserverQubit(1) << std::endl;



	a.PositiveGenerateBellState();

	a.Hadamard(1);
	a.CNOT(0,1);

	//a.CCNOT(0,1,2);
	//a.CSWAP(0,1,2);
	a.OuputEntangledQubitSet(0);
	
	*/
