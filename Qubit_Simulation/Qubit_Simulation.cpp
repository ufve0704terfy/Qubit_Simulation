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

 // ── 閘的種類 ─────────────────────────────────────────
enum class GateType {
	Identity,
	// 單 qubit
	H, X, Y, Z, S, Sdg, T, Tdg,
	Rx, Ry, Rz,
	// 雙 qubit
	CNOT, CZ, SWAP, iSWAP, SqrtSWAP,
	// 三 qubit
	CCNOT, CSWAP, Deutsch,
	// 特殊
	Measure,
	Custom
};

// ── 表格中一格的內容 ──────────────────────────────────
struct GateCell {
	GateType            type = GateType::Identity;
	std::vector<double> params = {};   // Rx/Ry/Rz 的角度，Deutsch 的角度
	int                 link_qubit = -1;   // 多 qubit 閘的另一個/控制 qubit
	int                 link_qubit2 = -1;   // 三 qubit 閘的第三個 qubit
	bool                is_primary = true; // 多 qubit 閘只由 primary 觸發執行
	std::string         custom_name = "";   // 自訂閘名稱
};

// ── 每個 time step 執行後的快照 ───────────────────────
struct StepSnapshot {
	int	step;
	// set_id → {qubit列表, 狀態向量}

	std::map<int, std::pair<std::vector<int>,
	std::vector<std::pair<double, double>>>> entangled_sets;

	// qubit → 測量結果（-1 代表尚未測量）

	std::map<int, int>                           measurement;

	// qubit → {P(0), P(1)}
	std::map<int, std::pair<double, double>>      probabilities;
};

// ── 整份電路的執行結果 ────────────────────────────────

struct CircuitResult {
	int                     num_qubits;
	int                     num_steps;
	std::vector<StepSnapshot> snapshots;   // 每個 time step 一份
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

		for(int count=0;count<qubits_in_set.size();count++)
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


};

class CircuitTable {

public:

	// table[qubit][time_step]
	std::vector<std::vector<GateCell>> table;
	int num_qubits, num_steps;

	CircuitTable(int qubits, int steps)
		: num_qubits(qubits), num_steps(steps),
		table(qubits, std::vector<GateCell>(steps)) {
	}

	// ── 單 qubit 閘 ────────────────────────────────
	void SetGate(int qubit, int time, GateType type,
		std::vector<double> params = {}) {
		table[qubit][time] = { type, params };
	}

	// ── 雙 qubit 閘：ctrl 是 primary ───────────────
	void SetTwoQubitGate(int primary, int secondary,
		int time, GateType type,
		std::vector<double> params = {}) {
		table[primary][time] = { type, params, secondary, -1, true };
		table[secondary][time] = { type, params, primary,  -1, false };
	}

	// ── 三 qubit 閘：q0 是 primary ─────────────────
	void SetThreeQubitGate(int q0, int q1, int q2,
		int time, GateType type,
		std::vector<double> params = {}) {
		table[q0][time] = { type, params, q1, q2, true };
		table[q1][time] = { type, params, q0, q2, false };
		table[q2][time] = { type, params, q0, q1, false };
	}

	// ── 測量 ────────────────────────────────────────
	void SetMeasure(int qubit, int time) {
		table[qubit][time] = { GateType::Measure };
	}

	// ── 執行並回傳結果 ───────────────────────────────
	CircuitResult Execute(Qubit_Simulation& sim) const {

		CircuitResult result{ num_qubits, num_steps };

		for (int t = 0; t < num_steps; t++) {

			// 執行這個 time step 的所有閘
			for (int q = 0; q < num_qubits; q++) {
				const GateCell& cell = table[q][t];
				if (!cell.is_primary) continue;  // 多 qubit 閘只執行一次
				ApplyGate(sim, q, cell);
			}

			// 記錄這個 time step 的快照
			result.snapshots.push_back(TakeSnapshot(sim, t));
		}

		return result;
	}

private:

	void ApplyGate(Qubit_Simulation& sim,
		int qubit, const GateCell& cell) const {
		switch (cell.type) {
		case GateType::H:       sim.HadamardGate(qubit); break;
		case GateType::X:       sim.XGate(qubit);        break;
		case GateType::Y:       sim.YGate(qubit);        break;
		case GateType::Z:       sim.ZGate(qubit);        break;
		case GateType::S:       sim.SGate(qubit);        break;
		case GateType::Sdg:     sim.SdgGate(qubit);      break;
		case GateType::T:       sim.TGate(qubit);        break;
		case GateType::Tdg:     sim.TdgGate(qubit);      break;
		case GateType::Rx:      sim.RxGate(qubit, cell.params[0]); break;
		case GateType::Ry:      sim.RyGate(qubit, cell.params[0]); break;
		case GateType::Rz:      sim.RzGate(qubit, cell.params[0]); break;

		case GateType::CNOT:    sim.CNOTGate(qubit, cell.link_qubit);  break;
		case GateType::CZ:      sim.CZGate(qubit, cell.link_qubit);    break;
		case GateType::SWAP:    sim.SWAPGate(qubit, cell.link_qubit);  break;
		case GateType::iSWAP:   sim.ISWAPGate(qubit, cell.link_qubit); break;
		case GateType::SqrtSWAP:sim.SqrtSWAPGate(qubit, cell.link_qubit); break;

		case GateType::CCNOT:
			sim.CCNOTGate(qubit, cell.link_qubit, cell.link_qubit2); break;
		case GateType::CSWAP:
			sim.CSWAPGate(qubit, cell.link_qubit, cell.link_qubit2); break;
		case GateType::Deutsch:
			sim.DeutschGate(qubit, cell.link_qubit,
				cell.link_qubit2, cell.params[0]);       break;

		case GateType::Measure: sim.ObserverQubit(qubit); break;
		default: break;
		}
	}

	StepSnapshot TakeSnapshot(Qubit_Simulation& sim, int step) const {
		// 這裡需要 Qubit_Simulation 開放幾個查詢介面
		// 目前你的 OuputEntangledQubitSet 可以先用
		// 建議補充以下三個 getter：
		//   GetEntangledSets()   → 所有 set 的資料
		//   GetProbabilities(q)  → {P0, P1}
		//   GetMeasurement(q)    → 結果或 -1
		StepSnapshot snap;
		snap.step = step;
		// ... 填入資料
		return snap;
	}
};


int main() {





	Qubit_Simulation a = Qubit_Simulation(3);


	/*
	Qubit_Simulation::Matrix A;


	A.InputMatrix();
	//Qubit_Simulation::Complex::OutputComplex(Qubit_Simulation::QMath::DeterminantsOfMatrix(A));
	Qubit_Simulation::Matrix::OutputMatrix(Qubit_Simulation::QMath::GenerateQRFactorization(A).first);

	*/


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



	/*
		a.PositiveGenerateBellState();

		a.Hadamard(1);
		a.CNOT(0,1);

		//a.CCNOT(0,1,2);
		//a.CSWAP(0,1,2);
		a.OuputEntangledQubitSet(0);*/


}