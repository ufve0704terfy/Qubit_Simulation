#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <map>
#include <cmath>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>
#include <sstream>

using namespace boost::multiprecision;

static constexpr const double Double_Epsilon = 1e-9;
static constexpr const long long int Fixed_Point = (1LL << 62);
static constexpr const int Fixed_shift = 62;
static constexpr const double Root_half = 0.70710678118;
static constexpr const double Pi = 3.141592658979;
static constexpr const unsigned long long int Fixed_Pi = 14488038916154245120ULL;
static constexpr const long long int Tolerent_Round = 1LL << 10;

// ══════════════════════════════════════════════════════════════
// 快照與誤差相關結構
// ══════════════════════════════════════════════════════════════

struct EntangledGroupSnapshot {
    std::vector<int>                       qubit_ids;
    std::vector<std::pair<double, double>>  amplitudes;  // {real, imag}
};

enum class ErrorType { None, X, Y, Z };

struct ErrorRecord {
    bool      occurred = false;
    ErrorType type = ErrorType::None;
    int       affected_qubit = -1;
    double    gate_error_prob = 0.0;
};

struct StepSnapshot {
    int                                    step_index = 0;
    std::string                            step_label = "";
    std::vector<EntangledGroupSnapshot>    groups = {};
    std::map<int, std::pair<double, double>>probabilities = {};
    std::map<int, int>                     measurements = {};
    ErrorRecord  error_record;
    int          cumulative_errors = 0;
    double       theoretical_error_accum = 0.0;
};

struct SimulationConfig {
    double error_rate = 0.03;
    bool   apply_errors = false;
    bool   track_errors = true;
    double three_qubit_error_mult = 1.5;
};

struct CircuitResult {
    int num_qubits = 0;
    std::vector<StepSnapshot> snapshots = {};
    SimulationConfig config; 
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

// ══════════════════════════════════════════════════════════════
// 單一指令
// ══════════════════════════════════════════════════════════════
struct CircuitInstruction {
    GateType                        type = GateType::Identity;
    std::vector<int>                qubits = {};
    double                          param = 0.0;
    int                             expected_val = 0;
    std::string                     label = "";
    std::vector<CircuitInstruction> sub_insts = {};
};

// ══════════════════════════════════════════════════════════════
// FixedComplex 與 Complex 類別（保持原樣）
// ══════════════════════════════════════════════════════════════

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
};

// ══════════════════════════════════════════════════════════════
// Qubit_Simulation 類別（保持原樣，已包含所有功能）
// ══════════════════════════════════════════════════════════════

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
        if (Entangled_Qubit_Set[qubit]->first.size() == 1) return true;

        int q_pos = GetSituation(qubit);
        std::vector<FixedComplex>& state = Entangled_Qubit_Set[qubit]->second;
        size_t half = state.size() >> 1;

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

        for (size_t i = 1; i < qubits.size(); i++)
            if (Entangled_Qubit_Set[qubits[0]]->first != Entangled_Qubit_Set[qubits[i]]->first)
                CombineEntangledQubitSet(qubits[0], qubits[i]);

        auto& state = Entangled_Qubit_Set[qubits[0]]->second;

        std::vector<int> positions;
        for (int q : qubits)
            positions.push_back(GetSituation(q));

        long long P_even = 0, P_odd = 0;
        for (size_t idx = 0; idx < state.size(); idx++) {
            int parity = 0;
            for (int pos : positions)
                parity ^= ((idx >> pos) & 1);

            int128_t abs_val = QMath::AbsoluteValue(state[idx]);
            long long prob = static_cast<long long>((abs_val * abs_val) >> Fixed_shift);
            (parity == 0 ? P_even : P_odd) += prob;
        }

        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<long long> dist(0, Fixed_Point);
        int result = (dist(gen) <= P_even) ? 0 : 1;
        long long norm = (result == 0) ? P_even : P_odd;
        long long norm_factor = QMath::NewtonSqrt((int128_t)norm << Fixed_shift);

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
        if (IsObservered(qubit)) return { 0.0, 0.0 };
        if (Entangled_Qubit_Set.find(qubit) == Entangled_Qubit_Set.end()) return { 0.0, 0.0 };

        auto& state = Entangled_Qubit_Set[qubit]->second;
        int q_pos = GetSituation(qubit);
        size_t half = state.size() >> 1;
        double P0 = 0.0, P1 = 0.0;

        for (size_t c = 0; c < half; c++) {
            auto a0 = state[BitAdd(c, q_pos, 0ULL)];
            auto a1 = state[BitAdd(c, q_pos, 1ULL)];

            double r0 = FixedPointToDouble(a0.Real), i0 = FixedPointToDouble(a0.Imaginary);
            double r1 = FixedPointToDouble(a1.Real), i1 = FixedPointToDouble(a1.Imaginary);

            P0 += r0 * r0 + i0 * i0;
            P1 += r1 * r1 + i1 * i1;
        }

        return { P0, P1 };
    }

    std::vector<EntangledGroupSnapshot> GetEntangledGroups() {
        std::vector<EntangledGroupSnapshot> result;
        std::vector<bool> visited(Qubit_Amount, false);

        for (int q = 0; q < Qubit_Amount; q++) {
            if (visited[q] || IsObservered(q)) continue;
            if (Entangled_Qubit_Set.find(q) == Entangled_Qubit_Set.end()) continue;

            auto* ptr = Entangled_Qubit_Set[q];
            if (!ptr) continue;
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
            case GateType::X:
                std::swap(s0, s1);
                break;
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
            case GateType::Rz:
                s0 = s0 * FixedComplex(FixedPoint(c), FixedPoint(-s));
                s1 = s1 * FixedComplex(FixedPoint(c), FixedPoint(s));
                break;
            case GateType::Y:
                std::swap(s0, s1);
                tmp = s0;
                s0.Real = -tmp.Imaginary; s0.Imaginary = tmp.Real;
                tmp = s1;
                s1.Real = tmp.Imaginary; s1.Imaginary = -tmp.Real;
                break;
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
                s10.Real = -tmp.Imaginary; s10.Imaginary = tmp.Real;
                tmp = s01;
                s01.Real = -tmp.Imaginary; s01.Imaginary = tmp.Real;
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
        for (int count = 0; count < (int)qubits_in_set.size(); count++)
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
        for (int count = 0; count < (int)qubits_in_set.size(); count++)
            ExtractSeparatedQubit(qubits_in_set[count]);
    }

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
    void CNOTGate(int c, int t) { ApplyTwoQubitGate(c, t, GateType::CNOT); }
    void CZGate(int c, int t) { ApplyTwoQubitGate(c, t, GateType::CZ); }
    void SWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::SWAP); }
    void ISWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::iSWAP); }
    void SqrtSWAPGate(int a, int b) { ApplyTwoQubitGate(a, b, GateType::SqrtSWAP); }
    void CCNOTGate(int c1, int c2, int t) { ApplyThreeQubitGate(c1, c2, t, GateType::CCNOT); }
    void CSWAPGate(int c, int a, int b) { ApplyThreeQubitGate(c, a, b, GateType::CSWAP); }
    void DeutschGate(int c1, int c2, int t, double angle) {
        ApplyThreeQubitGate(c1, c2, t, GateType::Deutsch, angle);
    }
};

// ══════════════════════════════════════════════════════════════
// CircuitExpression 與執行（含誤差追蹤）
// ══════════════════════════════════════════════════════════════

class CircuitExpression {
public:
    int                             num_qubits = 0;
    std::vector<CircuitInstruction> instructions = {};

    static CircuitExpression Parse(const std::string& expr) {
        CircuitExpression result;
        std::string s;
        for (char ch : expr)
            if (!std::isspace(static_cast<unsigned char>(ch))) s += ch;
        size_t fb = s.find('{');
        if (fb == std::string::npos)
            throw std::runtime_error("CircuitExpression: missing {}");
        result.num_qubits = std::stoi(s.substr(0, fb - 1));
        for (const auto& b : SplitTopLevelBlocks(s.substr(fb)))
            result.instructions.push_back(ParseBlock(b));
        return result;
    }

    void Execute(Qubit_Simulation& sim,
        const SimulationConfig& cfg = SimulationConfig()) const {
        sim.ResetQubitSet();
        for (int count = 0; count < num_qubits; count++)
            sim.GenerateQubit();

        auto lm = BuildLabelMap();
        std::mt19937_64 rng(std::random_device{}());
        int ip = 0;
        while (ip < static_cast<int>(instructions.size())) {
            ErrorRecord dummy;
            int jump = ExecuteOne(instructions[ip], sim, lm, cfg, rng, dummy);
            ip = (jump >= 0) ? jump : ip + 1;
        }
    }

    CircuitResult ExecuteWithCapture(
        Qubit_Simulation& sim,
        const SimulationConfig& cfg = SimulationConfig()) const {

        sim.ResetQubitSet();
        for (int count = 0; count < num_qubits; count++)
            sim.GenerateQubit();

        CircuitResult result;
        result.num_qubits = num_qubits;
        result.config = cfg;

        auto lm = BuildLabelMap();
        std::mt19937_64 rng(std::random_device{}());

        int    cumulative_errors = 0;
        double theoretical_accum = 0.0;

        int ip = 0;
        while (ip < static_cast<int>(instructions.size())) {
            const auto& inst = instructions[ip];
            ErrorRecord err;
            int jump = ExecuteOne(inst, sim, lm, cfg, rng, err);

            if (cfg.track_errors && err.gate_error_prob > 0.0)
                theoretical_accum = 1.0
                - (1.0 - theoretical_accum) * (1.0 - err.gate_error_prob);
            if (cfg.track_errors && err.occurred)
                cumulative_errors++;

            StepSnapshot snap = TakeSnapshot(inst, ip, sim);
            snap.error_record = err;
            snap.cumulative_errors = cumulative_errors;
            snap.theoretical_error_accum = theoretical_accum;

            result.snapshots.push_back(snap);
            ip = (jump >= 0) ? jump : ip + 1;
        }
        return result;
    }

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

    static std::string InstructionToString(const CircuitInstruction& inst) {
        std::ostringstream ss;
        ss << "{" << GateTypeName(inst.type);
        if (!inst.label.empty())   ss << "," << inst.label;
        for (int q : inst.qubits) ss << ",Q" << q;
        if (inst.param != 0.0)
            ss << "," << std::fixed << std::setprecision(1) << inst.param;
        if (inst.type == GateType::Conditional)
            ss << ",exp=" << inst.expected_val;
        for (const auto& s : inst.sub_insts) ss << " " << InstructionToString(s);
        ss << "}";
        return ss.str();
    }

private:

    // ── 解析 ──────────────────────────────────────────────────

    static GateType ParseGateName(const std::string& name) {
        static const std::unordered_map<std::string, GateType> T = {
            {"H",GateType::H},{"X",GateType::X},{"Y",GateType::Y},
            {"Z",GateType::Z},{"S",GateType::S},{"Sdg",GateType::Sdg},
            {"T",GateType::T},{"Tdg",GateType::Tdg},
            {"Rx",GateType::Rx},{"Ry",GateType::Ry},{"Rz",GateType::Rz},
            {"CNOT",GateType::CNOT},{"CZ",GateType::CZ},
            {"SWAP",GateType::SWAP},{"iSWAP",GateType::iSWAP},
            {"SqrtSWAP",GateType::SqrtSWAP},
            {"CCNOT",GateType::CCNOT},{"CSWAP",GateType::CSWAP},
            {"Deutsch",GateType::Deutsch},
            {"M",GateType::Measure},{"Measure",GateType::Measure},
            {"Parity",GateType::Parity},
            {"If",GateType::Conditional},
            {"Label",GateType::Label},{"Goto",GateType::Goto},
        };
        auto it = T.find(name);
        return (it != T.end()) ? it->second : GateType::Identity;
    }

    static std::pair<int, int> GateSignature(GateType t) {
        switch (t) {
        case GateType::H:  case GateType::X:  case GateType::Y:
        case GateType::Z:  case GateType::S:  case GateType::Sdg:
        case GateType::T:  case GateType::Tdg:
        case GateType::Measure: return { 1,0 };
        case GateType::Rx: case GateType::Ry: case GateType::Rz:
            return { 1,1 };
        case GateType::CNOT:  case GateType::CZ:
        case GateType::SWAP:  case GateType::iSWAP:
        case GateType::SqrtSWAP: return { 2,0 };
        case GateType::CCNOT: case GateType::CSWAP: return { 3,0 };
        case GateType::Deutsch: return { 3,1 };
        case GateType::Parity:  return { -1,0 };
        default: return { 0,0 };
        }
    }

    static std::vector<std::string> SplitTopLevelBlocks(const std::string& s) {
        std::vector<std::string> res;
        int depth = 0; std::string cur;
        for (char ch : s) {
            if (ch == '{') { if (depth > 0) cur += ch; depth++; }
            else if (ch == '}') {
                depth--;
                if (depth > 0) cur += ch;
                else { if (!cur.empty()) res.push_back(cur); cur.clear(); }
            }
            else { if (depth > 0) cur += ch; }
        }
        return res;
    }

    static std::vector<std::string> SplitTokens(const std::string& s) {
        std::vector<std::string> toks; std::string cur;
        for (size_t i = 0; i <= s.size(); i++) {
            if (i == s.size() || s[i] == ',') { if (!cur.empty()) toks.push_back(cur); cur.clear(); }
            else cur += s[i];
        }
        return toks;
    }

    static CircuitInstruction ParseBlock(const std::string& block) {
        size_t fc = block.find(',');
        std::string gname = (fc == std::string::npos) ? block : block.substr(0, fc);
        CircuitInstruction inst;
        inst.type = ParseGateName(gname);
        if (fc == std::string::npos) return inst;
        std::string rest = block.substr(fc + 1);

        if (inst.type == GateType::Label || inst.type == GateType::Goto)
        {
            inst.label = rest; return inst;
        }

        if (inst.type == GateType::Conditional) {
            size_t fb = rest.find('{');
            std::string hdr = (fb == std::string::npos) ? rest : rest.substr(0, fb - 1);
            auto htok = SplitTokens(hdr);
            if (htok.size() >= 2) {
                inst.qubits.push_back(std::stoi(htok[0]));
                inst.expected_val = std::stoi(htok[1]);
            }
            if (fb != std::string::npos)
                for (const auto& s : SplitTopLevelBlocks(rest.substr(fb)))
                    inst.sub_insts.push_back(ParseBlock(s));
            return inst;
        }

        std::pair<int, int> sig = GateSignature(inst.type);
        int nq = sig.first, np = sig.second;
        if (nq == -1) {
            for (const auto& t : SplitTokens(rest)) inst.qubits.push_back(std::stoi(t));
        }
        else {
            auto toks = SplitTokens(rest);
            for (int i = 0; i < nq && i < (int)toks.size(); i++) inst.qubits.push_back(std::stoi(toks[i]));
            if (np > 0 && nq < (int)toks.size()) inst.param = std::stod(toks[nq]);
        }
        return inst;
    }

    // ── 誤差模型 ──────────────────────────────────────────────

    static bool IsTwoQubitGate(GateType t) {
        return t == GateType::CNOT || t == GateType::CZ || t == GateType::SWAP ||
            t == GateType::iSWAP || t == GateType::SqrtSWAP;
    }
    static bool IsThreeQubitGate(GateType t) {
        return t == GateType::CCNOT || t == GateType::CSWAP || t == GateType::Deutsch;
    }

    static ErrorType ApplyRandomPauli(Qubit_Simulation& sim, int q, std::mt19937_64& rng) {
        std::uniform_int_distribution<int> d(0, 2);
        switch (d(rng)) {
        case 0: sim.ApplySingleQubitGate(q, GateType::X); return ErrorType::X;
        case 1: sim.ApplySingleQubitGate(q, GateType::Z); return ErrorType::Z;
        default:sim.ApplySingleQubitGate(q, GateType::Y); return ErrorType::Y;
        }
    }

    static void HandleError(const CircuitInstruction& inst, Qubit_Simulation& sim,
        const SimulationConfig& cfg, std::mt19937_64& rng,
        ErrorRecord& rec) {
        bool is2 = IsTwoQubitGate(inst.type);
        bool is3 = IsThreeQubitGate(inst.type);
        if (!is2 && !is3) return;

        double p = cfg.error_rate;
        if (is3) p = std::min(1.0, p * cfg.three_qubit_error_mult);
        rec.gate_error_prob = p;

        if (!cfg.apply_errors) return;

        std::uniform_real_distribution<double> roll(0.0, 1.0);
        if (roll(rng) < p) {
            std::uniform_int_distribution<int> pq(0, (int)inst.qubits.size() - 1);
            int tgt = inst.qubits[pq(rng)];
            rec.occurred = true;
            rec.affected_qubit = tgt;
            rec.type = ApplyRandomPauli(sim, tgt, rng);
        }
    }

    // ── 執行 ──────────────────────────────────────────────────

    std::unordered_map<std::string, int> BuildLabelMap() const {
        std::unordered_map<std::string, int> m;
        for (int i = 0; i < (int)instructions.size(); i++)
            if (instructions[i].type == GateType::Label)
                m[instructions[i].label] = i;
        return m;
    }

    static int ExecuteOne(const CircuitInstruction& inst, Qubit_Simulation& sim,
        const std::unordered_map<std::string, int>& lm,
        const SimulationConfig& cfg, std::mt19937_64& rng,
        ErrorRecord& err) {
        switch (inst.type) {
        case GateType::H:  case GateType::X:  case GateType::Y:
        case GateType::Z:  case GateType::S:  case GateType::Sdg:
        case GateType::T:  case GateType::Tdg:
            sim.ApplySingleQubitGate(inst.qubits[0], inst.type); return -1;
        case GateType::Rx: case GateType::Ry: case GateType::Rz:
            sim.ApplySingleQubitGate(inst.qubits[0], inst.type, inst.param); return -1;
        case GateType::CNOT:  case GateType::CZ:
        case GateType::SWAP:  case GateType::iSWAP: case GateType::SqrtSWAP:
            sim.ApplyTwoQubitGate(inst.qubits[0], inst.qubits[1], inst.type);
            HandleError(inst, sim, cfg, rng, err); return -1;
        case GateType::CCNOT: case GateType::CSWAP:
            sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1], inst.qubits[2], inst.type);
            HandleError(inst, sim, cfg, rng, err); return -1;
        case GateType::Deutsch:
            sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1], inst.qubits[2], inst.type, inst.param);
            HandleError(inst, sim, cfg, rng, err); return -1;
        case GateType::Measure:  sim.ObserverQubit(inst.qubits[0]); return -1;
        case GateType::Parity:   sim.MeasureParity(inst.qubits);    return -1;
        case GateType::Label:    return -1;
        case GateType::Goto: {
            auto it = lm.find(inst.label);
            return (it != lm.end()) ? it->second : -1;
        }
        case GateType::Conditional: {
            int meas = sim.ObserverQubit(inst.qubits[0]);
            if (meas == inst.expected_val)
                for (const auto& sub : inst.sub_insts) {
                    ErrorRecord se;
                    int jump = ExecuteOne(sub, sim, lm, cfg, rng, se);
                    if (se.occurred && !err.occurred) err = se;
                    err.gate_error_prob += se.gate_error_prob;
                    if (jump >= 0) return jump;
                }
            return -1;
        }
        default: return -1;
        }
    }

    static StepSnapshot TakeSnapshot(const CircuitInstruction& inst,
        int step, Qubit_Simulation& sim) {
        StepSnapshot snap;
        snap.step_index = step;
        snap.step_label = InstructionToString(inst);
        snap.groups = sim.GetEntangledGroups();
        int nq = sim.GetQubitAmount();
        for (int q = 0; q < nq; q++) {
            int meas = sim.GetMeasurementResult(q);
            snap.measurements[q] = meas;
            if (meas == -1) snap.probabilities[q] = sim.GetProbabilities(q);
        }
        return snap;
    }
};

// ══════════════════════════════════════════════════════════════
// CircuitPrinter（合併版本：完整列印 + 誤差追蹤）
// ══════════════════════════════════════════════════════════════

class CircuitPrinter {
public:

    static void Print(const CircuitExpression& circuit,
        const CircuitResult& result,
        std::ostream& out = std::cout) {
        PrintHeader(circuit, result, out);
        PrintCircuitDiagram(circuit, out);
        out << Row('=', W) << "\n";
        for (const auto& snap : result.snapshots)
            PrintSnapshot(snap, result.num_qubits, result, out);

        if (result.config.track_errors)
            PrintErrorSummary(result, out);

        out << Row('=', W) << "\n";
    }

    static void PrintDiagram(const CircuitExpression& circuit,
        std::ostream& out = std::cout) {
        PrintHeader(circuit, CircuitResult(), out);
        PrintCircuitDiagram(circuit, out);
    }

private:

    static const int W = 70;

    // ── 格式化工具 ──────────────────────────────────────────────

    static std::string Pad(const std::string& s, int width) {
        if ((int)s.size() >= width) return s;
        return s + std::string(width - s.size(), ' ');
    }
    static std::string PadR(const std::string& s, int width) {
        if ((int)s.size() >= width) return s;
        return std::string(width - s.size(), ' ') + s;
    }
    static std::string Row(char ch, int n) { return std::string(n, ch); }

    static std::string ProbBar(double p1, int width = 20) {
        int fill = static_cast<int>(std::round(p1 * width));
        fill = std::max(0, std::min(width, fill));
        return "[" + std::string(fill, '#') + std::string(width - fill, ' ') + "]";
    }

    static std::string FmtPct(double p) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(2) << p * 100.0 << "%";
        return ss.str();
    }

    static std::string FmtComplex(double re, double im) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << (re >= 0 ? " " : "") << re;
        ss << (im >= 0 ? "+" : "-");
        ss << std::abs(im) << "i";
        return ss.str();
    }

    static std::string StateBin(int idx, int n) {
        std::string s(n, '0');
        for (int i = 0; i < n; i++) s[n - 1 - i] = ((idx >> i) & 1) ? '1' : '0';
        return s;
    }

    static std::string ErrorTypeName(ErrorType t) {
        switch (t) {
        case ErrorType::X: return "X-error";
        case ErrorType::Y: return "Y-error";
        case ErrorType::Z: return "Z-error";
        default:           return "none";
        }
    }

    // ── 標題 ──────────────────────────────────────────────────

    static void PrintHeader(const CircuitExpression& circuit,
        const CircuitResult& result,
        std::ostream& out) {
        out << "\n" << Row('=', W) << "\n";
        out << " Quantum Circuit Simulation\n";
        out << " Qubits : " << circuit.num_qubits
            << "   Steps : " << circuit.instructions.size() << "\n";
        if (!result.snapshots.empty()) {
            out << " Error  : rate=" << FmtPct(result.config.error_rate)
                << "  tracking=" << (result.config.track_errors ? "on" : "off") << "\n";
        }
        out << Row('=', W) << "\n";
    }

    // ── 電路圖 ────────────────────────────────────────────────

    static void PrintCircuitDiagram(const CircuitExpression& circuit,
        std::ostream& out) {
        int nq = circuit.num_qubits;
        const int CW = 9;

        struct Slot {
            std::vector<bool>               occ;
            std::vector<CircuitInstruction> cell;
            Slot(int n) : occ(n, false), cell(n) {}
        };
        std::vector<Slot> slots;

        for (const auto& inst : circuit.instructions) {
            std::vector<int> tgts = CollectQubits(inst, nq);
            bool standalone = (inst.type == GateType::Label ||
                inst.type == GateType::Goto ||
                inst.type == GateType::Conditional);
            int si = -1;

            if (!standalone && !tgts.empty()) {
                int max_s = -1;
                // 找出這些目標 qubit 最後被佔用的 Slot (時間點)
                for (int q : tgts) {
                    for (int s = (int)slots.size() - 1; s >= 0; s--) {
                        if (slots[s].occ[q]) {
                            max_s = std::max(max_s, s);
                            break; // 找到該 qubit 最近被佔用的 slot 即可
                        }
                    }
                }
                // 必須排在最後被佔用之後的下一個 Slot，確保因果順序
                si = max_s + 1;
            }

            // 若超出現有 Slot 數量或不需要排入時間軸，則新增 Slot
            if (si == -1 || si == (int)slots.size()) {
                slots.push_back(Slot(nq));
                if (si == -1) si = (int)slots.size() - 1;
            }

            for (int q : tgts) { slots[si].occ[q] = true; slots[si].cell[q] = inst; }
        }

        // ... 底下印出 T0, T1 以及 Q[x] 的邏輯維持原樣不變 ...
        out << "\n";
        out << Pad("", 7);
        for (int s = 0; s < (int)slots.size(); s++)
            out << Pad("T" + std::to_string(s), CW);
        out << "\n";

        for (int q = 0; q < nq; q++) {
            out << Pad("Q[" + std::to_string(q) + "]", 7);
            for (const auto& sl : slots)
                out << Pad(sl.occ[q] ? GateSym(sl.cell[q], q) : "-------", CW);
            out << "\n";
        }
        out << "\n";
    }

    static std::vector<int> CollectQubits(const CircuitInstruction& inst, int nq) {
        std::vector<int> r;
        for (int q : inst.qubits) if (q >= 0 && q < nq) r.push_back(q);
        for (const auto& s : inst.sub_insts)
            for (int q : s.qubits) if (q >= 0 && q < nq) r.push_back(q);
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        return r;
    }

    static std::string GateSym(const CircuitInstruction& inst, int q) {
        switch (inst.type) {
        case GateType::H:      return "[H]";
        case GateType::X:      return "[X]";
        case GateType::Y:      return "[Y]";
        case GateType::Z:      return "[Z]";
        case GateType::S:      return "[S]";
        case GateType::Sdg:    return "[Sdg]";
        case GateType::T:      return "[T]";
        case GateType::Tdg:    return "[Tdg]";
        case GateType::Rx:     return "[Rx]";
        case GateType::Ry:     return "[Ry]";
        case GateType::Rz:     return "[Rz]";
        case GateType::Measure:return "[M]";
        case GateType::Parity: return "[Par]";
        case GateType::Label:  return "[LBL]";
        case GateType::Goto:   return "[GO]";
        case GateType::Conditional: return "[If]";
        case GateType::SWAP:   return "[x]";
        case GateType::iSWAP:  return "[ix]";
        case GateType::SqrtSWAP: return "[Vx]";
        case GateType::CZ:     return "[CZ]";
        case GateType::CNOT:
            return (!inst.qubits.empty() && q == inst.qubits[0]) ? "[ctrl]" : "[NOT]";
        case GateType::CCNOT:
            return (inst.qubits.size() > 2 && q == inst.qubits[2]) ? "[NOT]" : "[ctrl]";
        case GateType::CSWAP:
            return (!inst.qubits.empty() && q == inst.qubits[0]) ? "[ctrl]" : "[x]";
        case GateType::Deutsch:
            return (inst.qubits.size() > 2 && q == inst.qubits[2]) ? "[D]" : "[ctrl]";
        default: return "[?]";
        }
    }

    // ── 快照與糾纏組 ──────────────────────────────────────────

    static void PrintSnapshot(const StepSnapshot& snap,
        int                   num_qubits,
        const CircuitResult& result,
        std::ostream& out) {

        out << Row('-', W) << "\n";
        out << Pad("Step " + std::to_string(snap.step_index + 1), 12)
            << "| " << snap.step_label << "\n";
        out << Row('-', W) << "\n";

        // 誤差記錄（如果 track_errors 打開）
        if (result.snapshots[0].error_record.gate_error_prob > 0.0) {
            out << " [Error] "
                << "gate_prob=" << FmtPct(snap.error_record.gate_error_prob)
                << "  accum_theory=" << FmtPct(snap.theoretical_error_accum)
                << "  total_hits=" << snap.cumulative_errors;
            if (snap.error_record.occurred)
                out << "  *** "
                << ErrorTypeName(snap.error_record.type)
                << " on Q" << snap.error_record.affected_qubit << " ***";
            out << "\n\n";
        }

        // 測量結果
        bool any_meas = false;
        for (int q = 0; q < num_qubits; q++) {
            auto it = snap.measurements.find(q);
            if (it != snap.measurements.end() && it->second != -1) {
                if (!any_meas) {
                    out << " Measurement Results\n";
                    out << " " << Pad("Qubit", 8) << "Result\n";
                    out << " " << Row('-', 16) << "\n";
                    any_meas = true;
                }
                out << " " << Pad("Q[" + std::to_string(q) + "]", 8) << it->second << "\n";
            }
        }

        // 量子態
        if (!snap.groups.empty()) {
            if (any_meas) out << "\n";
            out << " Quantum State\n";
            for (const auto& g : snap.groups)
                PrintGroup(g, snap.probabilities, out);
        }
    }

    static void PrintGroup(
        const EntangledGroupSnapshot& g,
        const std::map<int, std::pair<double, double>>& probs,
        std::ostream& out) {

        int n = (int)g.qubit_ids.size();

        // 組標題
        std::string title = (n == 1)
            ? "Q[" + std::to_string(g.qubit_ids[0]) + "] (independent)"
            : "{" + [&] {
            std::string s;
            for (int i = 0; i < n; i++) { if (i)s += ","; s += "Q[" + std::to_string(g.qubit_ids[i]) + "]"; }
            return s;
            }() + "} (entangled)";

        out << "  " << title << "\n";

        // 欄標
        out << "  " << Pad("State", n + 4)
            << Pad("Amplitude", 20)
            << PadR("Prob", 8) << "\n";
        out << "  " << Row('-', n + 4 + 20 + 8) << "\n";

        // 各 state
        for (int idx = 0; idx < (int)g.amplitudes.size(); idx++) {
            double re = g.amplitudes[idx].first, im = g.amplitudes[idx].second;
            double p = re * re + im * im;
            if (p < 1e-9) continue;
            out << "  |" << StateBin(idx, n) << ">  "
                << Pad(FmtComplex(re, im), 20)
                << PadR(FmtPct(p), 8)
                << "\n";
        }

        // 邊際機率
        if (n == 1) {
            int q = g.qubit_ids[0];
            auto it = probs.find(q);
            if (it != probs.end())
                out << "  " << Pad("", n + 4)
                << "P(0)=" << FmtPct(it->second.first)
                << "  P(1)=" << FmtPct(it->second.second) << "\n";
        }
        else {
            out << "  Marginal probabilities:\n";
            for (int qi = 0; qi < n; qi++) {
                int q = g.qubit_ids[qi];
                auto it = probs.find(q);
                if (it != probs.end())
                    out << "    " << Pad("Q[" + std::to_string(q) + "]", 6)
                    << "P(0)=" << FmtPct(it->second.first)
                    << "  P(1)=" << FmtPct(it->second.second) << "\n";
            }
        }
        out << "\n";
    }

    static void PrintErrorSummary(const CircuitResult& result, std::ostream& out) {
        out << Row('-', W) << "\n";
        out << " Error Summary\n";
        out << Row('-', W) << "\n";

        if (result.snapshots.empty()) return;

        int total_errors = result.snapshots.back().cumulative_errors;
        double final_accum = result.snapshots.back().theoretical_error_accum;

        out << " Total error occurrences: " << total_errors << "\n";
        out << " Accumulated error prob:  " << FmtPct(final_accum) << "\n";
        out << "\n";
    }
};

// ══════════════════════════════════════════════════════════════
// Console 控制台類別 (更新版：完整指令手冊)
// ══════════════════════════════════════════════════════════════

class Console {
private:
    Qubit_Simulation sim;
    SimulationConfig config;
    CircuitExpression last_circuit;
    CircuitResult last_result;
    bool has_result = false;

    // 輔助函式：去除字串前後空白
    std::string Trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        size_t end = s.find_last_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        return s.substr(start, end - start + 1);
    }

    void PrintIntro() {
        std::cout << "======================================================\n";
        std::cout << "         歡迎使用量子電路模擬器 (Quantum Console)       \n";
        std::cout << "======================================================\n";
        std::cout << "初始狀態已載入。輸入 'help' 可以查看基本指令與操作方式。\n\n";
    }

    void CmdHelp() {
        std::cout << "\n[可用指令列表]\n";
        std::cout << "  help           : 顯示此指令列表\n";
        std::cout << "  manual         : 顯示量子電路詳細操作手冊、所有支援的閘與範例\n";
        std::cout << "  config         : 設定參數 (如: config error_rate 0.05)\n";
        std::cout << "  run \"<expr>\"   : 執行電路，表達式必須用雙引號包裝\n";
        std::cout << "                   (例如: run \"3, {H,0}, {CNOT,0,1}\")\n";
        std::cout << "  show           : 顯示最後一次執行的電路圖與狀態快照\n";
        std::cout << "  exit / quit    : 離開控制台\n\n";
    }

    void CmdManual() {
        std::cout << "\n================ 量子電路操作手冊 ================\n";
        std::cout << "[基本架構與語法]\n";
        std::cout << "本模擬器用於模擬量子位元 (Qubits) 的演化[cite: 1]。\n";
        std::cout << "指令語法必須用雙引號 \"\" 包起來，格式如下：\n";
        std::cout << "  \"<總量子位元數>, {量子閘, 目標1, 目標2...}, {量子閘...}\"\n\n";

        std::cout << "[所有可用的量子電路組件 (GateTypes)]\n";
        std::cout << "1. 基本單量子位元閘 (需 1 個目標位元)[cite: 1]:\n";
        std::cout << "   - H    : Hadamard 閘 (創造均勻疊加態)\n";
        std::cout << "   - X    : Pauli-X 閘 (量子 NOT 閘，0變1，1變0)\n";
        std::cout << "   - Y    : Pauli-Y 閘\n";
        std::cout << "   - Z    : Pauli-Z 閘 (相位反轉)\n";
        std::cout << "   - S    : S 閘 (繞 Z 軸旋轉 90 度)\n";
        std::cout << "   - Sdg  : S-dagger 閘 (S 閘的共軛轉置)\n";
        std::cout << "   - T    : T 閘 (繞 Z 軸旋轉 45 度)\n";
        std::cout << "   - Tdg  : T-dagger 閘 (T 閘的共軛轉置)\n\n";

        std::cout << "2. 參數化單量子位元閘 (需 1 個目標位元 + 1 個角度參數)[cite: 1]:\n";
        std::cout << "   (角度單位為度數 degree)\n";
        std::cout << "   - Rx   : 繞 X 軸旋轉 (例如: {Rx,0,90})\n";
        std::cout << "   - Ry   : 繞 Y 軸旋轉\n";
        std::cout << "   - Rz   : 繞 Z 軸旋轉\n\n";

        std::cout << "3. 雙量子位元閘 (需 2 個目標位元)[cite: 1]:\n";
        std::cout << "   - CNOT     : 控制反轉閘 (Control-NOT，參數順序: 控制位元, 目標位元)\n";
        std::cout << "   - CZ       : 控制相位閘 (Control-Z)\n";
        std::cout << "   - SWAP     : 交換閘 (交換兩個位元狀態)\n";
        std::cout << "   - iSWAP    : 虛數交換閘\n";
        std::cout << "   - SqrtSWAP : 平方根交換閘\n\n";

        std::cout << "4. 三量子位元閘 (需 3 個目標位元)[cite: 1]:\n";
        std::cout << "   - CCNOT    : Toffoli 閘 (雙控制 NOT，順序: 控1, 控2, 目標)\n";
        std::cout << "   - CSWAP    : Fredkin 閘 (控制交換閘)\n";
        std::cout << "   - Deutsch  : Deutsch 閘 (需額外 1 個角度參數，例如: {Deutsch,0,1,2,45})\n\n";

        std::cout << "5. 測量與控制流 (Measurement & Flow)[cite: 1]:\n";
        std::cout << "   - M / Measure: 測量單一量子位元並坍縮 (例如: {M,0})\n";
        std::cout << "   - Parity     : 測量多個量子位元的宇稱 (例如: {Parity,0,1,2})\n";
        std::cout << "   - Label      : 設定跳轉標籤 (例如: {Label,loop_start})\n";
        std::cout << "   - Goto       : 跳轉至特定標籤 (例如: {Goto,loop_start})\n";
        std::cout << "   - If         : 條件分支，若測量等於預期則執行內部指令。\n";
        std::cout << "                  格式: {If, 目標位元, 預期值, {指令...}}\n";
        std::cout << "                  (例如: {If,0,1,{X,1}}，若 Q0=1 則對 Q1 做 X 閘)\n\n";

        std::cout << "[經典示範 (Examples)]\n";
        std::cout << "▶ 範例 1：建立兩量子位元的貝爾態 (Bell State，最大糾纏態)\n";
        std::cout << "  輸入: run \"2, {H,0}, {CNOT,0,1}\"\n";
        std::cout << "  說明: 準備 2 個位元，先對 Q0 用 H 閘變成 0 與 1 的疊加態，再用 CNOT 閘將\n";
        std::cout << "        Q0 與 Q1 綁定。最終狀態為 |00> 與 |11> 各 50% 機率。\n\n";

        std::cout << "▶ 範例 2：量子遙傳 (Quantum Teleportation) 簡化結構與動態控制\n";
        std::cout << "  輸入: run \"3, {H,1}, {CNOT,1,2}, {CNOT,0,1}, {H,0}, {M,0}, {M,1}, {If,1,1,{X,2}}, {If,0,1,{Z,2}}\"\n";
        std::cout << "  說明: 此電路將未知的量子態從 Q0 傳遞到 Q2。過程包含建立糾纏(Q1,Q2)、\n";
        std::cout << "        貝爾測量(Q0,Q1)，以及透過 'If' 條件控制對 Q2 進行 X 或 Z 閘修正。\n";
        std::cout << "==================================================\n\n";
    }

    void CmdConfig(const std::string& args) {
        if (args.empty()) {
            std::cout << "\n[當前設定]\n";
            std::cout << "  Error Rate (error_rate)    : " << config.error_rate << "\n";
            std::cout << "  Apply Errors (apply_errors): " << (config.apply_errors ? "true" : "false") << "\n";
            std::cout << "  Track Errors (track_errors): " << (config.track_errors ? "true" : "false") << "\n\n";
            return;
        }

        std::stringstream ss(args);
        std::string key;
        ss >> key;

        if (key == "error_rate") {
            double val;
            if (ss >> val) { config.error_rate = val; std::cout << "已設定 error_rate = " << val << "\n"; }
        }
        else if (key == "apply_errors") {
            std::string val;
            if (ss >> val) { config.apply_errors = (val == "true" || val == "1"); std::cout << "已設定 apply_errors = " << config.apply_errors << "\n"; }
        }
        else if (key == "track_errors") {
            std::string val;
            if (ss >> val) { config.track_errors = (val == "true" || val == "1"); std::cout << "已設定 track_errors = " << config.track_errors << "\n"; }
        }
        else {
            std::cout << "未知的設定參數: " << key << "\n";
        }
    }

    void CmdRun(const std::string& args) {
        std::string expr = Trim(args);

        if (expr.empty()) {
            std::cout << "錯誤：請提供電路表達式。\n";
            std::cout << "提示：表達式必須使用雙引號 \"\" 包裝。\n";
            std::cout << "例如: run \"2, {H,0}, {CNOT,0,1}\"\n";
            return;
        }

        // 檢查開頭和結尾是否都是雙引號
        if (expr.front() == '"' && expr.back() == '"' && expr.length() >= 2) {
            // 擷取雙引號內部的字串
            std::string parsed_expr = expr.substr(1, expr.length() - 2);

            try {
                last_circuit = CircuitExpression::Parse(parsed_expr);
                std::cout << "解析成功！開始執行量子模擬...\n";
                last_result = last_circuit.ExecuteWithCapture(sim, config);
                has_result = true;
                std::cout << "模擬完成。輸入 'show' 來查看結果。\n";
            }
            catch (const std::exception& e) {
                std::cout << "解析或執行時發生錯誤: " << e.what() << "\n";
            }
        }
        else {
            // 如果沒有用雙引號包好，給予錯誤提示
            std::cout << "錯誤：語法不正確！電路表達式必須使用雙引號 \"\" 完整包裝。\n";
            std::cout << "正確範例: run \"3, {H,0}, {CNOT,0,1}\"\n";
            std::cout << "您的輸入: run " << args << "\n";
        }
    }

    void CmdShow() {
        if (!has_result) {
            std::cout << "目前沒有可顯示的執行結果。請先使用 'run' 指令執行電路。\n";
            return;
        }
        CircuitPrinter::Print(last_circuit, last_result, std::cout);
    }

public:
    Console() {
        // 初始化預設 config
        config.track_errors = true;
        config.apply_errors = false;
        config.error_rate = 0.03;
    }

    void Start() {
        PrintIntro();
        std::string input;

        while (true) {
            std::cout << "QConsole> ";
            if (!std::getline(std::cin, input)) break;

            input = Trim(input);
            if (input.empty()) continue;

            size_t space_pos = input.find(' ');
            std::string cmd = input.substr(0, space_pos);
            std::string args = (space_pos != std::string::npos) ? Trim(input.substr(space_pos + 1)) : "";

            if (cmd == "exit" || cmd == "quit") {
                std::cout << "正在結束控制台...\n";
                break;
            }
            else if (cmd == "help") {
                CmdHelp();
            }
            else if (cmd == "manual") {
                CmdManual();
            }
            else if (cmd == "config") {
                CmdConfig(args);
            }
            else if (cmd == "run") {
                CmdRun(args);
            }
            else if (cmd == "show") {
                CmdShow();
            }
            else {
                std::cout << "未知指令: '" << cmd << "'。輸入 'help' 以查看可用指令。\n";
            }
        }
    }
};

int main() {
    // 啟動互動式控制台
    Console console;
    console.Start();

    return 0;
}