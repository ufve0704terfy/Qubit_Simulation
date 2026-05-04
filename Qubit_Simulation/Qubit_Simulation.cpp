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
static constexpr const long long int Fixed_Epsilon = 500000LL;

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
    std::vector<double> qubit_error_accum = {};
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
    Conditional, Label, Goto,
    // 重置
    Reset,
    // 初始化特定 Qubit
    InitQubit
};

// ══════════════════════════════════════════════════════════════
// 單一指令
// ══════════════════════════════════════════════════════════════
struct CircuitInstruction {
    GateType                        type = GateType::Identity;
    std::vector<int>                qubits = {};
    double                          param = 0.0;
    int                             expected_val = 1;  // 預設為 1
    bool                            condition_is_measure = false; // 紀錄條件是否包含測量動作
    bool                            condition_is_parity = false; // 紀錄條件是否為宇稱測量
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

    static inline bool IsZero(const FixedComplex& C) {
        return (AbsoluteValue(C.Real) < Fixed_Epsilon) && (AbsoluteValue(C.Imaginary) < Fixed_Epsilon);
    }

    static inline long long int AbsoluteValue(const long long int& D) {
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

        if (Entangled_Qubit_Set[situation1]->second.empty() || Entangled_Qubit_Set[situation2]->second.empty()) {
            return;
        }

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
                if (!QMath::IsZero(First_Qubit_Set[x]) && !QMath::IsZero(Second_Qubit_Set[y]))
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

    void SplitEntangledGroup(const int root_situation) {
        if (Entangled_Qubit_Set.find(root_situation) == Entangled_Qubit_Set.end()) return;
        auto* group_ptr = Entangled_Qubit_Set[root_situation];
        if (!group_ptr) return;

        std::vector<int> current_qubits = group_ptr->first;
        int N = current_qubits.size();
        if (N <= 1) return;

        const std::vector<FixedComplex>& state = group_ptr->second;

        // 1. 找到 reference index (找出任何一個不為0的基礎態當作基準點)
        size_t ref_idx = 0;
        for (size_t idx = 0; idx < state.size(); idx++) {
            if (!QMath::IsZero(state[idx])) {
                ref_idx = idx;
                break;
            }
        }

        // 定義一個 Lambda 函式，利用 XOR 將基準態平移到 000...0，方便數學因式分解
        auto Amp = [&](size_t x) -> FixedComplex {
            return state[x ^ ref_idx];
            };

        FixedComplex c_0 = Amp(0);

        // 2. 初始化並查集 (DSU)，一開始假設所有 Qubit 都是獨立的
        std::vector<int> parent(N);
        for (int i = 0; i < N; i++) parent[i] = i;

        auto find = [&](int i) {
            int root = i;
            while (root != parent[root]) root = parent[root];
            int curr = i;
            while (curr != root) {
                int nxt = parent[curr];
                parent[curr] = root;
                curr = nxt;
            }
            return root;
            };

        auto merge = [&](int i, int j) {
            int root_i = find(i);
            int root_j = find(j);
            if (root_i != root_j) {
                parent[root_i] = root_j;
            }
            };

        // 3. 測試所有的狀態分佈，檢查張量積獨立性
        for (size_t x = 1; x < (1ULL << N); x++) {
            // 找出狀態 x 觸發了哪些目前假定的獨立 Component
            std::vector<int> comps;
            for (int i = 0; i < N; i++) {
                if ((x >> i) & 1) {
                    int c = find(i);
                    if (std::find(comps.begin(), comps.end(), c) == comps.end()) {
                        comps.push_back(c);
                    }
                }
            }

            // 如果該狀態跨越了兩個以上的 Component，檢查它們是否真的可以數學分解
            if (comps.size() > 1) {
                int comp_A = comps[0];
                size_t x_A = 0;
                for (int i = 0; i < N; i++) {
                    if (((x >> i) & 1) && find(i) == comp_A) {
                        x_A |= (1ULL << i);
                    }
                }
                size_t x_B = x ^ x_A; // 剩下的部份

                FixedComplex c_x = Amp(x);
                FixedComplex c_A = Amp(x_A);
                FixedComplex c_B = Amp(x_B);

                // 核心檢驗：如果獨立，必定滿足 c_x * c_0 == c_A * c_B
                FixedComplex lhs = c_x * c_0;
                FixedComplex rhs = c_A * c_B;

                FixedComplex diff = lhs - rhs;
                if (QMath::AbsoluteValue(diff) > Tolerent_Round) {
                    // 發生糾纏！因式分解失敗，必須把這些 Component 黏合在一起
                    for (size_t i = 1; i < comps.size(); i++) {
                        merge(comps[0], comps[i]);
                    }
                }
            }
        }

        // 4. 計算最終到底切出了幾個獨立的群組
        std::vector<int> comp_id(N, -1);
        int num_comps = 0;
        std::unordered_map<int, int> root_to_comp;
        for (int i = 0; i < N; i++) {
            int r = find(i);
            if (root_to_comp.find(r) == root_to_comp.end()) {
                root_to_comp[r] = num_comps++;
            }
            comp_id[i] = root_to_comp[r];
        }

        // 如果全部都在同一個連通塊，代表不可分離，直接結束
        if (num_comps <= 1) return;

        // 5. 將每個獨立群組的振幅各自萃取出來
        struct ComponentData {
            std::vector<int> names;
            std::vector<FixedComplex> state;
        };
        std::vector<ComponentData> extracted(num_comps);

        for (int k = 0; k < num_comps; k++) {
            std::vector<int> old_bit_pos;
            for (int count = 0; count < N; count++) {
                int b = N - 1 - count;
                if (comp_id[b] == k) {
                    extracted[k].names.push_back(current_qubits[count]);
                    old_bit_pos.push_back(b);
                }
            }

            int M = extracted[k].names.size();
            extracted[k].state.resize(1ULL << M, FixedComplex());

            int128_t norm_sq = 0;
            for (size_t local_idx = 0; local_idx < (1ULL << M); local_idx++) {
                size_t full_idx = ref_idx;
                for (int c = 0; c < M; c++) {
                    int new_b = M - 1 - c;
                    int old_b = old_bit_pos[c];
                    // 將 local_idx 映射回完整的大矩陣索引中
                    if ((local_idx & (1ULL << new_b)) != 0) {
                        full_idx |= (1ULL << old_b);
                    }
                    else {
                        full_idx &= ~(1ULL << old_b);
                    }
                }
                FixedComplex amp = state[full_idx];
                extracted[k].state[local_idx] = amp;

                int128_t abs_val = QMath::AbsoluteValue(amp);
                norm_sq += (abs_val * abs_val) >> Fixed_shift;
            }

            // 重新正規化 (Normalize)
            long long norm_factor = QMath::NewtonSqrt(norm_sq << Fixed_shift);
            for (auto& amp : extracted[k].state) {
                if (norm_factor > 0) amp = amp / norm_factor;
            }
        }

        // 6. 將分離好的群組正式寫回模擬器系統中
        for (int k = 0; k < num_comps; k++) {
            int root = extracted[k].names[0];
            Entangled_Qubit_Set_Pointer[root] = { extracted[k].names, extracted[k].state };
            for (int q : extracted[k].names) {
                Entangled_Qubit_Set[q] = &Entangled_Qubit_Set_Pointer[root];
            }
        }
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

    // 將所有 Qubit 重置回 |0> 初始狀態，保留 Qubit 數量不變
    void ReinitializeAllQubits() {
        int n = Qubit_Amount;
        Entangled_Qubit_Set.clear();
        Qubit_Set_Observation.clear();
        Entangled_Qubit_Set_Pointer.clear();
        for (int count = 0; count < n; count++) {
            Entangled_Qubit_Set_Pointer[count] = { {count}, {FixedComplex(Fixed_Point), FixedComplex()} };
            Qubit_Set_Observation.push_back({ 0, -1 });
            Entangled_Qubit_Set[count] = &Entangled_Qubit_Set_Pointer[count];
        }
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
        if (Probability0 <= Fixed_Epsilon) {
            Qubit_Set_Observation[situation] = { 1,1 };
            Probability1 = Fixed_Point;
        }
        else if (Probability0 >= Fixed_Point - Fixed_Epsilon) {
            Qubit_Set_Observation[situation] = { 1,0 };
            Probability0 = Fixed_Point;
        }
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

        std::vector<FixedComplex> Qubit_Set_Value = Entangled_Qubit_Set[situation]->second;

        size_t new_size = Qubit_Set_Value.size() >> 1;
        if (new_size == 0) new_size = 1;
        std::vector<FixedComplex> New_Qubit_Set_Value(new_size, FixedComplex());

        std::vector<int> Qubit_Set = Entangled_Qubit_Set[situation]->first;
        int Qubit_Situation = GetSituation(situation);

        // 2. 迴圈條件改成直接使用 new_size，不要再減去 (Qubit_Set_Value.size() == 2)
        for (unsigned long long int count = 0; count < new_size; count++) {
            if (!QMath::IsZero(Qubit_Set_Value[BitAdd(count, Qubit_Situation, State)])) {
                New_Qubit_Set_Value[count] = Qubit_Set_Value[BitAdd(count, Qubit_Situation, State)] / QMath::NewtonSqrt((int128_t)Probability << Fixed_shift);
            }
        }

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
            if (!QMath::IsZero(Qubit_Set_Value[BitAdd(count, Qubit_Situation, 0ULL)])) {
                int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Situation, 0ULL)]));
                Caculate *= Caculate;
                Caculate >>= Fixed_shift;
                Probability0 += static_cast<long long>(Caculate);
            }

            if (!QMath::IsZero(Qubit_Set_Value[BitAdd(count, Qubit_Situation, 1ULL)])) {
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
                if (!QMath::IsZero(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 0ULL)])) {
                    int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 0ULL)]));
                    Caculate *= Caculate;
                    Caculate >>= Fixed_shift;
                    Probability0 += static_cast<long long>(Caculate);
                }

                if (!QMath::IsZero(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 1ULL)])) {
                    int128_t Caculate = static_cast<int128_t>(QMath::AbsoluteValue(Qubit_Set_Value[BitAdd(count, Qubit_Set.size() - setting - 1, 1ULL)]));
                    Caculate *= Caculate;
                    Caculate >>= Fixed_shift;
                    Probability1 += static_cast<long long>(Caculate);
                }
            }

            if (Probability0 == 0 || Probability1 == 0)
                CutQubitSet(Qubit_Set[setting], Probability0, Probability1);
        }

        std::vector<int> remaining = Qubit_Set_Pointer->first;
        if (!remaining.empty() && Entangled_Qubit_Set[remaining[0]] != nullptr) {
            SplitEntangledGroup(remaining[0]);
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

        SplitEntangledGroup(q1);
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

        SplitEntangledGroup(q1);
    }

    // 將指定的 qubit 強制初始化為 |0> (target_state=0) 或 |1> (target_state=1)
    // 若該 qubit 處於糾纏狀態，會先從糾纏群組中投影分離，再設定目標狀態
    void InitializeQubit(int qubit, int target_state) {
        if (qubit < 0 || qubit >= Qubit_Amount) return;
        target_state = (target_state != 0) ? 1 : 0;

        // 若 qubit 目前在一個糾纏群組中，先將其投影分離
        if (Entangled_Qubit_Set.find(qubit) != Entangled_Qubit_Set.end() &&
            Entangled_Qubit_Set[qubit] != nullptr) {

            auto* group_ptr = Entangled_Qubit_Set[qubit];

            if (group_ptr->first.size() > 1) {
                int q_pos = GetSituation(qubit);
                std::vector<FixedComplex>& state = group_ptr->second;
                size_t half = state.size() >> 1;

                // 計算投影到 target_state 後的正規化係數
                int128_t norm_sq = 0;
                for (size_t c = 0; c < half; c++) {
                    FixedComplex amp = state[BitAdd(c, q_pos, (unsigned long long)target_state)];
                    int128_t av = QMath::AbsoluteValue(amp);
                    norm_sq += (av * av) >> Fixed_shift;
                }

                std::vector<FixedComplex> new_state(half, FixedComplex());

                if (norm_sq > 0) {
                    // target_state 有非零機率，正常投影
                    long long norm_factor = QMath::NewtonSqrt((int128_t)norm_sq << Fixed_shift);
                    for (size_t c = 0; c < half; c++) {
                        FixedComplex amp = state[BitAdd(c, q_pos, (unsigned long long)target_state)];
                        new_state[c] = amp / norm_factor;
                    }
                }
                else {
                    // target_state 在目前狀態下機率為零，改用另一個基底投影以保全剩餘 qubit 正規化
                    int other = 1 - target_state;
                    int128_t other_norm_sq = 0;
                    for (size_t c = 0; c < half; c++) {
                        FixedComplex amp = state[BitAdd(c, q_pos, (unsigned long long)other)];
                        int128_t av = QMath::AbsoluteValue(amp);
                        other_norm_sq += (av * av) >> Fixed_shift;
                    }
                    if (other_norm_sq > 0) {
                        long long norm_factor = QMath::NewtonSqrt((int128_t)other_norm_sq << Fixed_shift);
                        for (size_t c = 0; c < half; c++) {
                            FixedComplex amp = state[BitAdd(c, q_pos, (unsigned long long)other)];
                            new_state[c] = amp / norm_factor;
                        }
                    }
                }

                // 更新群組的狀態向量（移除本 qubit 的維度）
                group_ptr->second = new_state;

                // 從群組 qubit 名單中移除本 qubit
                auto& qubit_names = group_ptr->first;
                qubit_names.erase(
                    std::remove(qubit_names.begin(), qubit_names.end(), qubit),
                    qubit_names.end()
                );
            }

            // 切斷本 qubit 與原群組的連結
            Entangled_Qubit_Set[qubit] = nullptr;
        }

        // 清除已觀測標記
        Qubit_Set_Observation[qubit] = { 0, -1 };

        // 建立新的獨立狀態：|0> = {1,0}，|1> = {0,1}
        if (target_state == 0) {
            Entangled_Qubit_Set_Pointer[qubit] = { {qubit}, {FixedComplex(Fixed_Point), FixedComplex()} };
        }
        else {
            Entangled_Qubit_Set_Pointer[qubit] = { {qubit}, {FixedComplex(), FixedComplex(Fixed_Point)} };
        }
        Entangled_Qubit_Set[qubit] = &Entangled_Qubit_Set_Pointer[qubit];
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
        std::vector<double> dummy_errs(num_qubits, 0.0); // 提供不追蹤快照時的空白陣列

        while (ip < static_cast<int>(instructions.size())) {
            ErrorRecord dummy;
            int jump = ExecuteOne(instructions[ip], sim, lm, cfg, rng, dummy, dummy_errs);
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

        int cumulative_errors = 0;
        std::vector<double> qubit_error_accum(num_qubits, 0.0); // 針對每個 qubit 獨立追蹤誤差

        int ip = 0;
        while (ip < static_cast<int>(instructions.size())) {
            const auto& inst = instructions[ip];
            ErrorRecord err;
            int jump = ExecuteOne(inst, sim, lm, cfg, rng, err, qubit_error_accum);

            if (cfg.track_errors && err.occurred)
                cumulative_errors++;

            StepSnapshot snap = TakeSnapshot(inst, ip, sim);
            snap.error_record = err;
            snap.cumulative_errors = cumulative_errors;
            snap.qubit_error_accum = qubit_error_accum; // 儲存當下各 qubit 的狀態拷貝

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
        case GateType::Reset:       return "Reset";
        case GateType::InitQubit:   return "Init";
        default:                    return "?";
        }
    }

    static std::string InstructionToString(const CircuitInstruction& inst) {
        std::ostringstream ss;
        ss << "{" << GateTypeName(inst.type);
        if (!inst.label.empty())   ss << "," << inst.label;

        if (inst.type == GateType::Conditional) {
            if (inst.condition_is_parity) {
                ss << ",{Parity";
                for (int q : inst.qubits) ss << "," << q;
                if (inst.expected_val != 1) ss << "," << inst.expected_val;
                ss << "}";
            }
            else {
                ss << ",{" << (inst.condition_is_measure ? "M," : "")
                    << (inst.qubits.empty() ? -1 : inst.qubits[0]);
                if (inst.expected_val != 1) ss << "," << inst.expected_val;
                ss << "}";
            }
            for (const auto& s : inst.sub_insts) ss << " " << InstructionToString(s);
        }
        else {
            for (int q : inst.qubits) ss << ",Q" << q;
            if (inst.param != 0.0)
                ss << "," << std::fixed << std::setprecision(1) << inst.param;
            for (const auto& s : inst.sub_insts) ss << " " << InstructionToString(s);
        }
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
            {"Reset",GateType::Reset},{"RST",GateType::Reset},
            {"Init",GateType::InitQubit},{"INIT",GateType::InitQubit},
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
        case GateType::Reset:   return { 0,0 };
        case GateType::InitQubit: return { 1,1 };  // 1 qubit, 1 param (目標狀態 0 或 1)
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

        auto trim_str = [](const std::string& s) {
            size_t start = s.find_first_not_of(" \t\r\n");
            size_t end = s.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) return std::string("");
            return s.substr(start, end - start + 1);
            };

        if (inst.type == GateType::Label || inst.type == GateType::Goto) {
            inst.label = trim_str(rest); return inst;
        }

        if (inst.type == GateType::Conditional) {
            auto blocks = SplitTopLevelBlocks(rest);
            if (blocks.empty()) return inst;

            // 1. 解析條件區塊 blocks[0] (例如 "{M, 0}" 或 "{0, 0}")
            std::string cond_str = trim_str(blocks[0]);
            if (cond_str.front() == '{') cond_str.erase(0, 1);
            if (cond_str.back() == '}') cond_str.pop_back();

            auto cond_toks = SplitTokens(cond_str);
            if (!cond_toks.empty()) {
                std::string t0 = trim_str(cond_toks[0]);
                if (t0 == "M" || t0 == "Measure") {
                    inst.condition_is_measure = true;
                    if (cond_toks.size() > 1) inst.qubits.push_back(std::stoi(trim_str(cond_toks[1])));
                    if (cond_toks.size() > 2) inst.expected_val = std::stoi(trim_str(cond_toks[2]));
                    else inst.expected_val = 1; // 預設值
                }
                else if (t0 == "Parity") {
                    // 解析宇稱條件：所有數字為 qubit ID，
                    // 若末尾數字為 0 或 1 且移除後至少仍有 2 個 qubit，則視為 expected_val
                    inst.condition_is_parity = true;
                    inst.expected_val = 1; // 預設奇宇稱
                    std::vector<std::string> parity_toks(cond_toks.begin() + 1, cond_toks.end());
                    if (parity_toks.size() >= 3) { // 至少 3 個 token 才考慮末尾為 expected_val
                        int last = std::stoi(trim_str(parity_toks.back()));
                        if (last == 0 || last == 1) {
                            inst.expected_val = last;
                            parity_toks.pop_back();
                        }
                    }
                    for (const auto& t : parity_toks)
                        inst.qubits.push_back(std::stoi(trim_str(t)));
                    if (inst.qubits.size() < 2)
                        throw std::runtime_error("IF Parity 條件錯誤：至少需要 2 個 Qubit！");
                }
                else {
                    inst.condition_is_measure = false;
                    inst.qubits.push_back(std::stoi(t0));
                    if (cond_toks.size() > 1) inst.expected_val = std::stoi(trim_str(cond_toks[1]));
                    else inst.expected_val = 1; // 預設值
                }
            }

            // 2. 解析後續的子指令 blocks[1...N]
            for (size_t i = 1; i < blocks.size(); i++) {
                inst.sub_insts.push_back(ParseBlock(blocks[i]));
            }
            return inst;
        }

        std::pair<int, int> sig = GateSignature(inst.type);
        int nq = sig.first, np = sig.second;
        if (nq == -1) {
            for (const auto& t : SplitTokens(rest)) inst.qubits.push_back(std::stoi(trim_str(t)));
        }
        else {
            auto toks = SplitTokens(rest);
            for (int i = 0; i < nq && i < (int)toks.size(); i++) inst.qubits.push_back(std::stoi(trim_str(toks[i])));
            if (np > 0 && nq < (int)toks.size()) inst.param = std::stod(trim_str(toks[nq]));
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
        ErrorRecord& err, std::vector<double>& q_errs) {

        // ── 經典控制流（直接處理並回傳，不牽涉量子誤差） ──
        if (inst.type == GateType::Label) return -1;
        if (inst.type == GateType::Goto) {
            auto it = lm.find(inst.label);
            return (it != lm.end()) ? it->second : -1;
        }
        if (inst.type == GateType::Conditional) {
            if (inst.qubits.empty())
                throw std::runtime_error("IF 條件錯誤：未指定目標 Qubit！");

            // ── 宇稱條件：當場執行 MeasureParity，再判斷結果 ──
            if (inst.condition_is_parity) {
                int parity_result = sim.MeasureParity(inst.qubits);
                if (parity_result == -1)
                    throw std::runtime_error("IF Parity 條件錯誤：目標 Qubit 不可操作或已坍縮！");
                // 宇稱測量後，所有參與 qubit 的誤差歸零
                for (int pq : inst.qubits) q_errs[pq] = 0.0;
                if (parity_result == inst.expected_val) {
                    for (const auto& sub : inst.sub_insts) {
                        ErrorRecord se;
                        int jump = ExecuteOne(sub, sim, lm, cfg, rng, se, q_errs);
                        if (se.occurred && !err.occurred) err = se;
                        err.gate_error_prob += se.gate_error_prob;
                        if (jump >= 0) return jump;
                    }
                }
                return -1;
            }

            int q = inst.qubits[0];

            if (inst.condition_is_measure) {
                // 如果條件區塊要求當下測量，則執行測量並歸零誤差
                sim.ObserverQubit(q);
                q_errs[q] = 0.0;
            }
            else {
                // 如果只提供 qubit 標籤，則必須已經坍縮
                if (!sim.IsObservered(q)) {
                    throw std::runtime_error("IF 條件錯誤：目標 Qubit 尚未坍縮，且未在 IF 條件中指定測量！");
                }
            }

            int meas = sim.GetMeasurementResult(q);
            if (meas == inst.expected_val) {
                for (const auto& sub : inst.sub_insts) {
                    ErrorRecord se;
                    int jump = ExecuteOne(sub, sim, lm, cfg, rng, se, q_errs);
                    if (se.occurred && !err.occurred) err = se;
                    err.gate_error_prob += se.gate_error_prob;
                    if (jump >= 0) return jump;
                }
            }
            return -1;
        }

        // ── 量子操作與誤差前置計算 ──
        bool is2 = IsTwoQubitGate(inst.type);
        bool is3 = IsThreeQubitGate(inst.type);

        double p = cfg.error_rate;
        if (is3) p = std::min(1.0, p * cfg.three_qubit_error_mult);

        double e_comb = 0.0;
        if (is2 || is3) {
            err.gate_error_prob = p;
            if (is2 && inst.qubits.size() >= 2) {
                int q1 = inst.qubits[0], q2 = inst.qubits[1];
                e_comb = 1.0 - (1.0 - q_errs[q1]) * (1.0 - q_errs[q2]) * (1.0 - p);
            }
            else if (is3 && inst.qubits.size() >= 3) {
                int q1 = inst.qubits[0], q2 = inst.qubits[1], q3 = inst.qubits[2];
                e_comb = 1.0 - (1.0 - q_errs[q1]) * (1.0 - q_errs[q2]) * (1.0 - q_errs[q3]) * (1.0 - p);
            }
        }

        // ── 執行量子邏輯閘 ──
        switch (inst.type) {
        case GateType::H:  case GateType::X:  case GateType::Y:
        case GateType::Z:  case GateType::S:  case GateType::Sdg:
        case GateType::T:  case GateType::Tdg:
            sim.ApplySingleQubitGate(inst.qubits[0], inst.type); break;
        case GateType::Rx: case GateType::Ry: case GateType::Rz:
            sim.ApplySingleQubitGate(inst.qubits[0], inst.type, inst.param); break;
        case GateType::CNOT:  case GateType::CZ:
        case GateType::SWAP:  case GateType::iSWAP: case GateType::SqrtSWAP:
            sim.ApplyTwoQubitGate(inst.qubits[0], inst.qubits[1], inst.type);
            HandleError(inst, sim, cfg, rng, err); break;
        case GateType::CCNOT: case GateType::CSWAP:
            sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1], inst.qubits[2], inst.type);
            HandleError(inst, sim, cfg, rng, err); break;
        case GateType::Deutsch:
            sim.ApplyThreeQubitGate(inst.qubits[0], inst.qubits[1], inst.qubits[2], inst.type, inst.param);
            HandleError(inst, sim, cfg, rng, err); break;
        case GateType::Measure:
            sim.ObserverQubit(inst.qubits[0]); break;
        case GateType::Parity:
            sim.MeasureParity(inst.qubits);    break;
        case GateType::Reset:
            // 將所有 Qubit 重置回 |0> 並清除所有誤差累積
            sim.ReinitializeAllQubits();
            std::fill(q_errs.begin(), q_errs.end(), 0.0);
            return -1;
        case GateType::InitQubit:
            // 將指定 Qubit 初始化為 |0> 或 |1>，並清除該 Qubit 的誤差
            if (!inst.qubits.empty()) {
                sim.InitializeQubit(inst.qubits[0], static_cast<int>(inst.param));
                q_errs[inst.qubits[0]] = 0.0;
            }
            break;
        default: break;
        }

        // ── 誤差擴散與脫離糾纏歸零機制 ──
        auto groups = sim.GetEntangledGroups();
        for (const auto& g : groups) {
            // 規則一：如果該 Qubit 脫離糾纏成為獨立狀態 (Size == 1)，其誤差立刻歸零
            if (g.qubit_ids.size() == 1) {
                q_errs[g.qubit_ids[0]] = 0.0;
            }
            // 規則二：若剛才執行多量子閘，且該群組包含剛才被操作的目標位元，將新誤差「擴散」給該群組內「所有」位元
            else if (is2 || is3) {
                bool hit = false;
                for (int t : inst.qubits) {
                    if (std::find(g.qubit_ids.begin(), g.qubit_ids.end(), t) != g.qubit_ids.end()) {
                        hit = true; break;
                    }
                }
                if (hit) {
                    for (int q : g.qubit_ids) {
                        q_errs[q] = e_comb;
                    }
                }
            }
        }

        // 規則三：測量歸零 (包含單點測量與宇稱測量)
        if (inst.type == GateType::Parity) {
            for (int q : inst.qubits) {
                q_errs[q] = 0.0; // 參與宇稱測量的 Qubit 誤差強制歸零
            }
        }

        for (int q = 0; q < sim.GetQubitAmount(); ++q) {
            if (sim.IsObservered(q)) q_errs[q] = 0.0;
        }

        return -1;
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
        out << Row('=', W) << "\n";
    }

    // ── 電路圖 ────────────────────────────────────────────────

    static void PrintCircuitDiagram(const CircuitExpression& circuit,
        std::ostream& out) {
        int nq = circuit.num_qubits;
        const int CW = 9;

        // 更新：加入 is_exclusive 標記
        struct Slot {
            std::vector<bool>               occ;
            std::vector<CircuitInstruction> cell;
            bool                            is_exclusive;
            Slot(int n) : occ(n, false), cell(n), is_exclusive(false) {}
        };
        std::vector<Slot> slots;

        for (const auto& inst : circuit.instructions) {
            std::vector<int> tgts = CollectQubits(inst, nq);
            bool standalone = (inst.type == GateType::Label ||
                inst.type == GateType::Goto ||
                inst.type == GateType::Conditional);
            int si = -1;

            if (!standalone && !tgts.empty()) {
                bool is_multi_qubit = (tgts.size() >= 2);
                int max_s = -1;

                if (is_multi_qubit) {
                    // 若是多量子邏輯閘：強制排在現有所有 Slot 的下一個，不與任何人並列
                    max_s = (int)slots.size() - 1;
                }
                else {
                    // 若是單量子邏輯閘：找出目標 qubit 被佔用的最後一個 Slot，
                    // 同時遇到被多量子邏輯閘「獨佔 (is_exclusive)」的 Slot 時視為一堵牆，不能往前插隊。
                    for (int q : tgts) {
                        for (int s = (int)slots.size() - 1; s >= 0; s--) {
                            if (slots[s].occ[q] || slots[s].is_exclusive) {
                                max_s = std::max(max_s, s);
                                break;
                            }
                        }
                    }
                }

                // 決定要放入的 Slot 索引
                si = max_s + 1;

                // 若超出現有 Slot 數量，新增 Slot
                while (si >= (int)slots.size()) {
                    slots.push_back(Slot(nq));
                }

                // 如果是多量子位元閘，就把這個時間點設為獨佔
                if (is_multi_qubit) {
                    slots[si].is_exclusive = true;
                }

                // 寫入佔用狀態與指令
                for (int q : tgts) {
                    slots[si].occ[q] = true;
                    slots[si].cell[q] = inst;
                }
            }
        }

        // 底下印出 T0, T1 以及 Q[x] 的邏輯維持原樣不變
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
        case GateType::Reset:  return "[RST]";
        case GateType::InitQubit: {
            int st = static_cast<int>(inst.param);
            return (st == 0) ? "[I|0>]" : "[I|1>]";
        }
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
        if (result.config.track_errors) {
            bool has_error_risk = (snap.error_record.gate_error_prob > 0.0);
            bool has_accum_err = false;
            for (double e : snap.qubit_error_accum) { if (e > 1e-9) has_accum_err = true; }

            if (has_error_risk || has_accum_err) {
                out << " [Error]";
                if (has_error_risk) out << " gate_prob=" << FmtPct(snap.error_record.gate_error_prob);
                out << "  total_hits=" << snap.cumulative_errors << "\n";

                out << "  -> Qubit Errors: ";
                for (int q = 0; q < num_qubits; q++) {
                    if (snap.qubit_error_accum[q] > 1e-9) {
                        out << "Q" << q << "=" << FmtPct(snap.qubit_error_accum[q]) << "  ";
                    }
                }
                out << "\n";

                if (snap.error_record.occurred)
                    out << "  *** " << ErrorTypeName(snap.error_record.type)
                    << " applied on Q" << snap.error_record.affected_qubit << " ***\n";
                out << "\n";
            }
        }

        // 測量結果
        bool any_meas = false;
        // ... (下方這段印出測量結果與 Quantum State 的迴圈完全保持原樣不變)
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
        const auto& final_errs = result.snapshots.back().qubit_error_accum;

        out << " Total physical error occurrences: " << total_errors << "\n";
        out << " Final theoretical error probability per Qubit:\n";
        for (int q = 0; q < result.num_qubits; q++) {
            if (final_errs[q] > 1e-9) {
                out << "   Q[" << q << "]: " << FmtPct(final_errs[q]) << "\n";
            }
            else {
                out << "   Q[" << q << "]: 0.00% (Clean or Collapsed)\n";
            }
        }
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
        std::cout << "  config         : 查看目前設定\n";
        std::cout << "  run \"<expr>\"   : 執行電路，表達式必須用雙引號包裝\n";
        std::cout << "                   (例如: run \"3, {H,0}, {CNOT,0,1}\")\n";
        std::cout << "  show           : 顯示最後一次執行的電路圖與狀態快照\n";
        std::cout << "  exit / quit    : 離開控制台\n\n";
    }

    void CmdManual() {
        std::cout << "\n================ 量子電路操作手冊 ================\n";
        std::cout << "[基本架構與語法]\n";
        std::cout << "本模擬器用於模擬量子位元 (Qubits) 的演化。\n";
        std::cout << "指令語法必須用雙引號 \"\" 包起來，格式如下：\n";
        std::cout << "  \"<總量子位元數>, {量子閘, 目標1, 目標2...}, {量子閘...}\"\n\n";

        std::cout << "[所有可用的量子電路組件 (GateTypes)]\n";
        std::cout << "1. 基本單量子位元閘 (需 1 個目標位元):\n";
        std::cout << "   - H    : Hadamard 閘 (創造均勻疊加態)\n";
        std::cout << "   - X,Y,Z: Pauli 閘系列 (X為量子NOT，Z為相位反轉)\n";
        std::cout << "   - S,Sdg: 相位旋轉 90 度與其共軛\n";
        std::cout << "   - T,Tdg: 相位旋轉 45 度與其共軛\n\n";

        std::cout << "2. 參數化單量子位元閘 (需 1 個角度參數，單位: 度數):\n";
        std::cout << "   - Rx, Ry, Rz : 繞對應軸旋轉 (例如: {Rx,0,90})\n\n";

        std::cout << "3. 雙/三量子位元閘:\n";
        std::cout << "   - CNOT, CZ       : 雙位元控制閘\n";
        std::cout << "   - SWAP系列       : SWAP, iSWAP, SqrtSWAP\n";
        std::cout << "   - CCNOT, CSWAP   : 三位元控制閘\n";
        std::cout << "   - Deutsch        : 需額外 1 個角度 (例如: {Deutsch,0,1,2,45})\n\n";

        std::cout << "4. 測量與經典控制流 (Measurement & Flow):\n";
        std::cout << "   - M / Measure: 測量單一量子位元並坍縮。測量後該位元退出糾纏。\n";
        std::cout << "   - Parity     : 測量多個量子位元的宇稱。\n";
        std::cout << "   - Label, Goto: 跳轉標籤與執行跳轉。\n";
        std::cout << "   - Reset / RST: 將「全部」量子位元重置回 |0> 初始狀態，清除所有糾纏與測量紀錄。\n";
        std::cout << "                  格式: {Reset}  (無需任何參數)\n";
        std::cout << "                  常見用途: 在 IF 分支內確保 Qubit 回到已知狀態後再繼續操作。\n";
        std::cout << "                  範例: {If, {M, 0}, {Reset}, {H, 1}}  (若 Q0=1 則重置全部再做 H 閘)\n";
        std::cout << "   - Init / INIT : 將「指定」量子位元初始化為 |0> 或 |1>，清除其糾纏與測量紀錄。\n";
        std::cout << "                  格式: {Init, q, state}  (state = 0 表示 |0>，state = 1 表示 |1>)\n";
        std::cout << "                  若該 qubit 處於糾纏狀態，會先對群組做投影分離，其餘 qubit 狀態自動正規化。\n";
        std::cout << "                  範例1: {Init, 0, 0}      (將 Q0 初始化為 |0>)\n";
        std::cout << "                  範例2: {Init, 1, 1}      (將 Q1 初始化為 |1>)\n";
        std::cout << "                  範例3: {H,0},{CNOT,0,1},{Init,0,0}   (Bell 態後強制將 Q0 重設為 |0>，Q1 會坍縮)\n";
        std::cout << "   - If         : 條件分支 (條件可為動態測量，或檢查已坍縮的狀態)\n";
        std::cout << "                  格式: {If, {條件}, {指令1}, {指令2}, ...}\n";
        std::cout << "                  子指令可以是任意普通指令（H/X/CNOT/Reset 等），也可以再嵌套 If。\n";
        std::cout << "                  [條件寫法]\n";
        std::cout << "                  1. {M, q, val}          : 對 q 測量，預期值為 val (若省略 val 則預設為 1)\n";
        std::cout << "                  2. {q, val}             : 檢查已坍縮的 q，預期值為 val (若省略 val 則預設為 1)\n";
        std::cout << "                  3. {Parity, q0, q1, ...}: 當場測量多個 qubit 的宇稱，預期值預設為 1 (奇宇稱)\n";
        std::cout << "                     若末尾多一個 0 且 qubit 數 >= 2，則視為 expected_val=0 (偶宇稱)\n";
        std::cout << "                  範例1: {If, {M, 0}, {X, 1}}                       (測量 Q0，若=1 則對 Q1 做 X)\n";
        std::cout << "                  範例2: {If, {M, 0, 0}, {X, 1}}                    (測量 Q0，若=0 則對 Q1 做 X)\n";
        std::cout << "                  範例3: {M, 0}, {If, {0}, {X, 1}}                  (先測量 Q0，後判斷)\n";
        std::cout << "                  範例4: {If, {M,0}, {If, {M,1}, {X,2}}}            (嵌套 If)\n";
        std::cout << "                  範例5: {If, {M,0}, {Reset}, {H,1}}                (Q0=1 則重置後做 H)\n";
        std::cout << "                  範例6: {If, {Parity,0,1}, {X,2}}                  (Q0,Q1 宇稱為奇(1)則對 Q2 做 X)\n";
        std::cout << "                  範例7: {If, {Parity,0,1,2,0}, {Z,3}}              (Q0,Q1,Q2 宇稱為偶(0)則對 Q3 做 Z)\n\n";

        std::cout << "[經典示範 (Examples)]\n";
        std::cout << " 量子遙傳 (Quantum Teleportation)\n";
        std::cout << "  輸入: run \"3, {H,1}, {CNOT,1,2}, {CNOT,0,1}, {H,0}, {M,0}, {M,1}, {If,1,1,{X,2}}, {If,0,1,{Z,2}}\"\n";
        std::cout << "==================================================\n\n";
    }

    void CmdConfig(const std::string& args) {
        if (args.empty()) {
            std::cout << "\n[當前設定]\n";
            std::cout << "  (目前無可調整的設定項目)\n\n";
            return;
        }

        std::stringstream ss(args);
        std::string key;
        ss >> key;

        if (key == "error_rate" || key == "apply_errors" || key == "track_errors") {
            std::cout << "提示：誤差相關功能目前已停用，此設定不會產生任何效果。\n";
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
        // 初始化預設 config（誤差追蹤功能已停用）
        config.track_errors = false;
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