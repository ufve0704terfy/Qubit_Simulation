# Qubit Simulation Console (量子電路模擬器)

這是一個以 C++ 開發的高效能量子電路模擬器。本專案透過一個互動式的命令列控制台（QConsole），讓使用者可以動態輸入、解析並執行量子電路[cite: 2]。模擬器底層支援量子疊加、量子糾纏的狀態追蹤，並包含了自定義的複數與定點數（Fixed-Point）運算以確保機率幅計算的精度[cite: 2]。

## ✨ 核心特色 (Features)

* **互動式控制台 (Interactive CLI)**：提供友善的操作介面，即時解析並執行量子電路字串[cite: 2]。
* **豐富的量子閘支援**：內建包含 Pauli 閘、參數化旋轉閘、CNOT/SWAP 等雙位元閘，以及 Toffoli 等三位元閘[cite: 2]。
* **動態電路與控制流**：支援測量（Measurement）、宇稱檢查（Parity），以及基於動態測量或經典狀態讀取的條件執行（If）和標籤跳轉（Goto）[cite: 2]。
* **真實物理誤差模型 (Physical Error Diffusion)**：內建可配置的量子雜訊模擬[cite: 2]。支援誤差在糾纏群組中的局部擴散；當量子位元被測量或脫離糾纏成為獨立狀態時，其物理誤差將真實地歸零[cite: 2]。

## 🚀 快速開始 (Quick Start)

### 下載執行檔 (免編譯)
1. 前往本專案的 Releases 頁面[cite: 2]。
2. 下載最新的 `Qubit_Simulation.exe`[cite: 2]。
3. 雙擊執行即可進入 QConsole[cite: 2]。

### 從原始碼編譯 (Build from Source)
本專案使用標準 C++ 開發，推薦使用 **Visual Studio 2022** 進行編譯[cite: 2]：
1. 複製此儲存庫：`git clone <your-repo-url>`[cite: 2]
2. 使用 Visual Studio 2022 開啟專案資料夾[cite: 2]。
3. 將建置組態切換為 **Release** 與 **x64**[cite: 2]。
4. 為了確保執行檔的獨立性，建議在專案屬性中將 C/C++ -> 程式碼產生 -> 執行階段程式庫設為 `多執行緒 (/MT)`[cite: 2]。
5. 點擊「建置方案」，生成的 `.exe` 將位於 `x64/Release/` 目錄下[cite: 2]。

## 📖 控制台指令說明

進入 QConsole 後，你可以使用以下基本指令[cite: 2]：

* `help`：顯示基本指令列表[cite: 2]。
* `manual`：顯示完整的量子電路操作手冊與支援的元件清單[cite: 2]。
* `config`：查看或修改模擬參數（例如：`config error_rate 0.05` 或 `config apply_errors true`）[cite: 2]。
* `run "<expr>"`：執行量子電路（**注意：表達式必須包含在雙引號中**）[cite: 2]。
* `show`：顯示最後一次執行的電路圖與狀態快照[cite: 2]。
* `exit` / `quit`：離開程式[cite: 2]。

## 🛠️ 電路語法與量子閘列表

電路表達式的基本語法為：`"<總量子位元數>, {量子閘, 目標1, ...}, {量子閘...}"`[cite: 2]

### 支援的量子組件：
1. **單量子位元閘**[cite: 2]
   * `{H,0}`: Hadamard 閘[cite: 2]
   * `{X,0}`, `{Y,0}`, `{Z,0}`: Pauli 閘[cite: 2]
   * `{S,0}`, `{Sdg,0}`, `{T,0}`, `{Tdg,0}`: 相位閘[cite: 2]
   * `{Rx,0,90}`, `{Ry,0,90}`, `{Rz,0,90}`: 參數化旋轉閘（第三個參數為度數）[cite: 2]
2. **雙量子位元閘**[cite: 2]
   * `{CNOT,0,1}`: 控制反轉閘（控制位元, 目標位元）[cite: 2]
   * `{CZ,0,1}`: 控制相位閘[cite: 2]
   * `{SWAP,0,1}`, `{iSWAP,0,1}`, `{SqrtSWAP,0,1}`: 交換系列閘[cite: 2]
3. **三量子位元閘**[cite: 2]
   * `{CCNOT,0,1,2}`: Toffoli 閘（控1, 控2, 目標）[cite: 2]
   * `{CSWAP,0,1,2}`: Fredkin 閘[cite: 2]
   * `{Deutsch,0,1,2,45}`: Deutsch 閘（含角度參數）[cite: 2]
4. **測量與控制流**[cite: 2]
   * `{M,0}` 或 `{Measure,0}`: 測量特定位元（測量後退出糾纏，累積誤差歸零）[cite: 2]
   * `{Parity,0,1,2}`: 測量多位元的宇稱（參與測量的位元誤差歸零）[cite: 2]
   * `{Label,L1}`, `{Goto,L1}`: 標籤與無條件跳轉[cite: 2]
   * `{If, {條件}, {指令}}`: 條件執行分支。條件區塊支援以下兩種寫法[cite: 2]：
     * **當下測量**：`{M, 0, 1}` 對 Q0 進行測量，預期值為 1 (省略預期值則預設為 1)。
     * **讀取已坍縮狀態**：`{0, 1}` 檢查已坍縮的 Q0，預期值為 1 (省略預期值則預設為 1)。
     * *範例*：`{If, {M, 0}, {X, 1}}`（當下測量 Q0，若為 1 則對 Q1 執行 X 閘）[cite: 2]。

## 💡 經典示範 (Examples)

**1. 建立貝爾態 (Bell State - 最大糾纏態)**[cite: 2]
準備兩個量子位元，將其變為 $|00\rangle$ 與 $|11\rangle$ 各佔 50% 機率的完美糾纏狀態[cite: 2]：
```text
QConsole> run "2, {H,0}, {CNOT,0,1}"
QConsole> show