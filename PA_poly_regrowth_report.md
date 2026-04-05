# Memoryless Polynomial PA — 頻譜再生長分析報告

**日期**: 2026-04-05
**Script**: `PA_poly_regrowth.m`
**Figures**: `PA_poly_input.png`, `PA_poly_odd.png`, `PA_poly_all.png`, `PA_poly_grid.png`

---

## 1. 研究目的

用最簡單的 PA 模型 — **memoryless polynomial** — 搞清楚三件事：

1. 每個 polynomial order $k$ 造成的 spectral regrowth 落在哪個 frequency band？
2. Odd order 跟 even order 的行為有什麼差別？
3. 常聽到「even order 對 ACLR 沒貢獻」— 這句話到底對不對？什麼意義下才對？

---

## 2. Signal 與 PA 模型設定

### 2.1 Input Signal

- 取樣率 $F_s = 240$ MHz
- 訊號頻寬 $\text{BW} = 20$ MHz (complex baseband)
- 產生方式：在頻域放一個 $[-10, +10]$ MHz 的方形 mask，每個 bin 配**隨機相位**，然後做 IFFT。
  得到的 time-domain 訊號是 **bandlimited 隨機訊號**，PSD 在 in-band 平坦，out-of-band 是數值零。
- 正規化到 peak $= 1$，PAPR $\approx 9.35$ dB (類似 Rayleigh 分布)。

### 2.2 Memoryless Polynomial PA Model

```math
y(n) \;=\; \sum_{k=1}^{K} a_k \cdot x(n) \cdot |x(n)|^{k-1}
```

其中 $x(n)$ 是 complex baseband input，$|x(n)| = \sqrt{x(n)\,x^{\ast}(n)}$ 是 envelope magnitude。

**這次分析取 $a_k = 1$**，專注在每一項 basis function 的「形狀」與「位置」，不 fit 特定 PA。

第 $k$ 項的 basis function：

```math
y_k(n) \;\triangleq\; x(n) \cdot |x(n)|^{k-1}
```

---

## 3. 每個 Order 的 Bandwidth — Baseband 觀點

### 3.1 關鍵觀察：Odd 和 Even 是「不同類型」的 function

**Odd $k$** ($k = 1, 3, 5, \ldots$)：

```math
y_k = x \cdot |x|^{k-1} = x \cdot (x\,x^{\ast})^{\frac{k-1}{2}} = x^{\frac{k+1}{2}} \cdot (x^{\ast})^{\frac{k-1}{2}}
```

這是 $x$ 和 $x^{\ast}$ 的**純多項式** (polynomial)，joint degree $= k$。

**Even $k$** ($k = 2, 4, 6, \ldots$)：

```math
y_k = x \cdot |x|^{k-1} = x \cdot |x| \cdot (x\,x^{\ast})^{\frac{k-2}{2}}
```

含 $|x| = \sqrt{x\,x^{\ast}}$，這是一個**非多項式** (non-polynomial, has a square root)。

### 3.2 Bandwidth 推導

若 $x(t)$ 的 spectrum 支撐在 $[-\tfrac{B}{2}, \tfrac{B}{2}]$ ($B = \text{BW}$)，則：

- 時域乘法 $\longleftrightarrow$ 頻域捲積
- $x^a(t) \cdot (x^{\ast})^b(t)$ 的 spectrum 是 $a$ 個 $X(f)$ 跟 $b$ 個 $X^{\ast}(-f)$ 的 $(a+b)$-fold 捲積
- 捲積支撐相加：**bandwidth $= (a+b) \cdot B$**

**Odd $k$**： $a = \tfrac{k+1}{2}$, $b = \tfrac{k-1}{2}$, $a+b = k$

```math
\boxed{\text{BW}(y_k) \;=\; k \cdot \text{BW}, \qquad \text{支撐嚴格在 } \left[-\tfrac{k \cdot \text{BW}}{2},\, +\tfrac{k \cdot \text{BW}}{2}\right]}
```

**Even $k$**：含 $|x|$ 這個 sqrt，無法寫成有限項多項式

```math
\boxed{\text{BW}(y_k) \;=\; \infty \text{ (long tail)}, \qquad \text{但能量集中在 } \sim k \cdot \text{BW} \text{ 附近}}
```

> **這是本 report 最重要的一句話**：在 baseband polynomial 模型裡，odd 跟 even order 的 basis function **不是同一種數學物件**。Odd 是 polynomial，有精確 bandwidth；Even 含 sqrt，有 infinite tail。

---

## 4. 數值結果

### 4.1 Per-Order Power Table (from MATLAB)

| Order $k$ | $P_\text{in}$ (dB) | ACLR1 (dBc) | ACLR2 (dBc) | 預期 BW |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 0.00 | −39.06 | **−∞** | 20 MHz |
| **2** | −6.64 | −13.79 | **−37.48** | 40 MHz |
| 3 | −11.78 | −8.95 | **−∞** | 60 MHz |
| **4** | −15.85 | −6.42 | **−32.62** | 80 MHz |
| 5 | −19.13 | −4.78 | −25.46 | 100 MHz |
| **6** | −21.83 | −3.62 | −21.06 | 120 MHz |
| 7 | −24.11 | −2.76 | −17.90 | 140 MHz |
| **8** | −26.07 | −2.09 | −15.47 | 160 MHz |
| 9 | −27.80 | −1.55 | −13.53 | 180 MHz |

(粗體 row = even order; In-band 定義 $|f| \le 10$ MHz; ACLR1 是 $10 < |f| \le 30$ MHz; ACLR2 是 $30 < |f| \le 50$ MHz)

### 4.2 最關鍵的觀察

看 **ACLR2** 那一欄：

- **$k = 1$**： $-\infty$ ← 訊號只有 $\pm 10$ MHz，ACLR2 band (30~50 MHz) 完全在外
- **$k = 3$**： $-\infty$ ← $y_3 = x \cdot |x|^2 = x^2 \cdot x^{\ast}$ 是**純多項式**，bandwidth 嚴格 $= 60$ MHz，超出 $\pm 30$ MHz 之後 **exact zero** → ACLR2 band 真的什麼都沒有
- **$k = 2$**： $-37$ dBc (**非零！**) ← $y_2 = x \cdot |x|$ 含 sqrt，有 long tail，在 30~50 MHz 還有 $-37$ dB 的能量洩漏
- **$k = 4$**： $-32$ dBc (非零)
- **$k = 6, 8$**：都有明顯的 ACLR2 能量

這完全印證 §3 的數學推導：**odd order basis 有 hard bandwidth limit，even order basis 有 infinite tail**。

### 4.3 Figures 對應

- `PA_poly_odd.png` (Figure 1)：五條 odd order curves，normalize 到各自 peak。可以看到 $k=1, 3, 5, 7, 9$ 的 bandwidth 是 20, 60, 100, 140, 180 MHz，**每條都有清楚的 hard edge**。
- `PA_poly_all.png` (Figure 2)：九條 curves，odd=實線、even=虛線。仔細看就會發現 even (虛線) 的 skirt 拖得比相鄰 odd 的更遠。
- `PA_poly_grid.png` (Figure 3)：3×3 subplot，每格一個 order，ref 是 $k=1$ peak。Odd 的 subplot 是「梯形 + hard cutoff」，Even 的 subplot 是「bell curve + long tail」，差別一眼可見。

---

## 5. RF Passband 觀點 — 為什麼「Even Order 對 ACLR 沒貢獻」的 intuition 是對的

### 5.1 設定

實體 PA 吃的是 **real passband signal**：

```math
x_\text{RF}(t) \;=\; \mathrm{Re}\!\left\{ x(t)\,e^{j\omega_c t} \right\} \;=\; \tfrac{1}{2}\!\left( x(t)\, e^{j\omega_c t} + x^{\ast}(t)\, e^{-j\omega_c t} \right)
```

PA 是 real polynomial：

```math
y_\text{RF}(t) \;=\; \sum_{k=1}^{K} a_k \cdot x_\text{RF}(t)^k
```

### 5.2 Binomial 展開

對每個 $k$ 用 binomial theorem：

```math
x_\text{RF}(t)^k \;=\; \frac{1}{2^k} \sum_{m=0}^{k} \binom{k}{m}\, x^m (x^{\ast})^{k-m}\, e^{\,j(2m-k)\,\omega_c t}
```

每一項的 carrier frequency 是 $(2m - k)\,\omega_c$，$m = 0, 1, \ldots, k$。

所以 $x_\text{RF}^k$ 出現的載波頻率集合是：

```math
\{\,-k,\; -k+2,\; -k+4,\; \ldots,\; k-2,\; k\,\} \cdot \omega_c
```

### 5.3 Odd vs Even 的關鍵差別

| $k$ | 頻率集合 ($\times \omega_c$) | 有沒有 $\pm 1$？ | 經 BPF 後 in-band 還有東西嗎？ |
|:---:|:---|:---:|:---:|
| 1 | $\{-1, +1\}$ | ✓ | **有** |
| 2 | $\{-2, 0, +2\}$ | ✗ | **沒** |
| 3 | $\{-3, -1, +1, +3\}$ | ✓ | **有** |
| 4 | $\{-4, -2, 0, +2, +4\}$ | ✗ | **沒** |
| 5 | $\{-5, -3, -1, +1, +3, +5\}$ | ✓ | **有** |
| 6 | $\{-6, -4, -2, 0, +2, +4, +6\}$ | ✗ | **沒** |
| $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ |

**一般規則**：
- **Odd $k$** → $2m-k$ 是奇數 → 集合**包含** $\pm 1$ → 有 $\pm \omega_c$ 的 component
- **Even $k$** → $2m-k$ 是偶數 → 集合**只有** $0, \pm 2, \pm 4, \ldots$ → **沒有** $\pm \omega_c$

### 5.4 具體算 $k = 2$ 和 $k = 3$

**$k = 2$**（Even）：

```math
x_\text{RF}^2 \;=\; \tfrac{1}{4}\Big( x^2 e^{j 2\omega_c t} \;+\; 2|x|^2 \;+\; (x^{\ast})^2 e^{-j 2\omega_c t} \Big)
```

成分在 $\{0, \pm 2\omega_c\}$，**$\omega_c$ 附近完全沒東西**。
DC 項 $\tfrac{1}{2}|x|^2$ 是 envelope squared，是個 **baseband 低頻訊號**（頻寬 $0$ 到 $2\cdot\text{BW}$），**不在 fundamental band**。

**$k = 3$**（Odd）：

```math
x_\text{RF}^3 \;=\; \tfrac{1}{8}\Big( x^3 e^{j 3\omega_c t} \;+\; 3|x|^2 x\, e^{j\omega_c t} \;+\; 3|x|^2 x^{\ast}\, e^{-j\omega_c t} \;+\; (x^{\ast})^3 e^{-j 3\omega_c t} \Big)
```

在 $\omega_c$ 這邊的 term 是 $3|x|^2 x$。**這就是為什麼 baseband model 裡 $k=3$ 的 basis 是 $x \cdot |x|^2$** — 它來自 $x_\text{RF}^3$ 的 $e^{j\omega_c t}$ component。

### 5.5 RF 圖像

```
頻譜位置：        0       fc      2fc     3fc     4fc     5fc     ...
                 |       |       |       |       |       |
k=1  (odd)               ●                                         ← fundamental
k=2  (even)      ●               ●                                 ← DC + 2fc
k=3  (odd)               ●               ●                         ← fund + 3rd harm
k=4  (even)      ●               ●               ●                 ← DC + 2fc + 4fc
k=5  (odd)               ●               ●               ●         ← fund + 3rd + 5th
                         ↑
                         └── 只有這裡會被送到天線 (BPF passband)
```

**Even order 不產生 $\omega_c$ 附近的 component** → 通過 bandpass filter 後，只有 odd order 在 fundamental band 留下 distortion。這就是「even order 不貢獻 ACLR」的 physical intuition。

---

## 6. 是誰「濾掉」Even Order？

不是 PA 晶體管本身 — 晶體管非線性會**忠實產生**所有這些 even-order products (在 DC, $2f_c$, $4f_c$, …)。

**真正擋掉 even order 在輸出 spectrum 上出現的，是 PA 後面這條 RF 路徑上的 bandpass filtering**：

1. **Output matching network** — narrowband 設計，在 $f_c$ 附近共振，對 DC 跟 $2f_c$ 阻抗嚴重失配
2. **Harmonic termination / balun**
3. **TX band filter** (SAW / BAW / cavity)
4. **天線本身的 bandwidth**

所以精確的說法是：

> **「PA 輸出端到天線這整條路徑的 bandpass filtering，把 even-order 在 DC / $2f_c$ / $4f_c$ 的 products 擋掉了。在天線發射出去的 spectrum 裡，ACLR 只有 odd-order 的貢獻。」**

---

## 7. 重要 Caveat — Even Order 的間接影響

Even order **不是完全沒影響**，只是它的影響**不是**走「直接在 fundamental band 產生 distortion」這條路。真正的影響走三條間接路徑：

### 7.1 Envelope / Bias Memory Effect (最重要)

$x_\text{RF}^2$ 的 DC 項是 $\tfrac{1}{2}|x(t)|^2$，這是 envelope squared，頻譜在 $0 \sim 2\cdot\text{BW}$ 之間（baseband）。

如果 PA 的 bias feed 在這個頻率範圍 **video bandwidth 不夠**（bypass cap 不夠、bias inductor 有寄生），envelope 電壓就會**調變到 bias 上**，然後 re-modulate fundamental band，結果：

```math
V_\text{bias}(t) = V_\text{bias}^{(0)} + \alpha \cdot |x(t)|^2 \;\Longrightarrow\; \text{Gain}(t) = G_0 + \beta \cdot |x(t)|^2
```

→ 產生**類似 AM-AM 的動態失真**，而且**有 memory**（因為 bias network 有 RC time constant）。

**這是 wideband PA 最頭痛的一種 memory effect，根源就是 even order。**

### 7.2 Thermal Memory

$|x|^2$ 是瞬時功率。PA 的晶體管 junction 溫度跟著 $|x|^2$ 的低頻分量起伏（typical thermal time constant 是 $\mu$s 到 ms 級），溫度變化又改變 transconductance → 產生 thermal-induced distortion。

→ 這也來自 even order 的 DC/envelope 成分，但頻率更低，需要更長 memory。

### 7.3 Harmonic Reflection

Output matching 不完美 → $2f_c$ 的 reflection 回到 PA input → 再被放大一次 → 跟 fundamental 混合產生新的 in-band product。這在 high-power PA 或 wideband PA 特別明顯。

→ 還是 even order 在作怪，只是走 reflection 這條路徑。

---

## 8. 對 DPD 設計的 Implication

### 8.1 為什麼 DPD polynomial basis 只用 odd order

| DPD 方法 | Basis | 處理 even order 的方式 |
|:---|:---|:---|
| **Memoryless polynomial** | $x,\; x \cdot \lvert x \rvert^2,\; x \cdot \lvert x \rvert^4,\; \ldots$ | 不加 even basis — direct RF contribution = 0 |
| **Memory polynomial** | $x(n-q) \cdot \lvert x(n-q) \rvert^{k-1}$, $k$ odd | 用 memory terms $q$ 抓 envelope/bias/thermal memory |
| **GMP** | Memory poly + cross terms $x(n-q) \cdot \lvert x(n-q-l) \rvert^{k-1}$ | 用 cross terms 抓更複雜的 memory interaction |
| **Volterra** | 完整 kernel expansion | 理論上 complete 但參數爆炸 |

**核心觀念**：
- **Odd-order basis** 抓的是 PA 直接的 in-band polynomial distortion (AM-AM 的 odd 成分 + AM-PM)
- **Memory terms** 抓的是 even order 透過 bias/thermal/reflection 間接回來的 memory effects
- **加 even basis** 到 memoryless DPD model 只會 over-fit noise — 因為 direct path 本來就是零

### 8.2 實務建議

1. **窄頻訊號 (LTE 20 MHz, 甚至 5G 100 MHz)**：memoryless odd polynomial 就很夠
2. **寬頻訊號 (5G 200 MHz+, mmWave)**：memory polynomial 或 GMP，odd-only basis
3. **Wideband 高效率 PA (Class-F, Doherty)**：需要 GMP 或 Volterra，因為 bias memory / harmonic termination 效應強
4. **永遠不要在 DPD basis 裡手動加 even order terms** — 這是初學者最容易犯的錯

---

## 9. 總結

### 9.1 兩個答案

**問：Even order 會產生 ACLR 嗎？**

- **Baseband polynomial 數學模型**：**會**。
  因為 $y_k = x \cdot |x|^{k-1}$ 對 even $k$ 含 $|x|^\text{odd}$ 這個 sqrt 因子，不是 polynomial，spectrum 有 infinite tail 洩漏到 ACLR band。模擬結果表 4.1 明確顯示 $k=2, 4, 6, 8$ 的 ACLR 數字非零。

- **實體 RF PA (從天線看出去)**：**不會**。
  因為 $y_\text{RF} = \sum a_k x_\text{RF}^k$ 裡 even $k$ 的項 frequency 成分在 $0, \pm 2\omega_c, \pm 4\omega_c, \ldots$，**完全沒有 $\pm \omega_c$ 的成分**。經過 output matching + TX filter + 天線 bandpass 後，even order 在 fundamental band 是 exactly zero。

### 9.2 每個 Order 的 Regrowth 位置 (baseband polynomial model)

```math
\text{支撐範圍}(y_k) \;=\; \left[-\tfrac{k \cdot \text{BW}}{2},\; +\tfrac{k \cdot \text{BW}}{2}\right]
```

對 $\text{BW} = 20$ MHz：

| $k$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 主要能量範圍 (MHz) | ±10 | ±20 | ±30 | ±40 | ±50 | ±60 | ±70 | ±80 | ±90 |
| Polynomial? | ✓ | ✗ | ✓ | ✗ | ✓ | ✗ | ✓ | ✗ | ✓ |
| Hard cutoff? | ✓ | ✗ (tail) | ✓ | ✗ (tail) | ✓ | ✗ (tail) | ✓ | ✗ (tail) | ✓ |

### 9.3 一句話總結

> **Odd order basis 是 $x, x^{\ast}$ 的純多項式，有 exact $k \cdot \text{BW}$ bandwidth；Even order basis 含 $\sqrt{x x^{\ast}}$ 是 non-polynomial，有 infinite tail。實體 RF PA 的 bandpass 路徑把 even order 從 fundamental band 清乾淨，所以 DPD 只需要 odd basis；even order 的影響改走 memory effect 那條路，由 memory polynomial 的 delay terms 負責捕捉。**

---

## 附錄 A: 相關檔案

| 檔案 | 用途 |
|:---|:---|
| `PA_poly_regrowth.m` | 完整 MATLAB script |
| `PA_poly_input.png` | Input signal 驗證 (spectrum + amplitude distribution) |
| `PA_poly_odd.png` | Figure 1 — Odd orders overlay ($k = 1, 3, 5, 7, 9$) |
| `PA_poly_all.png` | Figure 2 — All orders overlay ($k = 1 \ldots 9$) |
| `PA_poly_grid.png` | Figure 3 — 3×3 per-order grid |
| `PA_poly_regrowth_report.md` | **本文件** |

## 附錄 B: 名詞對照

| 縮寫 / 符號 | 意義 |
|:---|:---|
| PA | Power Amplifier |
| DPD | Digital Pre-Distortion |
| BW | Bandwidth |
| ACLR | Adjacent Channel Leakage Ratio |
| PAPR | Peak-to-Average Power Ratio |
| BPF | Bandpass Filter |
| $f_c$, $\omega_c$ | Carrier frequency (linear / angular) |
| $x(t)$ | Complex baseband signal |
| $x^{\ast}$ | Complex conjugate of $x$ |
| $\lvert x \rvert$ | Magnitude (envelope) of $x$, $= \sqrt{x\,x^{\ast}}$ |
| $y_k$ | $k$-th order basis: $y_k = x \cdot \lvert x \rvert^{k-1}$ |
| GMP | Generalized Memory Polynomial |
