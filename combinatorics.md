- [生成函数](#生成函数)
  - [OGF / EGF 基础](#ogf--egf-基础)
  - [常用闭式](#常用闭式)
  - [运算解读](#运算解读)
  - [PGF（概率生成函数）](#pgf概率生成函数)
- [组合恒等式](#组合恒等式)
  - [基本恒等式](#基本恒等式)
  - [反演公式](#反演公式)
  - [迪利克雷卷积](#迪利克雷卷积)
  - [莫比乌斯反演](#莫比乌斯反演)
- [斯特林数](#斯特林数)
  - [定义与递推](#定义与递推)
  - [第二类·行 (固定 n)](#第二类行-固定-n)
  - [第二类·列 (固定 k)](#第二类列-固定-k)
  - [第一类·行 (固定 n)](#第一类行-固定-n)
  - [第一类·列 (固定 k)](#第一类列-固定-k)
- [经典组合数列](#经典组合数列)
  - [卡特兰数](#卡特兰数)
  - [错位排列](#错位排列)
  - [贝尔数](#贝尔数)
  - [分拆数](#分拆数)
- [模板使用说明](#模板使用说明)
  - [任意模数 NTT (mul_any)](#任意模数-ntt-mul_any)
  - [多项式除法 (divmod)](#多项式除法-divmod)
  - [多项式平移 (shift)](#多项式平移-shift)
  - [多点求值 / 快速插值](#多点求值--快速插值)
  - [Bostan–Mori (linear_recur)](#bostanmori-linear_recur)
  - [连续点值拉格朗日 (lagrange_continuous)](#连续点值拉格朗日-lagrange_continuous)
  - [FWT 三种卷积](#fwt-三种卷积)
  - [子集卷积](#子集卷积)
  - [斯特林数 行 / 列 选哪个](#斯特林数-行--列-选哪个)


## 生成函数

### OGF / EGF 基础

普通生成函数：$A(x) = \sum_{n \ge 0} a_n x^n$
指数生成函数：$\hat A(x) = \sum_{n \ge 0} a_n \dfrac{x^n}{n!}$

| | OGF 适用 | EGF 适用 |
|---|---|---|
| 对象 | 无标号、可重复（多重集、组合） | 有标号、不可重复（排列、标号树） |
| 乘法含义 | $c_n = \sum_k a_k b_{n-k}$，无序拼接 | $c_n = \sum_k \binom{n}{k} a_k b_{n-k}$，标号合并 |

### 常用闭式

OGF：

$$\frac{1}{1-x} = \sum_{n \ge 0} x^n,\qquad \frac{1}{(1-x)^k} = \sum_{n \ge 0} \binom{n+k-1}{k-1} x^n$$

$$\frac{1}{1-cx} = \sum_{n \ge 0} c^n x^n,\qquad \sum_{n \ge k} \binom{n}{k} x^n = \frac{x^k}{(1-x)^{k+1}}$$

$$(1+x)^n = \sum_{k} \binom{n}{k} x^k,\qquad \frac{1 - \sqrt{1-4x}}{2x} = \sum_n C_n x^n\ (\text{卡特兰})$$

EGF：

$$e^x = \sum_n \frac{x^n}{n!},\quad \frac{1}{1-x} \text{ EGF} \to a_n = n!,\quad \ln \frac{1}{1-x} \text{ EGF} \to a_n = (n-1)!$$

$$\frac{e^x + e^{-x}}{2} = \sum_n \frac{x^{2n}}{(2n)!},\quad \frac{e^x - e^{-x}}{2} = \sum_n \frac{x^{2n+1}}{(2n+1)!}$$

### 运算解读

| 运算 | OGF 含义 | EGF 含义 |
|---|---|---|
| $A + B$ | 不相交并（互斥选一） | 同左 |
| $A \cdot B$ | 拼接（有序对的计数和） | 标号合并（含 $\binom{n}{k}$） |
| $A^k / k!$（EGF）或 $A^k$（OGF） | $k$ 段拼接 | $k$ 个标号子结构 |
| $\dfrac{1}{1 - A}$（要求 $A(0) = 0$） | 任意多段拼接 | 同左 |
| $\exp A$（EGF, $A(0)=0$） | — | 任意多个标号子结构（无序） |
| $\ln \dfrac{1}{1 - A}$（EGF） | — | 标号子结构按"环"组合 |
| $A(x) \cdot \dfrac{1}{1 - x}$ | 前缀和 $b_n = \sum_{k \le n} a_k$ | — |
| $A(B(x))$（$B(0) = 0$） | 替换 / 复合 | 同左 |

### PGF（概率生成函数）

对取非负整数值的随机变量 $X$，其 PGF 为：

$$G_X(z) = \mathbb E[z^X] = \sum_{k \ge 0} P(X = k) \cdot z^k$$

本质是概率质量函数的 OGF，系数非负且 $\sum p_k = 1$。

基本性质：

| 性质 | 公式 |
|---|---|
| 归一化 | $G_X(1) = 1$ |
| 期望 | $\mathbb E[X] = G_X'(1)$ |
| 方差 | $\operatorname{Var}[X] = G_X''(1) + G_X'(1) - [G_X'(1)]^2$ |
| $k$ 阶阶乘矩 | $\mathbb E[X^{\underline{k}}] = G_X^{(k)}(1)$ |

独立随机变量之和的 PGF 等于各自 PGF 的乘积：$G_{X+Y}(z) = G_X(z) \cdot G_Y(z)$。这一性质让 PGF 在处理「独立同分布随机变量之和」时非常有力——只需对单个变量的 PGF 取幂。

与 OGF 的关系：PGF 就是归一化的 OGF（$\sum a_n = 1$）。若 $A(x)$ 是某组合类的 OGF，则 $A(x) / A(1)$ 给出「均匀随机选取一个对象时，大小的分布」的 PGF。反过来，若已知 PGF，乘以总对象数即得计数 OGF。

常见随机变量的 PGF：

| 分布 | $P(X = k)$ | $G_X(z)$ |
|---|---|---|
| 伯努利 $B(p)$ | $p^k (1-p)^{1-k},\ k \in \{0,1\}$ | $1 - p + pz$ |
| 二项 $B(n, p)$ | $\binom{n}{k} p^k (1-p)^{n-k}$ | $(1 - p + pz)^n$ |
| 泊松 $\operatorname{Pois}(\lambda)$ | $e^{-\lambda} \lambda^k / k!$ | $e^{\lambda(z - 1)}$ |
| 几何 $\operatorname{Geom}(p)$ | $(1-p)^{k-1} p,\ k \ge 1$ | $\dfrac{pz}{1 - (1-p)z}$ |

竞赛中的典型用法：求某个随机变量期望 / 方差时，构造 PGF 后求导代入 $z = 1$，避免显式求和。若只关心期望，通常直接线性期望即可；PGF 更适合「$k$ 次操作后某值的分布」「随机过程停时」等场景。

---

## 组合恒等式

### 基本恒等式

#### 组合恒等式

**递推式**
$$\binom{n}{k} = \binom{n}{n-k},\qquad \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}\ (\text{Pascal})$$

**吸收 / 提取：**

$$k \binom{n}{k} = n \binom{n-1}{k-1},\qquad (n-k)\binom{n}{k} = n \binom{n-1}{k}$$

**朱世杰（hockey-stick）：**

$$\sum_{i=k}^{n} \binom{i}{k} = \binom{n+1}{k+1},\qquad \sum_{i=0}^{n} \binom{r+i}{i} = \binom{r+n+1}{n}$$

**范德蒙德卷积：**

$$\sum_{k} \binom{m}{k} \binom{n}{r-k} = \binom{m+n}{r}$$

**行 / 列和：**

$$\sum_{k=0}^{n} \binom{n}{k} = 2^n,\qquad \sum_{k=0}^{n} (-1)^k \binom{n}{k} = [n = 0]$$

$$\sum_k k \binom{n}{k} = n \cdot 2^{n-1},\qquad \sum_k k^2 \binom{n}{k} = n(n+1) \cdot 2^{n-2}$$

**带平方：**$\sum_k \binom{n}{k}^2 = \binom{2n}{n}$（范德蒙德 $m = n, r = n$）。

**下降幂展开**：$x^n = \sum_k S(n,k) \cdot x^{\underline{k}}$，反向 $x^{\underline{n}} = \sum_k s(n,k) x^k$（$s$ 是带符号第一类）。

#### Min-Max 容斥 (最值反演)

Min-Max 容斥用于在集合的“最大值”与“最小值”之间进行转换。它不仅适用于离散数值，**更重要的是在期望下依然严格成立**。

**1. 基础形式**
对于任意有限全集 $S$ 中的数值集合：
$$\max(S) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} \min(T)$$
$$\min(S) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} \max(T)$$

**核心应用（期望化最值）：**
在随机游走或抽卡问题中，求“全集都被覆盖的期望时间（即最后覆盖那个元素的期望时间，$\max$）”往往极其困难。但求“集合中任意一个元素被覆盖的期望时间（即最先覆盖的期望时间，$\min$）”通常可以直接得出。由于期望的线性性质：
$$E(\max(S)) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} E(\min(T))$$
*注：要求元素非负或存在合理的测度空间。*

**2. 扩展形式（$k$ 大 / $k$ 小）**
当需要求集合中第 $k$ 大或第 $k$ 小的元素（规定最大为第 $1$ 大）时：
$$\mathrm{kmax}(S) = \sum_{T \subseteq S, |T| \ge k} (-1)^{|T|-k} \binom{|T|-1}{k-1} \min(T)$$
$$\mathrm{kmin}(S) = \sum_{T \subseteq S, |T| \ge k} (-1)^{|T|-k} \binom{|T|-1}{k-1} \max(T)$$
*常用于“收集到 $n-k+1$ 种卡片（即第 $k$ 大覆盖时间）”的期望计算。*

**3. 常见算法实现与优化套路**
Min-Max 容斥的公式是指数级的（需要枚举所有子集），在竞赛中通常结合以下技术进行优化：

*   **FMT / 高维前缀和 (SOS DP)：**
    如果题目需要求出**所有**集合 $S$ 的 $E(\max(S))$，可以先计算出所有的 $E(\min(T))$，将式子视为：
    $$F(S) = \sum_{T \subseteq S} (-1)^{|T|-1} E(\min(T))$$
    令 $g(T) = (-1)^{|T|-1} E(\min(T))$，对其做一次 FWT-OR 的正变换即可在 $O(n \cdot 2^n)$ 内得到所有答案。

*   **背包 DP（处理 $n$ 较大但 $k$ 较小的情况）：**
    当全集大小 $n$ 较大（如 $n \le 1000$ 或更大），但元素权重/概率只有有限几种，且无需区分具体元素时。
    设 $dp[i][j][\dots]$ 表示前 $i$ 个元素中，选出的子集大小为 $j$，且其 $\min(T)$ 的某种贡献值为 $\dots$ 的方案数或容斥系数和。
    转移时直接将容斥系数 $(-1)^{|T|-k} \binom{|T|-1}{k-1}$ 压入 DP 状态中，避免最后枚举集合。
    对于基础 Min-Max 容斥，常常只需维护加入元素后容斥系数变号的特性：$dp[i][\dots] = dp[i-1][\dots] - dp[i-1][\dots \text{加上当前元素}]$。

*   **生成函数 (EGF) 合并：**
    当要求期望时间且每个元素的概率为 $p_i$（$\sum p_i = 1$），$E(\min(T)) = \frac{1}{\sum_{i \in T} p_i}$ 时，通常将 $\frac{1}{\sum p_i}$ 写成积分形式 $\int_0^\infty e^{-t \sum p_i} dt$。此时可以将 Min-Max 容斥的指数级求和转化为多个多项式（或 EGF）的乘积，最终在 $O(n \log^2 n)$ 或分治 NTT 复杂度下解决。

### 反演公式

#### 二项式反演（最常用形）：

$$f(n) = \sum_{k=0}^{n} \binom{n}{k} g(k) \iff g(n) = \sum_{k=0}^{n} (-1)^{n-k} \binom{n}{k} f(k)$$

#### 子集反演（莫比乌斯）：

$$f(S) = \sum_{T \subseteq S} g(T) \iff g(S) = \sum_{T \subseteq S} (-1)^{|S \setminus T|} f(T)$$

#### 斯特林反演（$S$ 第二类，$s$ 带符号第一类）：

$$f(n) = \sum_k S(n,k) g(k) \iff g(n) = \sum_k s(n,k) f(k)$$

#### 子集和反演（高维前缀和的逆）：$f(S) = \sum_{T \supseteq S} g(T) \iff g(S) = \sum_{T \supseteq S} (-1)^{|T \setminus S|} f(T)$。这两个分别就是 `fwt_or` / `fwt_and` 的逆变换。

#### 最大 / 至多反演：若 $f(n)$ 是"至多 $n$"的方案，$g(n)$ 是"恰好 $n$"，则

$$f(n) = \sum_{k \le n} \binom{n}{k} g(k),\quad g(n) = \sum_{k \le n} (-1)^{n-k} \binom{n}{k} f(k)$$

#### GCD 反演（$\mu$）：

$$\sum_{i=1}^{n} \sum_{j=1}^{m} [\gcd(i,j) = 1] = \sum_{d=1}^{\min(n,m)} \mu(d) \lfloor n/d \rfloor \lfloor m/d \rfloor$$

### 迪利克雷卷积

对算术函数 $f, g: \mathbb Z_{>0} \to \mathbb C$，定义迪利克雷卷积：

$$(f * g)(n) = \sum_{d \mid n} f(d) \cdot g\!\left(\frac{n}{d}\right) = \sum_{ab = n} f(a) \cdot g(b)$$

**代数结构**：迪利克雷卷积构成交换幺半群——满足交换律、结合律，单位元为 $\varepsilon(n) = [n = 1]$。若 $f(1) \ne 0$，则 $f$ 存在迪利克雷逆。

**常用算术函数表**：

| 函数 | 定义 | 含义 |
|---|---|---|
| $\varepsilon(n)$ | $[n = 1]$ | 卷积单位元 |
| $1(n)$ | $1$ | 全 1 函数 |
| $\operatorname{id}_k(n)$ | $n^k$ | 幂函数，$\operatorname{id}_1$ 简记为 $\operatorname{id}$ |
| $\mu(n)$ | $(-1)^k$ 若 $n$ 是 $k$ 个不同质数之积，否则 $0$ | 莫比乌斯函数 |
| $\varphi(n)$ | $\sum_{i=1}^n [\gcd(i, n) = 1]$ | 欧拉函数 |
| $d(n)$ | $\sum_{d \mid n} 1$ | 除数个数 |
| $\sigma_k(n)$ | $\sum_{d \mid n} d^k$ | 除数和 |

**核心卷积恒等式**：

$$\mu * 1 = \varepsilon \qquad\text{(莫比乌斯函数的定义性恒等式)}$$

$$\varphi * 1 = \operatorname{id} \qquad (\sum_{d \mid n} \varphi(d) = n)$$

$$1 * 1 = d \qquad \operatorname{id}_k * 1 = \sigma_k$$

$$\mu * \operatorname{id} = \varphi \qquad (\varphi(n) = \sum_{d \mid n} \mu(d) \cdot n/d)$$

$\mu * 1 = \varepsilon$ 是最核心的一条：它说 $\mu$ 是 $1$ 的迪利克雷逆，也是莫比乌斯反演的代数本质——$g = f * 1 \iff f = g * \mu$。

### 莫比乌斯反演

**经典形式**（除数形式）：

$$g(n) = \sum_{d \mid n} f(d) \iff f(n) = \sum_{d \mid n} \mu(d) \cdot g\!\left(\frac{n}{d}\right) = \sum_{d \mid n} \mu\!\left(\frac{n}{d}\right) \cdot g(d)$$

**倍数形式**（常用于 $\gcd$ 求和）：

$$g(n) = \sum_{n \mid k \le N} f(k) \iff f(n) = \sum_{n \mid k \le N} \mu\!\left(\frac{k}{n}\right) g(k)$$

**代数理解**：$g = f * 1 \iff f = g * \mu$。卷积形式比和式更容易记忆和推导，所有反演本质都是乘以 $1$ 的逆元 $\mu$。

---

#### 应用场景

**1. GCD 条件提取**（最常用）

核心手法：用 $\mu * 1 = \varepsilon$ 将 $[\gcd = 1]$ 展开为和式：

$$[\gcd(i, j) = 1] = \sum_{d \mid \gcd(i, j)} \mu(d) = \sum_{d \mid i,\ d \mid j} \mu(d)$$

推广到 $[\gcd(i, j) = g]$：

$$[\gcd(i, j) = g] = \sum_{d \mid \gcd(i, j) / g} \mu(d) = \sum_{gd \mid i,\ gd \mid j} \mu(d)$$

**经典题**：求 $\sum_{i=1}^n \sum_{j=1}^m [\gcd(i, j) = 1]$。

$$\begin{aligned}
\sum_{i=1}^n \sum_{j=1}^m [\gcd(i, j) = 1]
&= \sum_{i=1}^n \sum_{j=1}^m \sum_{d \mid i,\ d \mid j} \mu(d) \\
&= \sum_{d=1}^{\min(n,m)} \mu(d) \sum_{i=1}^{\lfloor n/d \rfloor} \sum_{j=1}^{\lfloor m/d \rfloor} 1 \\
&= \sum_{d=1}^{\min(n,m)} \mu(d) \left\lfloor \frac{n}{d} \right\rfloor \left\lfloor \frac{m}{d} \right\rfloor
\end{aligned}$$

配上整除分块，$O(\sqrt n)$ 单次查询。

**2. 求 $\gcd$ 的和 / LCM 的和**

$$\sum_{i=1}^n \sum_{j=1}^m \gcd(i, j) = \sum_{d=1}^{\min(n,m)} d \cdot f\!\left(\left\lfloor \frac{n}{d} \right\rfloor, \left\lfloor \frac{m}{d} \right\rfloor\right)$$

其中 $f(n, m)$ 是 $\gcd = 1$ 的对数（上一条），枚举 $d = \gcd$ 后内层变为互质计数。

$$\sum_{i=1}^n \sum_{j=1}^m \operatorname{lcm}(i, j) = \sum_{d=1}^{\min(n,m)} \frac{1}{d} \sum_{i=1}^{\lfloor n/d \rfloor} \sum_{j=1}^{\lfloor m/d \rfloor} ij \cdot [\gcd(i, j) = 1] \cdot d^2$$

本质仍是提取 $\gcd$ 条件后处理。

**3. 乘积形式反演**

$$g(n) = \prod_{d \mid n} f(d) \iff f(n) = \prod_{d \mid n} g(d)^{\mu(n/d)}$$

出现在需要计算「某个函数值的积」且满足乘性关系时，如圆排列计数、Lyndon 词计数（用 Moreau 项链恒等式等价于 $\mu$ 反演的乘积形）。

**4. 周期序列与本质不同计数**

长度为 $n$ 的环 / 项链计数中，莫比乌斯反演用于去重：

$$\text{本质不同的周期为 } n \text{ 的序列数} = \frac{1}{n} \sum_{d \mid n} \mu(d) \cdot k^{n/d}$$

这就是 Burnside / Polya 定理中循环群的贡献，本质是倍数反演。

---

#### 使用技巧

**技巧 1：先写成卷积再翻译回和式。** 许多反演题容易写错下标。记住 $g = f * 1 \iff f = g * \mu$ 这一条，展开时 $d \mid n$ 和 $n/d$ 的位置就不会搞混。

**技巧 2：$\mu$ 预处理 + 整除分块。** 用线性筛 $O(N)$ 预处理 $\mu$ 的前缀和，查询时整除分块 $O(\sqrt n)$。常见组合：预处理 $\mu$ 前缀和 → 整除分块 → 内层 $O(1)$ 计算 $\lfloor n/d \rfloor \lfloor m/d \rfloor$。

**技巧 3：双层枚举交换顺序。** 套路：

$$\sum_{i=1}^n \sum_{j=1}^m f(\gcd(i, j)) = \sum_{d=1}^{\min(n,m)} f(d) \sum_{i=1}^n \sum_{j=1}^m [\gcd(i, j) = d]$$

令 $i = di', j = dj'$，内层变为 $[\gcd(i', j') = 1]$，再用技巧 1 展开 $\mu$。最后换回外层 $d$ 求和。

**技巧 4：$\varphi$ 替代 $\mu$。** 若只需要 $\gcd = 1$，有时用 $\varphi$ 更直接：

$$\sum_{i=1}^n \sum_{j=1}^n [\gcd(i, j) = 1] = 2 \sum_{i=1}^n \varphi(i) - 1$$

多出来的 $-1$ 是去掉 $(1,1)$ 重复计数。当 $n = m$ 且只需 $\gcd = 1$ 时，欧拉函数前缀和比 $\mu$ 整除分块常数更小。

**技巧 5：迪利克雷前缀和。** 对于 $g = f * 1$ 即 $g(n) = \sum_{d \mid n} f(d)$，给定 $f$ 求所有 $g(1..N)$ 可在 $O(N \log \log N)$ 内完成（用质数倍数枚举，类似埃筛）。这是除数形式的批量计算方式。反过来，已知 $g$ 求 $f$ 同理。

**技巧 6：$[\gcd = 1]$ 以外的情况。** 要解 $[\gcd = k]$ 时，先转化为 $\gcd = 1$：

$$[\gcd(i_1, \dots, i_m) = k] = [\gcd(i_1/k, \dots, i_m/k) = 1] \cdot [k \mid i_1] \cdots [k \mid i_m]$$

然后对新变量 $i_j' = i_j / k$ 做 $\mu$ 展开。

---

## 斯特林数

### 定义与递推

第二类 $S(n, k)$：把 $n$ 个有标号元素分到 $k$ 个非空无标号集合的方案数。

$$S(n, k) = k \cdot S(n-1, k) + S(n-1, k-1),\quad S(0,0) = 1$$

第一类（无符号）$\begin{bmatrix} n \\ k \end{bmatrix}$：$n$ 元置换分解成 $k$ 个轮换的方案数。

$$\begin{bmatrix} n \\ k \end{bmatrix} = (n-1) \begin{bmatrix} n-1 \\ k \end{bmatrix} + \begin{bmatrix} n-1 \\ k-1 \end{bmatrix}$$

带符号 $s(n, k) = (-1)^{n-k} \begin{bmatrix} n \\ k \end{bmatrix}$。

### 第二类·行 (固定 n)

$$S(n, k) = \frac{1}{k!} \sum_{i=0}^{k} (-1)^i \binom{k}{i} (k-i)^n
       = \sum_{i=0}^{k} \frac{(-1)^i}{i!} \cdot \frac{(k-i)^n}{(k-i)!}$$

是卷积 $a_i = (-1)^i / i!$ 与 $b_j = j^n / j!$，$O(n \log n)$。

```cpp
// S(n, 0..n)
static vector<int> stirling2_row(int n) {
    vector<int> fac(n + 1), ifac(n + 1);
    fac[0] = 1;
    for (int i = 1; i <= n; i++) fac[i] = 1ll * fac[i-1] * i % Poly::mod;
    ifac[n] = Poly::inv_mod(fac[n]);
    for (int i = n - 1; i >= 0; i--) ifac[i] = 1ll * ifac[i+1] * (i+1) % Poly::mod;
    vector<int> A(n + 1), B(n + 1);
    for (int i = 0; i <= n; i++) {
        A[i] = (i & 1) ? (Poly::mod - ifac[i]) % Poly::mod : ifac[i];
        B[i] = 1ll * Poly::pw(i, n) * ifac[i] % Poly::mod;
    }
    auto C = Poly::mul_p(A, B, Poly::mod);
    C.resize(n + 1);
    return C;
}
```

### 第二类·列 (固定 k)

EGF：$\sum_{n} S(n, k) \dfrac{x^n}{n!} = \dfrac{(e^x - 1)^k}{k!}$。$O(N \log N)$。

```cpp
// S(0..N, k)
static vector<int> stirling2_col(int N, int k) {
    vector<int> fac(N + 1), ifac(N + 1);
    fac[0] = 1;
    for (int i = 1; i <= N; i++) fac[i] = 1ll * fac[i-1] * i % Poly::mod;
    ifac[N] = Poly::inv_mod(fac[N]);
    for (int i = N - 1; i >= 0; i--) ifac[i] = 1ll * ifac[i+1] * (i+1) % Poly::mod;
    vector<int> ex(N + 1, 0);
    for (int i = 1; i <= N; i++) ex[i] = ifac[i];   // (e^x - 1)
    Poly p = Poly(move(ex)).pow(N + 1, k);          // (e^x - 1)^k
    int ifk = ifac[k];
    vector<int> r(N + 1);
    for (int i = 0; i <= N; i++)
        r[i] = 1ll * p[i] * fac[i] % Poly::mod * ifk % Poly::mod;
    return r;
}
```

### 第一类·行 (固定 n)

OGF：$\sum_k \begin{bmatrix} n \\ k \end{bmatrix} x^k = x(x+1)(x+2) \cdots (x+n-1) =: f_n(x)$。

倍增：$f_{2m}(x) = f_m(x) \cdot f_m(x + m)$（用 [shift](poly.md#多项式平移-fxc)），单乘 $(x + c)$ 处理奇数位。$O(n \log n)$。

```cpp
// |s(n, 0..n)|
static vector<int> stirling1_row(int n) {
    if (n == 0) return {1};
    Poly cur({0, 1});
    int cnt = 1;
    int hi = 31 - __builtin_clz(n);
    for (int i = hi - 1; i >= 0; i--) {
        cur = cur * cur.shift(cnt);
        cnt <<= 1;
        if ((n >> i) & 1) {
            cur.a.push_back(0);
            for (int j = (int)cur.size() - 1; j >= 1; j--)
                cur.a[j] = (cur.a[j-1] + 1ll * cnt * cur.a[j]) % Poly::mod;
            cur.a[0] = 1ll * cnt * cur.a[0] % Poly::mod;
            cnt += 1;
        }
    }
    return cur.a;
}
```

### 第一类·列 (固定 k)

EGF：$\sum_n \begin{bmatrix} n \\ k \end{bmatrix} \dfrac{x^n}{n!} = \dfrac{1}{k!} \left( \ln \dfrac{1}{1-x} \right)^k$。$O(N \log N)$。

```cpp
// |s(0..N, k)|
static vector<int> stirling1_col(int N, int k) {
    vector<int> iv(N + 2, 0); iv[1] = 1;
    for (int i = 2; i <= N; i++)
        iv[i] = 1ll * (Poly::mod - Poly::mod / i) * iv[Poly::mod % i] % Poly::mod;
    vector<int> v(N + 1, 0);
    for (int i = 1; i <= N; i++) v[i] = iv[i];      // -ln(1-x) = sum x^i / i
    Poly p = Poly(move(v)).pow(N + 1, k);
    vector<int> fac(N + 1), ifac(N + 1);
    fac[0] = 1;
    for (int i = 1; i <= N; i++) fac[i] = 1ll * fac[i-1] * i % Poly::mod;
    ifac[N] = Poly::inv_mod(fac[N]);
    for (int i = N - 1; i >= 0; i--) ifac[i] = 1ll * ifac[i+1] * (i+1) % Poly::mod;
    int ifk = ifac[k];
    vector<int> r(N + 1);
    for (int i = 0; i <= N; i++)
        r[i] = 1ll * p[i] * fac[i] % Poly::mod * ifk % Poly::mod;
    return r;
}
```

---

## 经典组合数列

### 卡特兰数

$$C_n = \frac{1}{n+1} \binom{2n}{n} = \binom{2n}{n} - \binom{2n}{n+1},\quad C_{n+1} = \sum_{i=0}^{n} C_i C_{n-i}$$

OGF：$C(x) = \dfrac{1 - \sqrt{1 - 4x}}{2x}$。常见模型：合法括号串、二叉树、不下降路径。

### 错位排列

$$D_n = n D_{n-1} + (-1)^n,\quad D_n = n! \sum_{k=0}^{n} \frac{(-1)^k}{k!}$$

EGF：$\dfrac{e^{-x}}{1-x}$。容斥经典：选 $k$ 个固定点的方案是 $\binom{n}{k} D_{n-k}$。

### 贝尔数

$$B_n = \sum_{k=0}^{n} S(n, k),\qquad B_{n+1} = \sum_{k=0}^{n} \binom{n}{k} B_k$$

EGF：$\hat B(x) = e^{e^x - 1}$，可用 `Poly::exp` 一次算出 $B_0,\dots,B_N$。

### 分拆数

$p(n)$：把 $n$ 拆为不计顺序的正整数之和的方案数。

OGF：$P(x) = \prod_{k \ge 1} \dfrac{1}{1 - x^k}$。

Euler 五边形定理给出 $O(n \sqrt n)$ 递推：

$$p(n) = \sum_{k \ge 1} (-1)^{k+1} \left[ p(n - g_k) + p(n - g_k') \right],\quad g_k = \frac{k(3k-1)}{2},\ g_k' = \frac{k(3k+1)}{2}$$

---

## 模板使用说明

### 任意模数 NTT (mul_any)

`Poly::mul_any(a, b, p)` 处理 $p$ 不是 NTT 友好质数（如 $10^9 + 7$）的情况。内部用三模数 NTT + CRT 直接对 $p$ 取模，输入系数无大小限制（只要每项 $< 10^9$ 即可）。

何时用：

- 题目模数是 $10^9 + 7$、$10^9 + 9$ 等非 NTT 模 → 用 `mul_any`
- 题目模数是 $998244353$ → 直接 `Poly::operator*`，更快
- 需要精确整数（不取模）且结果 $< 2^{63}$ → 用 `mul_exact`

注意：常数比单模 NTT 大 3 倍以上，能用就尽量避免。

### 多项式除法 (divmod)

主要用途：

1. **特征多项式取模**：在 $\mathbb F_p[x] / Q(x)$ 中工作时（线性递推、矩阵幂的极小多项式），频繁需要 $\bmod Q$。
2. **拉格朗日插值的辅助**：在 [快速插值](#多点求值--快速插值) 内部递归求 $f \bmod (x - x_i)$。
3. **GCD / 扩展欧几里得**：求 $\gcd(F, G)$ 时反复 `divmod`，但纯多项式 GCD 用得少。

输入要求 $G$ 末项非零，否则需要先 trim 高次零系数。

### 多项式平移 (shift)

`f.shift(c)` 返回 $f(x + c)$。常见用途：

- 已知 $f(0), f(1), \dots, f(n)$，求 $f(c), f(c+1), \dots, f(c+n)$ → 把 $f$ 拉格朗日插值得到再 shift 不太常用，更常见是用 [连续点值拉格朗日](#连续点值拉格朗日-lagrange_continuous) 单点查询。
- 第一类斯特林行的倍增，见上。
- 处理"把题目转化为 $x = y + c$"的形式变量替换。

### 多点求值 / 快速插值

| 用途 | 方法 |
|---|---|
| 给定多项式 $f$，求 $f(x_1),\dots,f(x_n)$ | `multi_eval` |
| 给定 $n$ 个点对，构造度 $< n$ 的 $f$ | `fast_interp` |
| 点数 $n \le 5000$ 左右 | 用朴素 [`lagrange`](poly.md#一般插值)，常数小 |
| 点数 $n \ge 10^5$ | 用 `fast_interp`，$O(n \log^2 n)$ |

注意 `fast_interp` 要求 $x_i$ 互异。

### Bostan–Mori (linear_recur)

求 $a_n$ 满足 $a_n = \sum_{i=1}^{k} c_i a_{n-i}$ 的第 $n$ 项。

输入约定：`a = [a_0, ..., a_{k-1}]`，`c = [c_1, c_2, ..., c_k]`（注意是 1-based 系数，传 0-based 数组），`n` 是查询位置（可达 $10^{18}$）。

常见误区：

- 与矩阵快速幂 $O(k^3 \log n)$ 比较，本算法 $O(k \log k \log n)$，$k \ge 30$ 起明显占优。
- 若递推不齐次（带常数项），先把数列减去其稳态使齐次化，或扩展为 $k+1$ 维。
- 对未知递推的数列，先用 Berlekamp–Massey 求出 $c$，再用本模板。

### 连续点值拉格朗日 (lagrange_continuous)

已知 $f(1), f(2), \dots, f(n)$（$f$ 是度 $\le n - 1$ 的多项式），求 $f(k)$。

适用场景：

- 自然数幂和 $\sum_{i=1}^{N} i^d$，关于 $N$ 是 $d+1$ 次多项式。先算前 $d + 2$ 项暴力求和，再 lagrange 单点求 $f(N)$。
- 任何"答案是 $n$ 的多项式但 $n$ 极大"的题目。

$k$ 接受 `long long`，会自动模归约。若需要多个 $k$，建议先 [`fast_interp`](#多点求值--快速插值) 还原系数再 `eval`。

### FWT 三种卷积

| 卷积 | 触发条件 | 典型题目 |
|---|---|---|
| OR | $h_S = \sum_{T \cup U = S}$ | 子集和 / 至少覆盖某些位 |
| AND | $h_S = \sum_{T \cap U = S}$ | 超集和 / 必须不含某些位 |
| XOR | $h_S = \sum_{T \oplus U = S}$ | 异或 DP（线性基题、Nim 游戏分析） |

OR 的正变换就是子集和（高维前缀和），AND 是超集和。XOR 是 Hadamard 变换。

数组长度必须是 $2^n$，下标范围 $[0, 2^n)$。多个数组同次变换可以共用 FWT 调用。

### 子集卷积

$h_S = \sum_{T \subseteq S} f_T \cdot g_{S \setminus T}$。比 OR 卷积多了一个**不相交**约束。

常见用途：

- 集合 DP 中"用一组互不相交的元素拼成全集"的转移。
- 染色 / 划分类问题中，按 popcount 分层处理。

代价：内存 $O(2^n n)$，时间 $O(2^n n^2)$。$n \le 20$ 范围内常用。

### 斯特林数 行 / 列 选哪个

| 需求 | 用 |
|---|---|
| 一行 $S(n, *)$ 或 $\begin{bmatrix} n \\ * \end{bmatrix}$，$n$ 单点 | 行版 |
| 一列 $S(*, k)$，$k$ 固定遍历 $n$ | 列版 |
| 单点 $S(n, k)$，$n, k$ 均 $\le 5000$ | $O(nk)$ 递推表更简单 |
| $n$ 极大但 $k$ 小（如 $k \le 100$） | 列版的 EGF 截到 $n$ 即可，注意此时 $n$ 受多项式长度限制；$n$ 真的极大时可能要"线性递推 + Bostan–Mori" |
