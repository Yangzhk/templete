- [生成函数](#生成函数)
  - [OGF / EGF 基础](#ogf--egf-基础)
  - [常用闭式](#常用闭式)
  - [运算解读](#运算解读)
- [组合恒等式](#组合恒等式)
  - [基本恒等式](#基本恒等式)
  - [反演公式](#反演公式)
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

---

## 组合恒等式

### 基本恒等式

$$\binom{n}{k} = \binom{n}{n-k},\qquad \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}\ (\text{Pascal})$$

吸收 / 提取：

$$k \binom{n}{k} = n \binom{n-1}{k-1},\qquad (n-k)\binom{n}{k} = n \binom{n-1}{k}$$

朱世杰（hockey-stick）：

$$\sum_{i=k}^{n} \binom{i}{k} = \binom{n+1}{k+1},\qquad \sum_{i=0}^{n} \binom{r+i}{i} = \binom{r+n+1}{n}$$

范德蒙德卷积：

$$\sum_{k} \binom{m}{k} \binom{n}{r-k} = \binom{m+n}{r}$$

行 / 列和：

$$\sum_{k=0}^{n} \binom{n}{k} = 2^n,\qquad \sum_{k=0}^{n} (-1)^k \binom{n}{k} = [n = 0]$$

$$\sum_k k \binom{n}{k} = n \cdot 2^{n-1},\qquad \sum_k k^2 \binom{n}{k} = n(n+1) \cdot 2^{n-2}$$

带平方：$\sum_k \binom{n}{k}^2 = \binom{2n}{n}$（范德蒙德 $m = n, r = n$）。

下降幂展开：$x^n = \sum_k S(n,k) \cdot x^{\underline{k}}$，反向 $x^{\underline{n}} = \sum_k s(n,k) x^k$（$s$ 是带符号第一类）。

### 反演公式

二项式反演（最常用形）：

$$f(n) = \sum_{k=0}^{n} \binom{n}{k} g(k) \iff g(n) = \sum_{k=0}^{n} (-1)^{n-k} \binom{n}{k} f(k)$$

子集反演（莫比乌斯）：

$$f(S) = \sum_{T \subseteq S} g(T) \iff g(S) = \sum_{T \subseteq S} (-1)^{|S \setminus T|} f(T)$$

斯特林反演（$S$ 第二类，$s$ 带符号第一类）：

$$f(n) = \sum_k S(n,k) g(k) \iff g(n) = \sum_k s(n,k) f(k)$$

子集和反演（高维前缀和的逆）：$f(S) = \sum_{T \supseteq S} g(T) \iff g(S) = \sum_{T \supseteq S} (-1)^{|T \setminus S|} f(T)$。这两个分别就是 `fwt_or` / `fwt_and` 的逆变换。

最大 / 至多反演：若 $f(n)$ 是"至多 $n$"的方案，$g(n)$ 是"恰好 $n$"，则

$$f(n) = \sum_{k \le n} \binom{n}{k} g(k),\quad g(n) = \sum_{k \le n} (-1)^{n-k} \binom{n}{k} f(k)$$

GCD 反演（$\mu$）：

$$\sum_{i=1}^{n} \sum_{j=1}^{m} [\gcd(i,j) = 1] = \sum_{d=1}^{\min(n,m)} \mu(d) \lfloor n/d \rfloor \lfloor m/d \rfloor$$

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
