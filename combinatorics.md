- [生成函数](#生成函数)
  - [常用闭式](#常用闭式)
- [组合恒等式](#组合恒等式)
  - [基本恒等式](#基本恒等式)
    - [组合恒等式](#组合恒等式-1)
    - [Min-Max 容斥](#min-max-容斥)
  - [反演公式](#反演公式)
    - [二项式反演（最常用形）：](#二项式反演最常用形)
    - [子集反演（莫比乌斯）：](#子集反演莫比乌斯)
    - [斯特林反演：](#斯特林反演)
    - [子集和反演：](#子集和反演)
    - [最大 / 至多反演：](#最大--至多反演)
    - [GCD 反演：](#gcd-反演)
  - [迪利克雷卷积](#迪利克雷卷积)
  - [莫比乌斯反演](#莫比乌斯反演)
    - [应用场景](#应用场景)
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

## 生成函数

### 常用闭式

**OGF：**

$$\frac{1}{1-x} = \sum_{n \ge 0} x^n,\qquad \frac{1}{(1-x)^k} = \sum_{n \ge 0} \binom{n+k-1}{k-1} x^n$$

$$\frac{1}{1-cx} = \sum_{n \ge 0} c^n x^n,\qquad \sum_{n \ge k} \binom{n}{k} x^n = \frac{x^k}{(1-x)^{k+1}}$$

$$(1+x)^n = \sum_{k} \binom{n}{k} x^k,\qquad \frac{1 - \sqrt{1-4x}}{2x} = \sum_n C_n x^n$$

**EGF：**

$$e^x = \sum_n \frac{x^n}{n!},\quad \frac{1}{1-x} \text{ EGF} \to a_n = n!,\quad \ln \frac{1}{1-x} \text{ EGF} \to a_n = (n-1)!$$

$$\frac{e^x + e^{-x}}{2} = \sum_n \frac{x^{2n}}{(2n)!},\quad \frac{e^x - e^{-x}}{2} = \sum_n \frac{x^{2n+1}}{(2n+1)!}$$


## 组合恒等式

### 基本恒等式

#### 组合恒等式

**递推式**
$$\binom{n}{k} = \binom{n}{n-k},\qquad \binom{n}{k} = \binom{n-1}{k} + \binom{n-1}{k-1}\ (\text{Pascal})$$

**吸收 / 提取：**

$$k \binom{n}{k} = n \binom{n-1}{k-1},\qquad (n-k)\binom{n}{k} = n \binom{n-1}{k}$$

**朱世杰恒等式：**

$$\sum_{i=k}^{n} \binom{i}{k} = \binom{n+1}{k+1},\qquad \sum_{i=0}^{n} \binom{r+i}{i} = \binom{r+n+1}{n}$$

**范德蒙德卷积：**

$$\sum_{k} \binom{m}{k} \binom{n}{r-k} = \binom{m+n}{r}$$

**行 / 列和：**

$$\sum_{k=0}^{n} \binom{n}{k} = 2^n,\qquad \sum_{k=0}^{n} (-1)^k \binom{n}{k} = [n = 0]$$

$$\sum_k k \binom{n}{k} = n \cdot 2^{n-1},\qquad \sum_k k^2 \binom{n}{k} = n(n+1) \cdot 2^{n-2}$$

**带平方：**

$$\sum_k \binom{n}{k}^2 = \binom{2n}{n}$$

**下降幂展开:**

$$x^n = \sum_k S(n,k) \cdot x^{\underline{k}}$$

$$x^{\underline{n}} = \sum_k s(n,k) x^k$$

#### Min-Max 容斥

Min-Max 容斥用于在集合的“最大值”与“最小值”之间进行转换。它不仅适用于离散数值，**更重要的是在期望下依然严格成立**。

**1. 基础形式**

对于任意有限全集 $S$ 中的数值集合：
$$\max(S) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} \min(T)$$
$$\min(S) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} \max(T)$$

**核心应用：**

在随机游走或抽卡问题中，求“全集都被覆盖的期望时间（即最后覆盖那个元素的期望时间，$\max$）”往往极其困难。但求“集合中任意一个元素被覆盖的期望时间（即最先覆盖的期望时间，$\min$）”通常可以直接得出。由于期望的线性性质：
$$E(\max(S)) = \sum_{T \subseteq S, T \neq \emptyset} (-1)^{|T|-1} E(\min(T))$$

**2. 扩展形式**

$$\mathrm{kmax}(S) = \sum_{T \subseteq S, |T| \ge k} (-1)^{|T|-k} \binom{|T|-1}{k-1} \min(T)$$
$$\mathrm{kmin}(S) = \sum_{T \subseteq S, |T| \ge k} (-1)^{|T|-k} \binom{|T|-1}{k-1} \max(T)$$


### 反演公式

#### 二项式反演（最常用形）：

$$f(n) = \sum_{k=0}^{n} \binom{n}{k} g(k) \iff g(n) = \sum_{k=0}^{n} (-1)^{n-k} \binom{n}{k} f(k)$$

#### 子集反演（莫比乌斯）：

$$f(S) = \sum_{T \subseteq S} g(T) \iff g(S) = \sum_{T \subseteq S} (-1)^{|S \setminus T|} f(T)$$

#### 斯特林反演：

$$f(n) = \sum_k S(n,k) g(k) \iff g(n) = \sum_k s(n,k) f(k)$$

#### 子集和反演：

$$f(S) = \sum_{T \supseteq S} g(T) \iff g(S) = \sum_{T \supseteq S} (-1)^{|T \setminus S|} f(T)$$

#### 最大 / 至多反演：

若 $f(n)$ 是"至多 $n$"的方案，$g(n)$ 是"恰好 $n$"，则

$$f(n) = \sum_{k \le n} \binom{n}{k} g(k),\quad g(n) = \sum_{k \le n} (-1)^{n-k} \binom{n}{k} f(k)$$

#### GCD 反演：

$$\sum_{i=1}^{n} \sum_{j=1}^{m} [\gcd(i,j) = 1] = \sum_{d=1}^{\min(n,m)} \mu(d) \lfloor n/d \rfloor \lfloor m/d \rfloor$$

### 迪利克雷卷积

对算术函数 $f, g: \mathbb Z_{>0} \to \mathbb C$，定义迪利克雷卷积：

$$(f * g)(n) = \sum_{d \mid n} f(d) \cdot g\!\left(\frac{n}{d}\right) = \sum_{ab = n} f(a) \cdot g(b)$$

**代数结构**：迪利克雷卷积构成交换幺半群——满足交换律、结合律，单位元为 $\varepsilon(n) = [n = 1]$。若 $f(1) \ne 0$，则 $f$ 存在迪利克雷逆。

**核心卷积恒等式**：

$$\mu * 1 = \varepsilon \qquad$$

$$\varphi * 1 = \operatorname{id} \qquad (\sum_{d \mid n} \varphi(d) = n)$$

$$1 * 1 = d \qquad \operatorname{id}_k * 1 = \sigma_k$$

$$\mu * \operatorname{id} = \varphi \qquad (\varphi(n) = \sum_{d \mid n} \mu(d) \cdot n/d)$$

$$g = f * 1 \iff f = g * \mu$$

### 莫比乌斯反演

**经典形式**：

$$g(n) = \sum_{d \mid n} f(d) \iff f(n) = \sum_{d \mid n} \mu(d) \cdot g\!\left(\frac{n}{d}\right) = \sum_{d \mid n} \mu\!\left(\frac{n}{d}\right) \cdot g(d)$$

**倍数形式**：

$$g(n) = \sum_{n \mid k \le N} f(k) \iff f(n) = \sum_{n \mid k \le N} \mu\!\left(\frac{k}{n}\right) g(k)$$

**代数理解**：$$g = f * 1 \iff f = g * \mu$$

---

#### 应用场景

**1. GCD 条件提取**

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

是 $a_i = (-1)^i / i!$ 与 $b_j = j^n / j!$ 的卷积，复杂度 $O(n \log n)$。

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

**EGF：**

$$\sum_{n} S(n, k) \dfrac{x^n}{n!} = \dfrac{(e^x - 1)^k}{k!}$$
复杂度 $O(N \log N)$

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

**OGF：**

$$\sum_k \begin{bmatrix} n \\ k \end{bmatrix} x^k = x(x+1)(x+2) \cdots (x+n-1) =: f_n(x)$$

**倍增:**

$$f_{2m}(x) = f_m(x) \cdot f_m(x + m)$$
单乘 $(x + c)$ 处理奇数位。复杂度 $O(n \log n)$。

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

**EGF：**

$$\sum_n \begin{bmatrix} n \\ k \end{bmatrix} \dfrac{x^n}{n!} = \dfrac{1}{k!} \left( \ln \dfrac{1}{1-x} \right)^k$$

复杂度 $O(N \log N)$

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

**OGF：**

$$C(x) = \dfrac{1 - \sqrt{1 - 4x}}{2x}$$

### 错位排列

$$D_n = n D_{n-1} + (-1)^n,\quad D_n = n! \sum_{k=0}^{n} \frac{(-1)^k}{k!}$$

**EGF：**
$$\dfrac{e^{-x}}{1-x}$$
举例：选 $k$ 个固定点的方案是 $\binom{n}{k} D_{n-k}$。

### 贝尔数

$$B_n = \sum_{k=0}^{n} S(n, k),\qquad B_{n+1} = \sum_{k=0}^{n} \binom{n}{k} B_k$$

**EGF**：
$$\hat B(x) = e^{e^x - 1}$$
可用 `Poly::exp` 一次算出 $B_0,\dots,B_N$

### 分拆数

**$p(n)$：** 把 $n$ 拆为不计顺序的正整数之和的方案数。

**OGF：**

$$P(x) = \prod_{k \ge 1} \dfrac{1}{1 - x^k}$$

Euler 五边形定理给出 $O(n \sqrt n)$ 递推：

$$p(n) = \sum_{k \ge 1} (-1)^{k+1} \left[ p(n - g_k) + p(n - g_k') \right]$$

$$\quad g_k = \frac{k(3k-1)}{2}   \ g_k' = \frac{k(3k+1)}{2}$$
