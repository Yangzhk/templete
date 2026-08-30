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
  - [莫比乌斯反演](#莫比乌斯反演)
- [斯特林数](#斯特林数)
  - [定义与递推](#定义与递推)
- [经典组合数列](#经典组合数列)
  - [卡特兰数](#卡特兰数)
  - [错位排列](#错位排列)
  - [贝尔数](#贝尔数)
  - [分拆数](#分拆数)

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

### 莫比乌斯反演

**经典形式**：

$$g(n) = \sum_{d \mid n} f(d) \iff f(n) = \sum_{d \mid n} \mu(d) \cdot g\!\left(\frac{n}{d}\right) = \sum_{d \mid n} \mu\!\left(\frac{n}{d}\right) \cdot g(d)$$

**倍数形式**：

$$g(n) = \sum_{n \mid k \le N} f(k) \iff f(n) = \sum_{n \mid k \le N} \mu\!\left(\frac{k}{n}\right) g(k)$$

**代数理解**：$$g = f * 1 \iff f = g * \mu$$

---

## 斯特林数

### 定义与递推

第二类 $S(n, k)$：把 $n$ 个有标号元素分到 $k$ 个非空无标号集合的方案数。

$$S(n, k) = k \cdot S(n-1, k) + S(n-1, k-1),\quad S(0,0) = 1$$

第一类（无符号）$\begin{bmatrix} n \\ k \end{bmatrix}$：$n$ 元置换分解成 $k$ 个轮换的方案数。

$$\begin{bmatrix} n \\ k \end{bmatrix} = (n-1) \begin{bmatrix} n-1 \\ k \end{bmatrix} + \begin{bmatrix} n-1 \\ k-1 \end{bmatrix}$$

带符号 $s(n, k) = (-1)^{n-k} \begin{bmatrix} n \\ k \end{bmatrix}$。

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
