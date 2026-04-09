- [模意义下运算类](#模意义下运算类)
- [矩阵运算类](#矩阵运算类)
- [莫队](#莫队)
- [字符串哈希](#字符串哈希)
- [后缀数组](#后缀数组)
  - [后缀数组（Suffix Array）能解决的问题](#后缀数组suffix-array能解决的问题)
  - [一、基础能力](#一基础能力)
    - [1. 子串排序 / 字典序比较](#1-子串排序--字典序比较)
    - [2. 快速判断子串是否出现](#2-快速判断子串是否出现)
    - [3. 统计子串出现次数](#3-统计子串出现次数)
  - [二、LCP（最长公共前缀）相关问题](#二lcp最长公共前缀相关问题)
    - [4. 任意两个后缀的 LCP](#4-任意两个后缀的-lcp)
    - [5. 最长重复子串](#5-最长重复子串)
    - [6. 最长公共子串（两个字符串）](#6-最长公共子串两个字符串)
  - [三、子串统计问题](#三子串统计问题)
    - [8. 不同子串个数](#8-不同子串个数)
    - [9. 第 $k$ 小子串](#9-第-k-小子串)
    - [10. 子串去重统计](#10-子串去重统计)
  - [四、重复结构与模式问题](#四重复结构与模式问题)
    - [11. 最长重复且不重叠子串](#11-最长重复且不重叠子串)
    - [12. 出现至少 $k$ 次的最长子串](#12-出现至少-k-次的最长子串)
    - [13. 字符串周期问题](#13-字符串周期问题)
  - [五、区间与排名问题](#五区间与排名问题)
    - [14. 子串排名](#14-子串排名)
    - [15. 区间不同子串数](#15-区间不同子串数)
- [kmp 前缀函数 (前缀的最长 border)](#kmp-前缀函数-前缀的最长-border)
- [xor 线性基](#xor-线性基)
- [异或线性基（XOR Basis）](#异或线性基xor-basis)
    - [数学形式](#数学形式)
  - [维度与性质](#维度与性质)
  - [能解决的问题类型](#能解决的问题类型)
    - [1. 最大 / 最小异或值](#1-最大--最小异或值)
    - [2. 判断可表示性](#2-判断可表示性)
    - [3. 统计不同异或结果个数](#3-统计不同异或结果个数)
    - [4. 第 $k$ 小异或值](#4-第-k-小异或值)
    - [5. 子集异或相关计数问题](#5-子集异或相关计数问题)
    - [6. 图论问题](#6-图论问题)
    - [7. 区间 / 前缀问题](#7-区间--前缀问题)
  - [总结](#总结)
- [mobius 函数 (反演)](#mobius-函数-反演)
  - [莫比乌斯函数（Möbius Function）](#莫比乌斯函数möbius-function)
  - [莫比乌斯反演（Möbius Inversion）](#莫比乌斯反演möbius-inversion)
    - [形式一（约数形式）](#形式一约数形式)
    - [形式二（对称形式）](#形式二对称形式)
    - [常见等价写法](#常见等价写法)
- [exgcd 求 线性同余方程](#exgcd-求-线性同余方程)
- [中国剩余定理](#中国剩余定理)
- [三模数 ntt](#三模数-ntt)
- [多项式乘法 (ntt)](#多项式乘法-ntt)
- [三维点类](#三维点类)
- [最大流](#最大流)
- [zkw 费用流](#zkw-费用流)
- [计算几何 (dls version)](#计算几何-dls-version)

### 模意义下运算类
```
#include <iostream>
#include <algorithm>
#include <cassert>

constexpr int P = 1e9 + 7;
using i64 = long long;
// assume -P <= x < 2P
int norm(int x) {
    if (x < 0) {
        x += P;
    }
    if (x >= P) {
        x -= P;
    }
    return x;
}
template<class T>
T power(T a, i64 b) {
    T res = 1;
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
}
struct Z {
    int x;
    Z(int x = 0) : x(norm(x)) {}
    Z(i64 x) : x(norm(x % P)) {}
    int val() const {
        return x;
    }
    Z operator-() const {
        return Z(norm(P - x));
    }
    Z inv() const {
        assert(x != 0);
        return power(*this, P - 2);
    }
    Z &operator*=(const Z &rhs) {
        x = i64(x) * rhs.x % P;
        return *this;
    }
    Z &operator+=(const Z &rhs) {
        x = norm(x + rhs.x);
        return *this;
    }
    Z &operator-=(const Z &rhs) {
        x = norm(x - rhs.x);
        return *this;
    }
    Z &operator/=(const Z &rhs) {
        return *this *= rhs.inv();
    }
    friend Z operator*(const Z &lhs, const Z &rhs) {
        Z res = lhs;
        res *= rhs;
        return res;
    }
    friend Z operator+(const Z &lhs, const Z &rhs) {
        Z res = lhs;
        res += rhs;
        return res;
    }
    friend Z operator-(const Z &lhs, const Z &rhs) {
        Z res = lhs;
        res -= rhs;
        return res;
    }
    friend Z operator/(const Z &lhs, const Z &rhs) {
        Z res = lhs;
        res /= rhs;
        return res;
    }
    friend std::istream &operator>>(std::istream &is, Z &a) {
        i64 v;
        is >> v;
        a = Z(v);
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, const Z &a) {
        return os << a.val();
    }
};
 
Z fac[N], inv[N], pow2[N];
 
Z qpow(Z a, int v = P - 2){
    Z res = 1;
    for(; v; a *= a, v >>= 1){
        if(v & 1) res *= a;
    }
    return res;
}
 
void init(){
    fac[0] = 1;
    for(int i = 1; i <= N - 10; i++){
        fac[i] = fac[i-1] * i;
    }
    inv[N - 10] = qpow(fac[N - 10]);
 
    for(int i = N - 11; i >= 0; i--) {
        inv[i] = inv[i+1] * (i+1);
    }
 
    pow2[0] = 1;
 
    for(int i = 1; i <= N - 10; i++) {
        pow2[i] = 2 * pow2[i-1];
    }

}

Z C(int x, int y){
    if(x < y || x < 0 || y < 0) return 0;
    return fac[x] * inv[y] * inv[x - y];
}
```

### 矩阵运算类
```
/*
重载运算: 矩阵加 / 减 / 乘 / 乘方 / 比较
*/
#include <bits/stdc++.h>
#define ll long long
using namespace std;
constexpr int sz = 500;
constexpr int mod = 998244353;

struct matrix{
    int a[sz][sz];

    matrix() { memset(a, 0, sizeof a); }

    matrix(int b[sz][sz]) {
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++)
                a[i][j] = b[i][j];
    }
    matrix operator - (const matrix& T) const {
        matrix res;
        for (int i = 0; i < sz; i++)
            for (int j = 0; j < sz; j++) {
                res.a[i][j] = (a[i][j] - T.a[i][j] + mod) % mod;
            }
        return res;
    }

    matrix operator + (const matrix& T) const {
        matrix res;
        for (int i = 0; i < sz; i++)
            for (int j = 0; j < sz; j++) {
                res.a[i][j] = (a[i][j] + T.a[i][j]) % mod;
            }
        return res;
    }

    matrix operator * (const matrix& T) const {
        matrix res;
        int r;
        for(int i = 0; i < sz; i++)
            for(int k = 0; k < sz; k++) {
                r = a[i][k];
                for(int j = 0; j < sz; j++)
                    res.a[i][j] = (res.a[i][j] + 1ll * r * T.a[k][j]) % mod;
            }
        return res;
    }

    matrix operator ^ (ll x) const {
        matrix res, bas;
        for (int i = 0; i < sz; i++) res.a[i][i] = 1;
        for (int i = 0; i < sz; i++)
            for (int j = 0; j < sz; j++) bas.a[i][j] = a[i][j] % mod;
        while(x) {
            if (x & 1) res = res * bas;
            bas = bas * bas;
            x >>= 1;
        }
        return res;
    }

    bool operator == (const matrix& T) const {
        bool same = true;
        for(int i = 0; i < sz; i++)
            for(int j = 0; j < sz; j++)
                if(a[i][j] != T.a[i][j])
                    same = false;
        return same;
    }

};
```

### 莫队

```
void add(int x) {
    cnt[x]++;
    if (cnt[x] == 1) cur++;  // 示例：不同数个数
}

void remove(int x) {
    cnt[x]--;
    if (cnt[x] == 0) cur--;
}

while (L > l) add(a[--L]);
while (R < r) add(a[++R]);
while (L < l) remove_(a[L++]);
while (R > r) remove_(a[R--]);

```

### 字符串哈希 
```
using ull = unsigned long long;
ull base = 131;
ull mod1 = 1e9 + 21, mod2 = 1e9 + 7;

ull get_hash1(std::string s) {
    int len = s.size();
    ull ans = 0;
    for(int i = 0; i < len; i++) ans = (ans * base + (ull)s[i]) % mod1;
    return ans;
}

ull get_hash2(std::string s) {
    int len = s.size();
    ull ans = 0;
    for(int i = 0; i < len; i++) ans = (ans * base + (ull)s[i]) % mod2;
    return ans;
}

bool cmp(const std::string s, const std::string t){ // 相同返回 0, 不同返回 1
    bool f1 = get_hash1(s) != get_hash1(t);
    bool f2 = get_hash2(s) != get_hash2(t);
    return f1 or f2;
}

std::map<std::pair<ull, ull>, int> hash_map;//用 map 存双哈希值

```

### 后缀数组

#### 后缀数组（Suffix Array）能解决的问题

后缀数组（SA）是一种将字符串所有后缀按字典序排序的数据结构，通常配合高度数组（LCP）使用，可以高效解决大量字符串问题。

#### 一、基础能力

##### 1. 子串排序 / 字典序比较

- 任意两个子串 $s[l_1 \dots r_1]$ 和 $s[l_2 \dots r_2]$ 的字典序比较
- 转化为比较对应后缀的排名 + LCP

##### 2. 快速判断子串是否出现

- 判断字符串 $t$ 是否为 $s$ 的子串
- 在 SA 上二分查找

时间复杂度：
$$
O(|t| \log n)
$$

##### 3. 统计子串出现次数

- 利用 SA + 二分找到区间
- 区间长度即出现次数

#### 二、LCP（最长公共前缀）相关问题

##### 4. 任意两个后缀的 LCP

- 利用 height 数组 + RMQ：
$$
\text{LCP}(i, j) = \min(\text{height}[l+1 \dots r])
$$

##### 5. 最长重复子串

- 即：
$$
\max(\text{height}[i])
$$

##### 6. 最长公共子串（两个字符串）

- 将两个串拼接：
$$
s + \# + t
$$
- 找相邻后缀分别来自两个串的最大 LCP
- 
#### 三、子串统计问题

##### 8. 不同子串个数

公式：
$$
\frac{n(n+1)}{2} - \sum_{i=2}^{n} \text{height}[i]
$$

##### 9. 第 $k$ 小子串

- 每个后缀贡献：
$$
(n - \text{sa}[i] + 1) - \text{height}[i]
$$
- 按顺序枚举定位第 $k$ 个

##### 10. 子串去重统计

- 利用 height 去掉重复前缀

#### 四、重复结构与模式问题

##### 11. 最长重复且不重叠子串

- 二分答案 + 检查 height 区间是否合法

##### 12. 出现至少 $k$ 次的最长子串

- 二分长度 + 检查是否有连续 height $\ge L$

##### 13. 字符串周期问题

- 利用 LCP 判断周期性

#### 五、区间与排名问题

##### 14. 子串排名

- 求某个子串在所有子串中的排名

##### 15. 区间不同子串数

- 需要结合数据结构（如线段树 / 离线处理）


```
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstdio>

using namespace std;

constexpr int N = 1e6 + 10;

char s[N];
int n;
int m, p, rk[N << 1], oldrk[N << 1], sa[N << 1], id[N], cnt[N], height[N];

void SA() {
        
    scanf("%s", s + 1);
    n = strlen(s + 1);
    m = 128;

    for(int i = 1; i <= n; i++) cnt[rk[i] = s[i]]++;//记录每个字符数量, 以及每个位置的初始排名
    for(int i = 1; i <= m; i++) cnt[i] += cnt[i-1];//记录每个字符的排名区间
    for(int i = n; i >= 1; i--) sa[cnt[rk[i]]--] = i;//记录初始每个排名的位置编号

    for(int z = 1;; z <<= 1, m = p) {
        int cur = 0;
        for(int i = n - z + 1; i <= n; i++) id[++cur] = i;
        for(int i = 1; i <= n; i++)
            if(sa[i] > z) id[++cur] = sa[i] - z;
        
        fill(cnt, cnt + m + 1, 0);
        for(int i = 1; i <= n; i++) cnt[rk[i]]++;
        for(int i = 1; i <= m; i++) cnt[i] += cnt[i-1];
        for (int i = n; i >= 1; i--) sa[cnt[rk[id[i]]]--] = id[i];

        p = 0;
        copy(rk, rk + n + 1, oldrk);

        for(int i = 1; i <= n; i++) {
            if(oldrk[sa[i]] == oldrk[sa[i-1]] && 
                oldrk[sa[i] + z] == oldrk[sa[i-1] + z]) rk[sa[i]] = p;
            else 
                rk[sa[i]] = ++p;
        }

        if(p == n) break;
    }

    for(int i = 1, k = 0; i <= n; i++) {
        if(rk[i] == 0) continue;
        if(k) --k;
        while(s[i + k] == s[sa[rk[i] - 1] + k]) ++k;
        height[rk[i]] = k;
    }
}

```

### kmp 前缀函数 (前缀的最长 border)
```
vector<int> prefix_function(string s) {
    //0-index
    int n = (int)s.length();
    vector<int> pi(n);
    for (int i = 1; i < n; i++) {
        int j = pi[i - 1];
        while (j > 0 && s[i] != s[j]) j = pi[j - 1];
        if (s[i] == s[j]) j++;
        pi[i] = j;
    }
    return pi;
}
```

### xor 线性基

### 异或线性基（XOR Basis）

##### 数学形式

考虑一个集合 $S = \{a_1, a_2, \dots, a_n\}$，其中每个元素为非负整数。

在二进制表示下，将每个数视为 $\mathbb{F}_2$ 上的向量（即每一位是 $0/1$，运算为按位异或）。

异或线性基是一组向量 $\{b_1, b_2, \dots, b_k\}$，满足：

1. **线性无关性**：
   $$
   \forall i,\ b_i \ne \bigoplus_{j \ne i} c_j b_j,\quad c_j \in \{0,1\}
   $$

2. **张成性（生成原集合的线性空间）**：
   $$
   \forall x \in \text{span}(S),\quad x = \bigoplus_{i=1}^k c_i b_i,\quad c_i \in \{0,1\}
   $$

其中 $\oplus$ 表示按位异或运算。

---

#### 维度与性质

- 线性基大小 $k \le \lfloor \log_2(\max a_i) \rfloor + 1$
- 可表示的不同数的个数为：
  $$
  2^k
  $$
- 表示方式唯一（在固定基的情况下）

---

#### 能解决的问题类型

##### 1. 最大 / 最小异或值

- 求某个集合中选取若干数异或后的最大值：
  $$
  \max \bigoplus_{i \in T} a_i
  $$

- 或最小非零值

---

##### 2. 判断可表示性

- 判断某个数 $x$ 是否可以由集合 $S$ 异或得到：
  $$
  x \in \text{span}(S)\ ?
  $$

---

##### 3. 统计不同异或结果个数

- 不同异或结果数量为：
  $$
  2^k
  $$

---

##### 4. 第 $k$ 小异或值

- 将线性基转为“标准形”（如高位优先的行最简形式）
- 按二进制枚举组合求第 $k$ 小

---

##### 5. 子集异或相关计数问题

- 如统计满足某种条件的子集异或值个数
- 常结合高斯消元与状态压缩

---

##### 6. 图论问题

- 在图上维护路径异或（如树上路径 xor、带环图）
- 求：
  - 两点路径 xor 最大值
  - 环空间的 xor 基

---

##### 7. 区间 / 前缀问题

- 前缀异或 + 线性基：
  - 区间最大 xor
  - 子数组 xor 相关问题

---

#### 总结

异或线性基本质是：

> 在 $\mathbb{F}_2$ 上对整数进行高斯消元，维护一个极大线性无关组，从而高效处理“异或可达性”和“最优值”问题。

```
long long p[N];

void insert(long long x) {
    // 0 <= x < 2^60
    for (int i = 60; ~i; i--) {
        if (!(x >> i)) continue;
        if (!p[i]) {
            p[i] = x;
            break;
        }
        x ^= p[i];
    }
}
```

### mobius 函数 (反演)

#### 莫比乌斯函数（Möbius Function）

莫比乌斯函数 $\mu(n)$ 定义为：

$$
\mu(n) =
\begin{cases}
1, & n = 1 \\
(-1)^k, & n \text{ 是 } k \text{ 个不同质数的乘积} \\
0, & n \text{ 含有平方因子}
\end{cases}
$$

等价地：
- 若 $n = p_1 p_2 \cdots p_k$（所有 $p_i$ 互不相同），则 $\mu(n) = (-1)^k$
- 若存在 $p^2 \mid n$，则 $\mu(n) = 0$

---

#### 莫比乌斯反演（Möbius Inversion）

##### 形式一（约数形式）

若对任意正整数 $n$，有：

$$
f(n) = \sum_{d \mid n} g(d)
$$

则有反演公式：

$$
g(n) = \sum_{d \mid n} \mu(d)\, f\left(\frac{n}{d}\right)
$$

---

##### 形式二（对称形式）

若：

$$
f(n) = \sum_{d \mid n} g\left(\frac{n}{d}\right)
$$

则：

$$
g(n) = \sum_{d \mid n} \mu(d)\, f\left(\frac{n}{d}\right)
$$

---

##### 常见等价写法

也可以写成：

$$
g(n) = \sum_{d \mid n} \mu\left(\frac{n}{d}\right)\, f(d)
$$

```
vector<int> mu(N), primes;
vector<bool> isprime(N, true);

void mobius(int n) { // O(n) 筛 mobius 函数
    mu[1] = 1;

    for (int i = 2; i <= n; i++) {
        if (isprime[i]) {
            primes.push_back(i);
            mu[i] = -1;
        }
        for (int p : primes) {
            if (i * p > n) break;
            isprime[i * p] = false;

            if (i % p == 0) {
                mu[i * p] = 0;
                break;
            } else {
                mu[i * p] = -mu[i];
            }
        }
    }
}

int getmu(int x) { // O(sqrt(n)) 求单个 mobius 函数
    int s = 1;
    for(int i = 2; i * i <= x; i++) {
        int cnt = 0;
        if(x % i != 0) continue;
        while(x % i == 0) x /= i, cnt++;
        if(cnt > 1) return 0;
        s *= -1;
    }
    if(x > 1) s *= -1;
    return s;
}
```

### exgcd 求 线性同余方程
**扩展欧几里得算法**

时间复杂度：$O(\log n)$  
用于求解 $ax + by = \gcd(a, b)$ 的一组特解 $(x, y)$

模拟辗转相除的过程：

$$
ax_1 + by_1 = \gcd(a, b)
$$

$$
bx_2 + (a \bmod b)\, y_2 = \gcd(b, a \bmod b)
$$

由于：
$$
\gcd(a, b) = \gcd(b, a \bmod b)
$$

可得：
$$
ax_1 + by_1 
= bx_2 + (a \bmod b)\, y_2
= bx_2 + \bigl(a - \lfloor \tfrac{a}{b} \rfloor b\bigr) y_2
$$

整理得到：
$$
x_1 = y_2
$$

$$
y_1 = x_2 - \lfloor \tfrac{a}{b} \rfloor y_2
$$

```

void Exgcd(int a, int b, int &x, int &y)
{
    if(b == 0) {
        x = 1, y = 0;
        return;
    }
    Exgcd(b, a % b, x, y);
    tie(x, y) = tuple{y, x - (a / b) * y};
}

```

### 中国剩余定理

上面算法计算所得的 $x$ 对于任意 $i = 1,2,\cdots,k$ 满足
$$
x \equiv a_i \pmod{n_i}.
$$

当 $i \ne j$ 时，有
$$
m_j \equiv 0 \pmod{n_i},
$$
故
$$
c_j \equiv m_j \equiv 0 \pmod{n_i}.
$$

又有
$$
c_i \equiv m_i \cdot (m_i^{-1} \bmod n_i) \equiv 1 \pmod{n_i}.
$$

所以我们有：
$$
\begin{aligned}
x 
&\equiv \sum_{j=1}^{k} a_j c_j \pmod{n_i} \\
&\equiv a_i c_i \pmod{n_i} \\
&\equiv a_i \cdot m_i \cdot (m_i^{-1} \bmod n_i) \pmod{n_i} \\
&\equiv a_i \pmod{n_i}.
\end{aligned}
$$

```
LL CRT(int k, LL* a, LL* r) {
    LL n = 1, ans = 0;
    for (int i = 1; i <= k; i++) n = n * r[i];
    for (int i = 1; i <= k; i++) {
        LL m = n / r[i], b, y;
        exgcd(m, r[i], b, y);  // b * m mod r[i] = 1
        ans = (ans + a[i] * m * b % n) % n;
    }
    return (ans % n + n) % n;
}

```

### 三模数 ntt

$998244353, 1004535809, 469762049$

取 $3$ 为原根

### 多项式乘法 (ntt)

```

constexpr int N = 4e6 + 10;
constexpr int mod = 998244353; // 1004535809, 469762049

int n, m;
int F[N], G[N], rev[N];

int Pow(int a, int k){
    int res = 1;
    for(; k; a = 1ll * a * a % mod, k >>= 1){
        if(k & 1) res = 1ll * res * a % mod;
    }
    return res;
}

void NTT(int len, int *F, bool type){
    for(int i = 0; i < len; i++) if(i < rev[i]) swap(F[i], F[rev[i]]);
    for(int k = 1; k < len; k <<= 1){
        int eps = Pow(type ? 3 : 332748118, (mod - 1) / (k << 1));
        for(int i = 0; i < len; i += (k << 1)) {
            for(int j = i, g = 1; j < i + k; j++, g = 1ll * g * eps % mod) {
                int tmp1 = F[j], tmp2 = 1ll * g * F[j + k] % mod;
                F[j] = tmp1 + tmp2 >= mod ? tmp1 + tmp2 - mod : tmp1 + tmp2;
                F[j + k] = tmp1 - tmp2 < 0 ? tmp1 - tmp2 + mod : tmp1 - tmp2;
            }
        }
    }
}

void polymul(int *F, int *G){
    int len, ln;
    for(len = 1, ln = 0; len <= 2 * (n-1); len <<= 1, ln++);
    for(int i = 0; i < len; i++) rev[i] = (rev[i >> 1] >> 1) + (i & 1) * (1 << ln - 1);

    NTT(len, F, true), NTT(len, G, true);
    rep(i, 0, len - 1) F[i] = 1ll * F[i] * G[i] % mod;
    NTT(len, F, false);

    int Inv = Pow(len, mod - 2);
    rep(i, 0, len - 1) F[i] = 1ll * F[i] * Inv % mod;
}

```

### 三维点类

```
const double EPS = 1e-12;
const double pi = 3.141593;

struct Vec {
    double x, y, z;
    double len;
    
    Vec() {}
    Vec(double x, double y, double z): x(x), y(y), z(z) {}
    
    Vec operator + (const Vec& b) const {
        return Vec(x+b.x, y+b.y, z+b.z);
    }
    Vec operator - (const Vec& b) const {
        return Vec(x-b.x, y-b.y, z-b.z);
    }
    Vec operator * (double k) const {
        return Vec(x*k, y*k, z*k);
    }
};

double dot(Vec a, Vec b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec cross(Vec a, Vec b) {
    return Vec(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

double norm(Vec a) {
    return sqrt(dot(a,a));
}
```

### 最大流

```
inline bool bfs() {
    fep(i, 1, n) cur[i] = head[i], d[i] = inf;
    queue<int> q;
    q.push(S);
    d[S] = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        For(i, u) if (e[i].c && d[u] + 1 < d[e[i].v]) {
            d[e[i].v] = d[u] + 1;
            q.push(e[i].v);
        }
    }
    return d[T] ^ inf;
}

inline int dfs(int u, int flow) {
    if (u == T || !flow) return flow;
    int Sum = 0;
    for (int &i = cur[u]; i; i = e[i].nxt) {
        if (e[i].c && d[u] + 1 == d[e[i].v]) {
            int f = dfs(e[i].v, min(flow - Sum, e[i].c));
            e[i].c -= f;
            e[i ^ 1].c += f;
            Sum += f;
            if (Sum == flow) return Sum;
        }
    }
    return Sum;
}

inline int dinic() {
    int maxflow = 0;
    while (bfs()) maxflow += dfs(S, inf);
    return maxflow;
}

```

### zkw 费用流

```

inline bool spfa() {
    for (re int i = 1; i <= n; ++i) vis[i] = false, deep[i] = INF;
    deque<int> q;
    q.push_back(t);
    vis[t] = true;
    deep[t] = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        vis[u] = false;
        for (re int i = head[u]; i; i = e[i].nxt) {
            if (e[i ^ 1].c && deep[u] + e[i ^ 1].w < deep[e[i].v])
            {
                deep[e[i].v] = deep[u] + e[i ^ 1].w;
                if (!vis[e[i].v]) {
                    vis[e[i].v] = true;
                    if (!q.empty() && deep[e[i].v] < deep[q.front()])
                        q.push_front(e[i].v);
                    else
                        q.push_back(e[i].v);
                }
            }
        }
    }
    return deep[s] != INF;
}

inline int dfs(int u, int flow) {
    int sum = 0;
    vis[u] = true;
    if (u == t || !flow) return flow;
    for (re int i = head[u]; i; i = e[i].nxt) {
        if (!vis[e[i].v] && e[i].c && deep[u] - e[i].w == deep[e[i].v]) {
            int res = dfs(e[i].v, min(flow, e[i].c));
            e[i].c -= res;
            e[i ^ 1].c += res;
            sum += res;
            flow -= res;
            mincost += res * e[i].w;
            if (!flow) break;
        }
    }
    return sum;
}

inline int zkw() {
    int flow = 0;
    while (spfa()) {
        vis[t] = true;
        while (vis[t]) {
            for (re int i = 1; i <= n; ++i) vis[i] = false;
            flow += dfs(s, INF);
        }
    }
    return flow;
}
```

### 计算几何 (dls version)

```
#include <bits/stdc++.h>

#define rep(i, a, b) for(int i = a; i <= b; i++)
#define db double

using namespace std;

const db EPS = 1e-9;

inline int sign(db a) { return a < -EPS ? -1 : a > EPS; }

inline int cmp(db a, db b) { return sign(a - b); }

struct P {
    db x, y;

    P() {}
    P(db _x, db _y) : x(_x), y(_y) {}
    
    P operator + (P p) { return {x + p.x, y + p.y}; }
    P operator - (P p) { return {x - p.x, y - p.y}; }
    P operator * (db d) { return {x * d, y * d}; }
    P operator / (db d) { return {x / d, y / d}; }

    bool operator < (P p) const {
        int c = cmp(x, p.x);
        if(c) return c == -1;
        return cmp(y, p.y) == -1;
    }

    /*极角排序
    bool operator < (P p) const {
        if(quad() != b.quad()) return quad() < p.quad();
        return sign(det(p)) == -1;
    } 
    */

    bool operator == (P o) const {
        return cmp(x, o.x) == 0 && cmp(y, o.y) == 0;
    }

    db distTo(P p) { return (*this - p).abs(); }
    db alpha() { return atan2(y, x); }
    void read() { cin >> x >> y; }
    void write() { cout << "(" << x << "," << y << ")" << "\n"; }
    db abs() { return sqrt(abs2()); }
    db abs2() { return x * x + y * y; }
    P rot90() { return P(-y, x); }
    P unit() { return *this / abs(); }
    int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }

    db dot(P p) { return x * p.x + y * p.y; }
    db det(P p) { return x * p.y - y * p.x; }

    P rot(db an) { return {x * cos(an) - y * sin(an), x * sin(an) + y * cos(an)}; }

};

db cross(P p1, P p2, P p3) { return ((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)); }
db crossOp(P p1, P p2, P p3) { return sign(cross(p1, p2, p3)); }

bool chkLL(P p1, P p2, P q1, P q2) {//直线平行
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return sign(a1 + a2) != 0;
}

P isLL(P p1, P p2, P q1, P q2) { // 直线交点
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return (p1 * a2 + p2 * a1) / (a1 + a2);
}

bool intersect(db l1, db r1, db l2, db r2) { // 跨立实验
    if(l1 > r1) swap(l1, r1); if(l2 > r2) swap(l2, r2);
    return !( cmp(r1, l2) == -1 || cmp(r2, l1) == -1 );
}

bool isSS(P p1, P p2, P q1, P q2) { // 线段相交
    return intersect(p1.x, p2.x, q1.x, q2.x) && intersect(p1.y, p2.y, q1.y, q2.y) &&
    crossOp(p1, p2, q1) * crossOp(p1, p2, q2) <= 0 && crossOp(q1, q2, p1)
            * crossOp(q1, q2, p2) <= 0;
}

bool isSS_strict(P p1, P p2, P q1, P q2) { // 是否严格相交
    return crossOp(p1, p2, q1) * crossOp(p1, p2, q2) < 0 && crossOp(q1, q2, p1)
            * crossOp(q1, q2, p2) < 0;
}

bool isMiddle(db a, db m, db b) {
    return sign(a - m) == 0 || sign(b - m) == 0 || (a < m != b < m);
}

bool isMiddle(P a, P m, P b) {
    return isMiddle(a.x, m.x, b.x) && isMiddle(a.y, m.y, b.y);
}

bool onSeg(P p1, P p2, P q) { // 点在线段上
    return crossOp(p1, p2, q) == 0 && isMiddle(p1, q, p2);
}

bool onSeg_strict(P p1, P p2, P q) {
    return crossOp(p1, p2, q) == 0 && sign((q - p1).dot(p1 - p2) * sign((q - p2).dot(p1 - p2)) < 0);
}

P proj(P p1, P p2, P q) { // 投影
    P dir = p2 - p1;
    return p1 + dir * (dir.dot(q - p1) / dir.abs2());
}

P reflect(P p1, P p2, P q) { // 对称点 
    return proj(p1, p2, q) * 2 - q;
}

db nearest(P p1, P p2, P q) { // 点到线段的最短距离
    P h = proj(p1, p2, q);
    if(isMiddle(p1, h, p2))
        return q.distTo(h);
    return min(p1.distTo(q), p2.distTo(q));
}

db disSS(P p1, P p2, P q1, P q2) { // 线段到线段的最短距离
    if(isSS(p1, p2, q1, q2)) return 0;
    return min(min(nearest(p1, p2, q1), nearest(p1, p2, q2)), min(nearest(q1, q2, p1), nearest(q1, q2, p2)));
}

db rad(P p1, P p2) { // 向量夹角
    return atan2l(p1.det(p2), p1.dot(p2));
}

db area(vector<P> ps) { // 多边形面积
    db ret = 0; rep(i, 0, ps.size() - 1) ret += ps[i].det(ps[(i + 1) % ps.size()]);

    return ret / 2;
}

int contain(vector<P> ps, P p) { // 2: inside, 1: on_seg, 0: outside
    int n = ps.size(), ret = 0;
    rep(i, 0, n-1) {
        P u = ps[i], v = ps[(i + 1) % n];
        if(onSeg(u, v, p)) return 1;
        if(cmp(u.y, v.y) <= 0) swap(u, v);
        if(cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0) continue;
        ret ^= crossOp(p, u, v) > 0;
    }
    return ret * 2;
}

vector<P> convexHull(vector<P> ps) {
    int n = ps.size(); if(n <= 1) return ps;
    sort(ps.begin(), ps.end());
    vector<P> qs(n * 2); int k = 0;
    for(int i = 0; i < n; qs[k++] = ps[i++])
        while(k > 1 && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0) k--;
    for(int i = n - 2, t = k; i >= 0; qs[k++] = ps[i--])
        while(k > t && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0) k--;
    qs.resize(k - 1);
    return qs;
}

db dist2(P a, P b) { return (a - b).abs2(); }

db convexDiameter(vector<P> ps) { // 旋转卡尺
    int n = ps.size();
    if (n <= 1) return 0;
    if (n == 2) return ps[0].distTo(ps[1]);
    
    db maxd2 = 0;
    // 旋转卡尺寻找最远点对
    for (int i = 0, j = 1; i < n; i++) {
        // 当 (i, i+1, j+1) 的面积大于 (i, i+1, j) 时，j 向后移
        // 这里使用 det (叉积) 计算平行四边形面积
        while (sign((ps[(i + 1) % n] - ps[i]).det(ps[(j + 1) % n] - ps[i]) - 
                   (ps[(i + 1) % n] - ps[i]).det(ps[j] - ps[i])) > 0) {
            j = (j + 1) % n;
        }
        maxd2 = max({maxd2, dist2(ps[i], ps[j]), dist2(ps[(i + 1) % n], ps[j])});
    }
    return sqrt(maxd2);
}

// 辅助函数：将凸包的顶点调整为从左下角开始，且为逆时针
void reorder(vector<P>& ps) {
    int pos = 0;
    for(int i = 1; i < ps.size(); i++) {
        if(ps[i].y < ps[pos].y || (ps[i].y == ps[pos].y && ps[i].x < ps[pos].x)) pos = i;
    }
    rotate(ps.begin(), ps.begin() + pos, ps.end());
}

vector<P> minkowskiSum(vector<P> A, vector<P> B) { // 闵可夫斯基和
    reorder(A); reorder(B);
    int n = A.size(), m = B.size();
    
    // 构造差分向量
    vector<P> v1(n), v2(m);
    for(int i = 0; i < n; i++) v1[i] = A[(i + 1) % n] - A[i];
    for(int i = 0; i < m; i++) v2[i] = B[(i + 1) % m] - B[i];
    
    vector<P> res;
    res.push_back(A[0] + B[0]);
    
    // 归并排序边向量（按极角）
    int i = 0, j = 0;
    while(i < n && j < m) {
        if(sign(v1[i].det(v2[j])) > 0) res.push_back(res.back() + v1[i++]);
        else if(sign(v1[i].det(v2[j])) < 0) res.push_back(res.back() + v2[j++]);
        else res.push_back(res.back() + v1[i++] + v2[j++]); // 平行向量合并
    }
    while(i < n) res.push_back(res.back() + v1[i++]);
    while(j < m) res.push_back(res.back() + v2[j++]);
    
    // 重新跑一遍凸包以去重或处理共线（可选，视具体题目精度要求而定）
    return convexHull(res); 
}
```
