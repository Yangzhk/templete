- [模意义下运算类](#模意义下运算类)
- [矩阵运算类](#矩阵运算类)
- [莫队](#莫队)
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
    - [7. 多字符串最长公共子串](#7-多字符串最长公共子串)
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
  - [六、进阶应用](#六进阶应用)
    - [16. 后缀自动机 / 后缀树替代问题](#16-后缀自动机--后缀树替代问题)
    - [17. 最小表示法（循环串）](#17-最小表示法循环串)
    - [18. 多模式匹配](#18-多模式匹配)
  - [七、典型问题总结](#七典型问题总结)
  - [总结](#总结)
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
  - [总结](#总结-1)
- [mobius 函数 (反演)](#mobius-函数-反演)
  - [莫比乌斯函数（Möbius Function）](#莫比乌斯函数möbius-function)
  - [莫比乌斯反演（Möbius Inversion）](#莫比乌斯反演möbius-inversion)
    - [形式一（约数形式）](#形式一约数形式)
    - [形式二（对称形式）](#形式二对称形式)
    - [常见等价写法](#常见等价写法)
- [exgcd 求 线性同余方程](#exgcd-求-线性同余方程)
- [字符串哈希](#字符串哈希)
- [多项式乘法 (ntt)](#多项式乘法-ntt)
- [计算几何 2-D](#计算几何-2-d)
- [三维点类](#三维点类)
- [最大流](#最大流)
- [zkw 费用流](#zkw-费用流)

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

##### 7. 多字符串最长公共子串

- 滑动窗口 + 单调队列维护 LCP 最小值

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

#### 六、进阶应用

##### 16. 后缀自动机 / 后缀树替代问题

- SA + LCP 可以替代：
  - 后缀树（节省空间）
  - 后缀自动机（部分问题）

##### 17. 最小表示法（循环串）

- 倍增字符串后用 SA 求最小字典序表示

##### 18. 多模式匹配

- 多个模式串统一在 SA 上二分

#### 七、典型问题总结

后缀数组特别适用于以下类型：

- 子串相关问题（查找、统计、排名）
- 重复结构问题（LCP、重复子串）
- 多字符串比较问题
- 字典序问题

#### 总结

> 后缀数组的核心作用：
>
> **将“子串问题”转化为“后缀排序 + 区间 LCP 问题”，从而用排序和 RMQ 高效解决复杂字符串问题。**

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

### 多项式乘法 (ntt)
```
constexpr int N = 4e6 + 10;
constexpr int mod = 998244353;

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

### 计算几何 2-D

```
/*
2-D Point / Segment
*/

// 符号判断：0表示相等，1表示正，-1表示负
constexpr double eps = 1e-9;
int dcmp(double x) {
    if (fabs(x) < eps) return 0;
    return x < 0 ? -1 : 1;
}

struct Point { // 如果输入的是 int 类型, 就改成 int, 包括 cross, dot 也改
    double x, y;
    Point(double x = 0, double y = 0) : x(x), y(y) {}
};

typedef Point Vector;

Vector operator + (Vector A, Vector B) { return Vector(A.x + B.x, A.y + B.y); }
Vector operator - (Point A, Point B) { return Vector(A.x - B.x, A.y - B.y); }
Vector operator * (Vector A, double p) { return Vector(A.x * p, A.y * p); }
bool operator == (Point a, Point b) { return dcmp(a.x - b.x) == 0 && dcmp(a.y - b.y) == 0; }

// 叉积运算
double cross(Point a, Point b) {
    return a.x * b.y - a.y * b.x;
}

// 计算向量 ab 和 ac 的位置关系
double cross(Point a, Point b, Point c) {
    return cross(b - a, c - a);
}

// 基础点积：两个从原点出发的向量
double dot(Point a, Point b) {
    return a.x * b.x + a.y * b.y;
}

// 向量内角点积：计算向量 AB 和 AC 的点积
double dot(Point a, Point b, Point c) {
    return dot(b - a, c - a);
}

// 长度与距离
double Length(Vector A) { 
    return sqrt(dot(A, A)); 
}
double Distance(Point A, Point B) {
    return Length(A - B); 
}

// 向量夹角 (弧度)
double Angle(Point a, Point b) {
    return acos(dot(a, b) / Length(a) / Length(b));
}

double dist_to_line(Point p, Point a, Point b) {
    Point v1 = b - a, v2 = p - a;
    return abs(cross(v1, v2)) / Length(v1);
}

// 投影点：点 p 在直线 ab 上的投影点
Point GetLineProjection(Point p, Point a, Point b) { // ab 不能是零向量
    if (dcmp(dot(v, v)) == 0) return a; // 或直接返回 a
    Point v = b - a;
    return a + v * (dot(v, p - a) / dot(v, v));
}

double DistToSegment(Point p, Point a, Point b) {
    if (a == b) return Length(p - a);
    Point v1 = b - a, v2 = p - a, v3 = p - b;
    if (dot(v1, v2) < 0) return Length(v2); // 投影在 a 左侧
    if (dot(v1, v3) > 0) return Length(v3); // 投影在 b 右侧
    // 投影在线段内，使用叉积求面积/底边长
    return abs(cross(v1, v2)) / Length(v1);
}

bool OnSegment(Point p, Point a, Point b) {
    // 叉积为 0 保证共线，点积 <= 0 保证 p 在 a, b 之间
    return dcmp(cross(a - p, b - p)) == 0 &&
       dcmp(dot(a - p, b - p)) <= 0;
}

struct Line {
    Point p;   // 直线上的一点（起点）
    Vector v;  // 方向向量（从 p 指向直线另一个点的向量）
    double angle; // 极角，用于半平面交排序

    Line() {}
    // 通过两个点定义直线：p1 为起点，p2-p1 为方向向量
    Line(Point p1, Point p2) : p(p1), v(p2 - p1) {
        angle = atan2(v.y, v.x);
    }
};

// 辅助定义：线段（Segment）
// 在竞赛中，线段通常直接复用 Line 类，但在逻辑上只取 t 在 [0, 1] 之间的部分
typedef Line Segment;

bool isParallel(Point a, Point b, Point c, Point d) {
    return dcmp(cross(b - a, d - c)) == 0;
}

// 直线 P1P2 和 Q1Q2 的交点 (前提是已确认不平行)
Point get_line_intersection(Point p1, Point p2, Point q1, Point q2) {
    double a1 = cross(q2 - q1, p1 - q1);
    double a2 = cross(q2 - q1, p2 - q1);
    // 利用面积比例关系 (叉积即平行四边形面积)
    return { (p1.x * a2 - p2.x * a1) / (a2 - a1), (p1.y * a2 - p2.y * a1) / (a2 - a1) };
}

// 判断线段是否相交（跨立实验）
bool is_segment_intersection(Point a, Point b, Point c, Point d) {
    // 检查两条线段是否互相跨立
    return max(a.x, b.x) >= min(c.x, d.x) && max(c.x, d.x) >= min(a.x, b.x) &&
           max(a.y, b.y) >= min(c.y, d.y) && max(c.y, d.y) >= min(a.y, b.y) &&
           dcmp(cross(a, b, c)) * dcmp(cross(a, b, d)) <= 0 && 
           dcmp(cross(c, d, a)) * dcmp(cross(c, d, b)) <= 0;
}

/*
凸多边形
*/

vector<Point> get_convex_hull(vector<Point>& p) {
    int n = p.size();
    if (n <= 2) return p;
    sort(p.begin(), p.end(), [](Point a, Point b) {
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    vector<Point> hull;
    // 求下凸壳
    for (int i = 0; i < n; ++i) {
        while (hull.size() > 1 && cross(hull[hull.size() - 2], hull.back(), p[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(p[i]);
    }
    // 求上凸壳
    int lower_size = hull.size();
    for (int i = n - 2; i >= 0; --i) {
        while (hull.size() > lower_size && cross(hull[hull.size() - 2], hull.back(), p[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(p[i]);
    }
    hull.pop_back(); // 起点被重复添加了两次
    return hull;
}

double get_area(const vector<Point>& hull) {
    double area = 0;
    int n = hull.size();
    for (int i = 0; i < n; ++i) {
        area += cross(hull[i], hull[(i + 1) % n]);
    }
    return abs(area) / 2.0;
}

double ConvexDiameter(const vector<Point>& hull) { // 旋转卡尺
    int n = hull.size();
    if (n <= 1) return 0;
    if (n == 2) return Length(hull[0] - hull[1]);
    
    double max_d = 0;
    for (int i = 0, j = 1; i < n; i++) {
        // 旋转寻找距离边 i -> i+1 最远的顶点 j
        // 比较的是三角形面积（叉积）
        while (abs(cross(hull[i], hull[(i + 1) % n], hull[(j + 1) % n])) > 
               abs(cross(hull[i], hull[(i + 1) % n], hull[j]))) {
            j = (j + 1) % n;
        }
        max_d = max({max_d, Length(hull[i] - hull[j]), Length(hull[(i + 1) % n] - hull[j])});
    }
    return max_d;
}

// 0: 外部, 1: 内部, 2: 边界上
int IsPointInPolygon(Point p, const vector<Point>& poly) {
    bool wn = false;
    int n = poly.size();
    for (int i = 0; i < n; i++) {
        Point s1 = poly[i], s2 = poly[(i + 1) % n];

        // 1. 首先判断点是否在边上 (使用你板子里的 OnSegment)
        if (OnSegment(p, s1, s2)) return 2;

        // 2. 射线法核心：判断射线与边的相交
        // 这里的逻辑处理了水平边和端点重合的特殊情况
        if (((dcmp(s1.y - p.y) <= 0 && dcmp(p.y - s2.y) < 0) || 
             (dcmp(s2.y - p.y) <= 0 && dcmp(p.y - s1.y) < 0)) &&
            (dcmp(p.x - (s2.x - s1.x) * (p.y - s1.y) / (s2.y - s1.y) - s1.x) < 0)) {
            wn = !wn;
        }
    }
    return wn ? 1 : 0;
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
