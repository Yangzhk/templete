### 常用的 poly 板子
```
#include <bits/stdc++.h>
using namespace std;

const int mod = 998244353, G = 3;
const int P1 = 998244353, P2 = 1004535809, P3 = 469762049;

int pw(int x, long long k, int p = mod) {
    assert(k >= 0 && p > 1);
    int r = 1;
    x = (x % p + p) % p;
    for (; k; k >>= 1, x = 1LL * x * x % p)
        if (k & 1) r = 1LL * r * x % p;
    return r;
}
int inv_mod(int x, int p = mod) {
    assert(x % p != 0);
    return pw(x, p - 2, p); // p 为素数
}

// p 只能是 P1 / P2 / P3，其原根均为 3。
void ntt(vector<int> &a, int n, bool inv, int p = mod) {
    assert(n > 0 && (n & (n - 1)) == 0);
    assert((p == P1 && n <= (1 << 23)) ||
           (p == P2 && n <= (1 << 21)) ||
           (p == P3 && n <= (1 << 26)));
    a.resize(n, 0);
    static vector<int> rev;
    static int last = -1;
    if (last != n) {
        rev.assign(n, 0); last = n;
        for (int i = 1; i < n; i++)
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? n >> 1 : 0);
    }
    for (int i = 0; i < n; i++)
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    for (int m = 1; m < n; m <<= 1) {
        int wn = pw(G, (p - 1) / (m << 1), p);
        if (inv) wn = inv_mod(wn, p);
        for (int i = 0; i < n; i += m << 1) {
            int w = 1;
            for (int j = 0; j < m; j++, w = 1LL * w * wn % p) {
                int x = a[i + j], y = 1LL * w * a[i + j + m] % p;
                a[i + j] = x + y >= p ? x + y - p : x + y;
                a[i + j + m] = x - y < 0 ? x - y + p : x - y;
            }
        }
    }
    if (inv) {
        int iv = inv_mod(n, p);
        for (int &x : a) x = 1LL * x * iv % p;
    }
}

// NTT 卷积，p 只能是三种 NTT 质数之一
vector<int> mul_p(vector<int> A, vector<int> B, int p = mod) {
    if (A.empty() || B.empty()) return {};
    int len = (int)A.size() + (int)B.size() - 1, n = 1;
    while (n < len) n <<= 1;
    for (int &x : A) x = (x % p + p) % p;
    for (int &x : B) x = (x % p + p) % p;
    ntt(A, n, false, p); ntt(B, n, false, p);
    for (int i = 0; i < n; i++) A[i] = 1LL * A[i] * B[i] % p;
    ntt(A, n, true, p);
    A.resize(len);
    return A;
}

// 重构三模 CRT 结果，返回 [0, P1*P2*P3) 中的代表元
__int128 crt128(int r1, int r2, int r3) {
    static const int i12 = inv_mod(P1 % P2, P2);
    static const long long M12 = 1LL * P1 * P2;
    static const int i123 = inv_mod((int)(M12 % P3), P3);
    long long t1 = 1LL * ((r2 - r1) % P2 + P2) % P2 * i12 % P2;
    long long x = r1 + 1LL * P1 * t1;
    long long t2 = ((r3 - x % P3) % P3 + P3) % P3 * i123 % P3;
    return (__int128)x + (__int128)M12 * t2;
}

int crt_mod(int r1, int r2, int r3, int p) {
    assert(p > 0);
    return (int)(crt128(r1, r2, r3) % p);
}

// 有符号精确结果，要求 |真实系数| < P1*P2*P3/2 且能放进 long long
long long crt(int r1, int r2, int r3) {
    __int128 x = crt128(r1, r2, r3);
    const __int128 M = (__int128)P1 * P2 * P3;
    if (x > M / 2) x -= M;
    assert(x >= LLONG_MIN && x <= LLONG_MAX);
    return (long long)x;
}

// 任意目标模数卷积：系数先规范到 [0,p)，并要求普通整数卷积系数 < P1*P2*P3
vector<int> mul_any(const vector<int> &a, const vector<int> &b, int p) {
    assert(p > 0);
    if (a.empty() || b.empty()) return {};
    vector<int> A = a, B = b;
    for (int &x : A) x = ((long long)x % p + p) % p;
    for (int &x : B) x = ((long long)x % p + p) % p;
    auto c1 = mul_p(A, B, P1), c2 = mul_p(A, B, P2), c3 = mul_p(A, B, P3);
    vector<int> r(c1.size());
    for (int i = 0; i < (int)r.size(); i++) r[i] = crt_mod(c1[i], c2[i], c3[i], p);
    return r;
}

vector<long long> mul_exact(const vector<int> &a, const vector<int> &b) {
    if (a.empty() || b.empty()) return {};
    auto c1 = mul_p(a, b, P1), c2 = mul_p(a, b, P2), c3 = mul_p(a, b, P3);
    vector<long long> r(c1.size());
    for (int i = 0; i < (int)r.size(); i++) r[i] = crt(c1[i], c2[i], c3[i]);
    return r;
}

// 普通二进制多项式快速幂：计算 A(x)^k mod x^n
vector<int> binpow(vector<int> a, long long k, int n) {
    assert(k >= 0 && n >= 0);
    if (n == 0) return {};
    if ((int)a.size() > n) a.resize(n);
    vector<int> r{1};
    while (k) {
        if (k & 1) {
            r = mul_p(r, a);
            if ((int)r.size() > n) r.resize(n);
        }
        k >>= 1;
        if (k) {
            a = mul_p(a, a);
            if ((int)a.size() > n) a.resize(n);
        }
    }
    r.resize(n, 0);
    return r;
}

// 原 Poly 的导数、积分、单点求值和多项式平移的 vector 版本
vector<int> poly_deriv(const vector<int> &a) {
    if (a.size() <= 1) return {0};
    vector<int> r(a.size() - 1);
    for (int i = 1; i < (int)a.size(); i++) r[i - 1] = 1LL * a[i] * i % mod;
    return r;
}
vector<int> poly_integ(const vector<int> &a) {
    int n = (int)a.size();
    if (n == 0) return {};
    assert(n < mod);
    vector<int> iv(n + 1), r(n + 1);
    iv[1] = 1;
    for (int i = 2; i <= n; i++) iv[i] = 1LL * (mod - mod / i) * iv[mod % i] % mod;
    for (int i = 0; i < n; i++) r[i + 1] = 1LL * a[i] * iv[i + 1] % mod;
    return r;
}
int poly_eval(const vector<int> &a, int x) {
    x = (x % mod + mod) % mod;
    int r = 0;
    for (int i = (int)a.size() - 1; i >= 0; i--)
        r = (1LL * r * x + a[i]) % mod;
    return r;
}
vector<int> poly_shift(const vector<int> &a, int c) { // f(x+c)
    int n = a.size();
    if (!n) return {};
    assert(n < mod);
    c = (c % mod + mod) % mod;
    vector<int> fac(n), ifac(n), A(n), C(n), r(n);
    fac[0] = 1;
    for (int i = 1; i < n; i++) fac[i] = 1LL * fac[i - 1] * i % mod;
    ifac[n - 1] = inv_mod(fac[n - 1]);
    for (int i = n - 2; i >= 0; i--) ifac[i] = 1LL * ifac[i + 1] * (i + 1) % mod;
    for (int i = 0; i < n; i++) A[n - 1 - i] = 1LL * a[i] * fac[i] % mod;
    for (int i = 0, ck = 1; i < n; i++, ck = 1LL * ck * c % mod)
        C[i] = 1LL * ck * ifac[i] % mod;
    auto h = mul_p(A, C);
    for (int j = 0; j < n; j++) r[j] = 1LL * h[n - 1 - j] * ifac[j] % mod;
    return r;
}

// 斯特林数整行 O(n log^2 n) / O(n log n)
vector<int> stirling1(int n) {
    assert(0 <= n && n < mod);
    if (n == 0) return {1};
    if (n == 1) return {0, 1};
    int m = n / 2;
    auto a = stirling1(m);
    auto b = poly_shift(a, m);
    auto c = mul_p(a, b);
    if (n & 1) c = mul_p(c, {n - 1, 1});
    c.resize(n + 1);
    return c;
}

// 第一类斯特林数整行（无符号），O(n log^2 n)
// 返回 a[k] = [n k]
vector<int> stirling1(int n) {
    auto dfs = [&](auto &&self, int l, int r) -> vector<int> {
        if (l == r) return {1};
        if (r - l == 1) return {l, 1};
        int mid = (l + r) / 2;
        auto a = self(self, l, mid);
        auto b = self(self, mid, r);
        return mul_p(a, b);
    };
    return dfs(dfs, 0, n);
}

vector<int> stirling2(int n) {
    assert(0 <= n && n < mod);
    vector<int> fac(n + 1), ifac(n + 1);
    fac[0] = 1;
    for (int i = 1; i <= n; i++) fac[i] = 1LL * fac[i - 1] * i % mod;
    ifac[n] = inv_mod(fac[n]);
    for (int i = n; i >= 1; i--) ifac[i - 1] = 1LL * ifac[i] * i % mod;
    vector<int> A(n + 1), B(n + 1);
    for (int i = 0; i <= n; i++) {
        A[i] = (i & 1) ? (ifac[i] ? mod - ifac[i] : 0) : ifac[i];
        B[i] = 1LL * pw(i, n) * ifac[i] % mod;
    }
    auto c = mul_p(A, B);
    c.resize(n + 1);
    return c;
}

/ 第一类无符号 / 第二类斯特林数整列
// first = true : 第一类无符号
// first = false: 第二类
// 返回 r[n] = S(n,k), 0 <= n <= N
vector<int> stirling_col(int N, int k, bool first) {
    assert(0 <= N && N < mod && k >= 0);
    vector<int> r(N + 1);
    if (k > N) return r;
    if (k == 0) {
        r[0] = 1;
        return r;
    }

    vector<int> fac(N + 1), ifac(N + 1);
    fac[0] = 1;
    for (int i = 1; i <= N; i++)
        fac[i] = 1LL * fac[i - 1] * i % mod;

    ifac[N] = inv_mod(fac[N]);
    for (int i = N; i >= 1; i--)
        ifac[i - 1] = 1LL * ifac[i] * i % mod;

    int m = N - k;
    vector<int> a(m + 1);

    for (int i = 0; i <= m; i++) {
        if (first)
            a[i] = 1LL * fac[i] * ifac[i + 1] % mod;
        else
            a[i] = ifac[i + 1];
    }

    auto b = binpow(a, k, m + 1);

    for (int n = k; n <= N; n++)
        r[n] = 1LL * fac[n] * ifac[k] % mod * b[n - k] % mod;

    return r;
}

vector<int> stirling1_col(int N, int k) {
    return stirling_col(N, k, true);
}

vector<int> stirling2_col(int N, int k) {
    return stirling_col(N, k, false);
}

// Bostan-Mori: a_n = sum_{i=1..k} c[i-1]*a_{n-i}
int linear_recur(vector<int> a, vector<int> c, long long n) {
    int k = c.size();
    assert(k > 0 && (int)a.size() == k && n >= 0);
    if (n < k) return a[(int)n];
    vector<int> Q(k + 1); Q[0] = 1;
    for (int i = 0; i < k; i++) Q[i + 1] = (mod - c[i]) % mod;
    vector<int> P = mul_p(a, Q);
    P.resize(k);
    while (n > 0) {
        vector<int> Qm = Q;
        for (int i = 1; i < (int)Qm.size(); i += 2) Qm[i] = (mod - Qm[i]) % mod;
        auto U = mul_p(P, Qm), V = mul_p(Q, Qm);
        for (int i = 0; i < (int)P.size(); i++) {
            int j = 2 * i + (int)(n & 1);
            P[i] = j < (int)U.size() ? U[j] : 0;
        }
        for (int i = 0; i < (int)Q.size(); i++) {
            int j = 2 * i;
            Q[i] = j < (int)V.size() ? V[j] : 0;
        }
        n >>= 1;
    }
    return 1LL * P[0] * inv_mod(Q[0]) % mod;
}

// 普通拉格朗日插值 O(n^2)，x[i] 两两不同（模 mod）
int lagrange_general(const vector<int> &x, const vector<int> &y, long long k) {
    int n = x.size();
    assert(n == (int)y.size() && n < mod);
    int t = (k % mod + mod) % mod;
    int ans = 0;
    for (int i = 0; i < n; i++) {
        int xi = (x[i] % mod + mod) % mod;
        if (xi == t) return (y[i] % mod + mod) % mod;
        long long a = 1, b = 1;
        for (int j = 0; j < n; j++) if (i != j) {
            int xj = (x[j] % mod + mod) % mod;
            a = a * ((t - xj + mod) % mod) % mod;
            b = b * ((xi - xj + mod) % mod) % mod;
        }
        assert(b != 0);
        ans = (ans + (long long)( (y[i] % mod + mod) % mod) * a % mod * inv_mod(b) % mod) % mod;
    }
    return ans;
}

// 连续点值拉格朗日插值 O(n)，y[i] = f(i+1)
int lagrange_continuous(const vector<int> &y, long long k) {
    int n = y.size();
    assert(n < mod);
    if (n == 0) return 0;
    int t = (k % mod + mod) % mod;
    if (1 <= t && t <= n) return (y[t - 1] % mod + mod) % mod;
    vector<int> fac(n), ifac(n), pre(n + 1), suf(n + 2);
    fac[0] = 1;
    for (int i = 1; i < n; i++) fac[i] = 1LL * fac[i - 1] * i % mod;
    ifac[n - 1] = inv_mod(fac[n - 1]);
    for (int i = n - 2; i >= 0; i--) ifac[i] = 1LL * ifac[i + 1] * (i + 1) % mod;
    pre[0] = suf[n + 1] = 1;
    for (int i = 1; i <= n; i++) pre[i] = 1LL * pre[i - 1] * ((t - i + mod) % mod) % mod;
    for (int i = n; i >= 1; i--) suf[i] = 1LL * suf[i + 1] * ((t - i + mod) % mod) % mod;
    int ans = 0;
    for (int i = 1; i <= n; i++) {
        int z = 1LL * ((y[i - 1] % mod + mod) % mod) * pre[i - 1] % mod * suf[i + 1] % mod
                    * ifac[i - 1] % mod * ifac[n - i] % mod;
        if ((n - i) & 1) z = (mod - z) % mod;
        ans = (ans + z) % mod;
    }
    return ans;
}
```


### 子集卷积 fwt

#### 1. 解决什么问题

给定两个长度均为 $N=2^k$ 的数组 $a,b$，求**位运算卷积**：

$$
c_s=\sum_{i\operatorname{\,op\,}j=s}a_i b_j,
\qquad \operatorname{op}\in\{\mathrm{OR},\mathrm{AND},\mathrm{XOR}\}.
$$

暴力枚举 $(i,j)$ 是 $O(N^2)$。FWT（广义地也包括 OR/AND 的子集 Zeta 变换）把卷积转成**变换域逐点乘法**：

$$
\boxed{\ c=T^{-1}\bigl(T(a)\cdot T(b)\bigr)\ }
$$

这里的 $\cdot$ 是逐项相乘，整体复杂度 $O(N\log N)=O(k2^k)$。

#### 2. 三种变换的含义

设 $S,T$ 是二进制掩码，$\subseteq$ 表示二进制位意义上的包含。

| 运算 | 正变换 $\widehat a_S$ | 单个二进制位的合并 $(x,y)\to(\cdots)$ | 逆变换 |
|---|---|---|---|
| **OR** | $\displaystyle\sum_{T\subseteq S}a_T$ | $(x,x+y)$ | $(x,y-x)$ |
| **AND** | $\displaystyle\sum_{T\supseteq S}a_T$ | $(x+y,y)$ | $(x-y,y)$ |
| **XOR** | $\displaystyle\sum_T(-1)^{\operatorname{popcount}(S\,\&\,T)}a_T$ | $(x+y,x-y)$ | 同样合并，最后整体除以 $N$ |

**为什么可以逐点乘？**

- **OR**：$i\mathbin{|}j\subseteq S$ 等价于 $i\subseteq S$ 且 $j\subseteq S$，所以 $\widehat c_S=\widehat a_S\widehat b_S$。
- **AND**：$i\mathbin{\&}j\supseteq S$ 等价于 $i\supseteq S$ 且 $j\supseteq S$，所以同理。
- **XOR**：字符 $\chi_S(x)=(-1)^{\operatorname{popcount}(S\,\&\,x)}$ 满足 $\chi_S(i\oplus j)=\chi_S(i)\chi_S(j)$，因此同理。

**直观理解：**每一位只有 $0,1$ 两种状态。三种变换都是对每一位执行一个 $2\times2$ 线性变换，然后依次处理全部 $k$ 位。XOR 的变换矩阵

$$
H=\begin{pmatrix}1&1\\1&-1\end{pmatrix},\qquad H^2=2I,
$$

所以每一位逆变换要除以 $2$，全部 $k$ 位合起来就是除以 $2^k=N$。

> 注：严格分类时 OR/AND 常称为 **子集 Zeta 变换及 Möbius 反演**；竞赛中习惯把三者合称 FWT。

#### 3. C++17 模板（模 $998244353$）

输入长度可不足 $2^k$，卷积函数会自动补零。`type`：`0=OR`、`1=AND`、`2=XOR`。

```cpp
#include <bits/stdc++.h>
using namespace std;
#define int long long
const int mod = 998244353;

int qpow(int a, int b) {
    int res = 1;
    for(; b; b >>= 1, a = a * a % mod)
        if(b & 1) res = res * a % mod;
    return res;
}

// OR: [x, y] -> [x, x + y]
void fwt_or(vector<int>& a, bool inv = false) {
    int n = a.size();
    for(int len = 1; len < n; len <<= 1)
        for(int i = 0; i < n; i += len << 1)
            for(int j = 0; j < len; j++) {
                int &x = a[i + j], &y = a[i + j + len];
                if(!inv) {
                    y += x;
                    if(y >= mod) y -= mod;
                } else {
                    y -= x;
                    if(y < 0) y += mod;
                }
            }
}

// AND: [x, y] -> [x + y, y]
void fwt_and(vector<int>& a, bool inv = false) {
    int n = a.size();
    for(int len = 1; len < n; len <<= 1)
        for(int i = 0; i < n; i += len << 1)
            for(int j = 0; j < len; j++) {
                int &x = a[i + j], &y = a[i + j + len];
                if(!inv) {
                    x += y;
                    if(x >= mod) x -= mod;
                } else {
                    x -= y;
                    if(x < 0) x += mod;
                }
            }
}

// XOR: [x, y] -> [x + y, x - y]
void fwt_xor(vector<int>& a, bool inv = false) {
    int n = a.size();
    for(int len = 1; len < n; len <<= 1)
        for(int i = 0; i < n; i += len << 1)
            for(int j = 0; j < len; j++) {
                int &x = a[i + j], &y = a[i + j + len];
                int u = x, v = y;
                x = (u + v) % mod;
                y = (u - v + mod) % mod;
            }
    if(inv) {
        int t = qpow(n, mod - 2);  // 除以 n
        for(auto &x : a) x = x * t % mod;
    }
}

// type: 0 = OR, 1 = AND, 2 = XOR
vector<int> bit_conv(vector<int> a, vector<int> b, int type) {
    int n = 1;
    while(n < max(a.size(), b.size())) n <<= 1;
    a.resize(n); b.resize(n);

    auto trans = [&](vector<int>& v, bool inv) {
        if(type == 0) fwt_or(v, inv);
        else if(type == 1) fwt_and(v, inv);
        else fwt_xor(v, inv);
    };

    trans(a, false);
    trans(b, false);
    for(int i = 0; i < n; i++) a[i] = a[i] * b[i] % mod;
    trans(a, true);
    return a;
}
```

调用：`auto c = bit_conv(a, b, 2);` 即计算 XOR 卷积。

**注意：**XOR 逆变换要求 $N$ 在模意义下有逆元；对模 $998244353$、通常的 $2^k$ 长度没有问题。三种卷积的结果索引均在 $[0,N-1]$ 内。

---

### 99% 用不到的 poly 板子
- **依赖**：原模板已有的 `mod`、`pw`、`inv_mod`、`mul_p`、`poly_deriv`、`poly_integ`。
- **调用**：新增函数直接定义在全局作用域，例如 `inv(a,n)`；为避免与 `std::norm` 重名，辅助函数使用 `norm_mod`。
- **长度**：FPS 系列都计算 `mod x^n`，返回**恰好 n 项**，尾项保留 0；涉及取模的辅助函数只在多点求值/插值内部使用。
- **限制**：涉及积分、阶乘的运算按原模数要求次数低于 `mod`；卷积还受原 NTT 最大长度限制。

#### 结论速览
| 函数 | 数学意义 / 前提 | 复杂度 |
|---|---|---|
| `inv(a,n)` | 求 $a(x)^{-1}\bmod x^n$，$a_0\ne0$ | $O(n\log n)$ |
| `ln(a,n)` | 求 $\ln a(x)\bmod x^n$，$a_0=1$ | $O(n\log n)$ |
| `exp(a,n)` | 求 $e^{a(x)}\bmod x^n$，$a_0=0$ | $O(n\log n)$ |
| `sqrt(a,n)` | 求 $b(x)^2\equiv a(x)\pmod{x^n}$；无解返回空数组 | $O(n\log n)$ |
| `pow(a,k,n)` | 求 $a(x)^k\bmod x^n$，$k\ge0$ | $O(n\log n)$ |
| `ProductTree(x).eval(f)` | 一次计算 $f(x_1),\ldots,f(x_n)$ | $O(n\log^2 n)$，假设 $\deg f<n$ |
| `ProductTree(x).interpolate(y)` | 已知 $n$ 个互异点，求次数 $<n$ 的插值多项式 | $O(n\log^2 n)$ |

#### 关键原理
##### 牛顿迭代（逆 / exp / sqrt）
多项式求逆已知 $b\equiv a^{-1}\pmod{x^m}$ 时，用

$$b\leftarrow b(2-ab)\pmod{x^{2m}}$$

一次把精度翻倍。

对数与指数利用

$$\ln a=\int\frac{a'}a\,dx,\qquad b\leftarrow b(1+a-\ln b)\pmod{x^{2m}}.$$

开方利用

$$b\leftarrow\frac12\left(b+\frac{a}{b}\right)\pmod{x^{2m}}.$$

若 $a(x)=x^{2t}h(x)$，先对 $h$ 开方再乘 $x^t$；首个非零系数的次数为奇数时无解。常数项开方使用 Tonelli–Shanks。

##### 多项式快速幂
设 $a=x^tc(1+u)$，其中 $u(0)=0$、$c\ne0$，则

$$a^k=x^{tk}c^k\exp\bigl(k\ln(1+u)\bigr).$$

如果 $tk\ge n$，答案直接为零。实现已处理 $a_0=0$ 和 $k=0$。

##### 多点求值与插值
以分治乘积树存储 $M(x)=\prod_i(x-x_i)$。求值时递归取模；插值时用

$$f(x)=\sum_i\frac{y_i}{M'(x_i)}\frac{M(x)}{x-x_i}.$$

通过同一棵乘积树做分治合并，避免朴素的 $O(n^2)$ 插值。求值过程必须用到多项式取余，因此仅保留一个**内部取余辅助函数**，不额外提供独立的多项式除法模板。

#### 完整补充代码
```cpp
using vi = vector<int>;

//取模规范化
int norm_mod(long long x) { x %= mod; if (x < 0) x += mod; return (int)x; }
//删除多项式末尾多余的零
void trim(vi &a) {
    for (int &x : a) x = norm_mod(x);
    while (!a.empty() && a.back() == 0) a.pop_back();
}
//截取前 n 项，不足补零
vi cut(const vi &a, int n) {
    vi b(n, 0);
    for (int i = 0; i < n && i < (int)a.size(); ++i) b[i] = norm_mod(a[i]);
    return b;
}
//多项式加法
vi add(const vi &a, const vi &b) {
    vi c(max(a.size(), b.size()), 0);
    for (int i = 0; i < (int)c.size(); ++i)
        c[i] = ((i < (int)a.size() ? (long long)a[i] : 0) +
                (i < (int)b.size() ? (long long)b[i] : 0)) % mod;
    return c;
}

// 1. FPS 求逆: a[0] != 0, 返回 1/a mod x^n, O(n log n)
vi inv(const vi &a, int n) {
    if (n == 0) return {};
    assert(!a.empty() && norm_mod(a[0]) != 0);
    vi b{inv_mod(norm_mod(a[0]))};
    for (int m = 2; (int)b.size() < n; m <<= 1) {
        int k = min(m, n);
        vi t = mul_p(cut(a, k), b);
        t.resize(k);
        for (int &x : t) x = x ? mod - x : 0;
        t[0] = (t[0] + 2) % mod;
        b = mul_p(b, t);
        b.resize(k);
    }
    return cut(b, n);
}

// 2. FPS ln: a[0] = 1, O(n log n)
vi ln(const vi &a, int n) {
    if (n == 0) return {};
    assert(!a.empty() && norm_mod(a[0]) == 1);
    if (n == 1) return {0};
    vi t = mul_p(poly_deriv(cut(a, n)), inv(a, n));
    t.resize(n - 1);
    return poly_integ(t);  // 结果长度 n
}

// 3. FPS exp: a[0] = 0, O(n log n)
vi exp(const vi &a, int n) {
    if (n == 0) return {};
    assert(a.empty() || norm_mod(a[0]) == 0);
    vi b{1};
    for (int m = 2; (int)b.size() < n; m <<= 1) {
        int k = min(m, n);
        vi t = cut(a, k), s = ln(b, k);
        for (int i = 0; i < k; ++i) t[i] = norm_mod((long long)t[i] - s[i]);
        t[0] = (t[0] + 1) % mod;
        b = mul_p(b, t);
        b.resize(k);
    }
    return cut(b, n);
}

// 模意义开方：返回较小的一个根；无根返回 -1（mod 为奇素数）
int mod_sqrt(int a) {
    a = norm_mod(a);
    if (a == 0) return 0;
    if (pw(a, (mod - 1) / 2) != 1) return -1;
    int q = mod - 1, s = 0;
    while (!(q & 1)) q >>= 1, ++s;
    int z = 2;
    while (pw(z, (mod - 1) / 2) != mod - 1) ++z;
    int c = pw(z, q), x = pw(a, (q + 1) / 2), t = pw(a, q), m = s;
    while (t != 1) {
        int i = 0, v = t;
        while (v != 1) { v = 1LL * v * v % mod; ++i; }
        assert(i < m);
        int b = pw(c, 1LL << (m - i - 1));
        x = 1LL * x * b % mod;
        c = 1LL * b * b % mod;
        t = 1LL * t * c % mod;
        m = i;
    }
    return min(x, mod - x);
}

// 4. FPS sqrt: 求 b^2 = a mod x^n；n>0 时无解返回空数组。
// 允许 a 前面有偶数个零；答案选定一个根。O(n log n)
vi sqrt(const vi &a, int n) {
    if (n == 0) return {};
    vi A = cut(a, n);
    int pos = 0;
    while (pos < n && A[pos] == 0) ++pos;
    if (pos == n) return vi(n, 0);
    if (pos & 1) return {};
    int sh = pos / 2, k = n - pos;
    vi f(A.begin() + pos, A.end());
    int c = mod_sqrt(f[0]);
    if (c == -1) return {};
    vi b{c};
    const int inv2 = (mod + 1) / 2;
    for (int m = 2; (int)b.size() < k; m <<= 1) {
        int len = min(m, k);
        vi t = mul_p(cut(f, len), inv(b, len));
        t.resize(len);
        vi nb(len);
        for (int i = 0; i < len; ++i)
            nb[i] = 1LL * ((i < (int)b.size() ? b[i] : 0) + (long long)t[i]) % mod * inv2 % mod;
        b.swap(nb);
    }
    vi ans(n, 0);
    for (int i = 0; i < (int)b.size(); ++i) ans[i + sh] = b[i];
    return ans;
}

// 5. FPS pow: 求 a^k mod x^n，k >= 0；允许 a[0] = 0。O(n log n)
vi pow(const vi &a, long long k, int n) {
    assert(k >= 0 && n >= 0);
    if (n == 0) return {};
    vi ans(n, 0);
    if (k == 0) { ans[0] = 1; return ans; }
    vi A = cut(a, n);
    int p = 0;
    while (p < n && A[p] == 0) ++p;
    if (p == n || p > (n - 1) / k) return ans;
    int sh = p * k, len = n - sh, c = A[p];
    vi b(len, 0);
    int ic = inv_mod(c);
    for (int i = 0; i < len && p + i < n; ++i) b[i] = 1LL * A[p + i] * ic % mod;
    vi t = ln(b, len);
    for (int &x : t) x = 1LL * x * (k % mod) % mod;
    vi u = exp(t, len);
    int ck = pw(c, k);
    for (int i = 0; i < len; ++i) ans[i + sh] = 1LL * u[i] * ck % mod;
    return ans;
}

// 内部辅助：对乘积树中的首一多项式 b 取余。
// 不是独立的多项式除法模板；供下面的 ProductTree 使用。
vi mod_monic(vi a, vi b) {
    trim(a); trim(b);
    assert(!b.empty() && b.back() == 1);
    if (a.size() < b.size()) return a;
    if (b.size() == 1) return {};
    int k = (int)a.size() - (int)b.size() + 1;
    if (min(k, (int)b.size()) <= 32) {
        // 小规模采用朴素消元，仅保留余数。
        for (int i = k - 1; i >= 0; --i) {
            int c = a[i + (int)b.size() - 1];
            for (int j = 0; j < (int)b.size(); ++j)
                a[i + j] = norm_mod((long long)a[i + j] - 1LL * c * b[j]);
        }
        a.resize(b.size() - 1);
        trim(a);
        return a;
    }
    vi ra = a, rb = b;
    reverse(ra.begin(), ra.end());
    reverse(rb.begin(), rb.end());
    vi q = mul_p(cut(ra, k), inv(rb, k));
    q.resize(k);
    reverse(q.begin(), q.end());
    vi qb = mul_p(q, b);
    vi r(b.size() - 1);
    for (int i = 0; i < (int)r.size(); ++i)
        r[i] = norm_mod((long long)a[i] - qb[i]);
    trim(r);
    return r;
}

// 6. 多点求值 & 快速插值: O(n log^2 n)；输入点可重复求值，插值需互不相同。
struct ProductTree {
    int n;
    vi x;
    vector<vi> tr;
    explicit ProductTree(const vi &points) : n(points.size()), x(points), tr(4 * max(1, n)) {
        for (int &v : x) v = norm_mod(v);
        if (n) build(1, 0, n);
    }
    void build(int u, int l, int r) {
        if (r - l == 1) { tr[u] = {norm_mod(-x[l]), 1}; return; }
        int m = (l + r) / 2;
        build(u * 2, l, m); build(u * 2 + 1, m, r);
        tr[u] = mul_p(tr[u * 2], tr[u * 2 + 1]);
    }
    void eval_dfs(int u, int l, int r, const vi &f, vi &ans) const {
        if (f.empty()) return;
        if (r - l == 1) { ans[l] = f[0]; return; }
        int m = (l + r) / 2;
        eval_dfs(u * 2, l, m, mod_monic(f, tr[u * 2]), ans);
        eval_dfs(u * 2 + 1, m, r, mod_monic(f, tr[u * 2 + 1]), ans);
    }
    vi eval(const vi &f) const {
        vi ans(n, 0);
        if (n) eval_dfs(1, 0, n, mod_monic(f, tr[1]), ans);
        return ans;
    }
    vi interp_dfs(int u, int l, int r, const vi &w) const {
        if (r - l == 1) return {w[l]};
        int m = (l + r) / 2;
        vi a = interp_dfs(u * 2, l, m, w);
        vi b = interp_dfs(u * 2 + 1, m, r, w);
        return add(mul_p(a, tr[u * 2 + 1]), mul_p(b, tr[u * 2]));
    }
    vi interpolate(const vi &y) const {
        assert((int)y.size() == n && n < mod);
        if (n == 0) return {};
        vi d = eval(poly_deriv(tr[1]));
        vi w(n);
        for (int i = 0; i < n; ++i) {
            assert(d[i] != 0); // 检查重根/重复插值点
            w[i] = 1LL * norm_mod(y[i]) * inv_mod(d[i]) % mod;
        }
        vi ans = interp_dfs(1, 0, n, w);
        ans.resize(n);
        return ans;
    }
};

```

#### 调用示例
```cpp
vector<int> a = {1, 2, 3};
auto ia = inv(a, 5);       // a^{-1} mod x^5
auto la = ln(a, 5);       // ln a mod x^5
auto ea = exp({0, 1}, 5); // exp(x) mod x^5
auto pa = pow(a, 10, 5);  // a^10 mod x^5

ProductTree tree({1, 2, 4});
auto ys = tree.eval({3, 2, 1});     // 求 f(1),f(2),f(4)
auto f = tree.interpolate(ys);      // 还原次数 < 3 的 f
```

#### 注意事项
1. `sqrt(a,n)`：当 `n>0` 时，**返回空数组表示无解**；系数全部为零时则返回长度 `n` 的零数组。只返回两个平方根中的一个（若存在）。
2. `ProductTree`：`eval` 可以有重复的输入点；`interpolate` 必须保证输入点两两不同（模 `mod`）。其复杂度以目标次数不超过点数为常见情形。
3. 你上传的原模板中，**`stirling1(int n)` 出现了两个同名实现**，合并成同一份 `.cpp` 前需二选一；另外注释行 `/ 第一类无符号 / 第二类斯特林数整列` 应改为 `// 第一类无符号 / 第二类斯特林数整列`，否则原文件本身不能编译。

**验证**：新增代码用随机测试核对 FPS 求逆、对数/指数、平方根、快速幂和多点求值/插值的结果。
