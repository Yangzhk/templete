- [NTT 与基础框架](#ntt-与基础框架)
- [三模数 NTT（精确卷积 / 任意模数）](#三模数-ntt精确卷积--任意模数)
- [基础运算](#基础运算)
  - [加 / 减 / 乘](#加--减--乘)
  - [求导 / 积分](#求导--积分)
- [多项式求逆](#多项式求逆)
- [多项式 ln / exp](#多项式-ln--exp)
- [多项式开根](#多项式开根)
- [多项式快速幂](#多项式快速幂)
- [多项式除法 / 取模](#多项式除法--取模)
- [多项式平移 f(x+c)](#多项式平移-fxc)
- [多点求值](#多点求值)
- [快速插值](#快速插值)
- [常系数齐次线性递推 (Bostan–Mori)](#常系数齐次线性递推-bostanmori)
- [单点求值](#单点求值)
- [拉格朗日插值](#拉格朗日插值)
  - [一般插值](#一般插值)
  - [连续点值插值](#连续点值插值)
- [位运算卷积 FWT](#位运算卷积-fwt)
- [子集卷积](#子集卷积)
- [复杂度总结](#复杂度总结)


## NTT 与基础框架

模数 $998244353 = 119 \times 2^{23} + 1$，原根 $g = 3$，支持长度 $\le 2^{23}$ 的 NTT。

设计要点：

- `pw` / `inv_mod` 接受任意素数模 `p`，默认 `mod`，消掉原 `qpow` / `qpow_mod` 的重复。
- `ntt` 用 `bool inv` 表示是否逆变换（true=逆变换，与教科书一致），并以 `p` 为参数，三模数共用。
- `rev` 数组按当前长度缓存（`thread_local static`），避免每次卷积都重建。

```cpp
#include <bits/stdc++.h>
using namespace std;

class Poly {
public:
    static const int mod = 998244353;
    static const int G   = 3;

    vector<int> a;

    Poly() = default;
    Poly(vector<int> v) : a(move(v)) {}
    Poly(initializer_list<int> il) : a(il) {}

    int  size() const          { return (int)a.size(); }
    void resize(int n)         { a.resize(n, 0); }
    int& operator[](int i)     { return a[i]; }
    int  operator[](int i) const { return i < size() ? a[i] : 0; }

    // 模幂：默认模数 mod，三模数 NTT 时传入 p
    static int pw(int x, long long k, int p = mod) {
        int r = 1; x = ((x % p) + p) % p;
        for (; k; x = 1ll * x * x % p, k >>= 1)
            if (k & 1) r = 1ll * r * x % p;
        return r;
    }
    static int inv_mod(int x, int p = mod) { return pw(x, p - 2, p); }

private:
    static int next_pow2(int x) { int n = 1; while (n < x) n <<= 1; return n; }

    // 缓存 rev，长度变化时重建
    static const vector<int>& rev_buf(int n) {
        static thread_local vector<int> rev;
        static thread_local int last = -1;
        if (last != n) {
            rev.assign(n, 0);
            for (int i = 1; i < n; i++)
                rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? n >> 1 : 0);
            last = n;
        }
        return rev;
    }

public:
    // n 必须为 2 的幂；inv=true 即逆变换
    static void ntt(vector<int>& a, int n, bool inv, int p = mod, int g = G) {
        a.resize(n, 0);
        const auto& rev = rev_buf(n);
        for (int i = 0; i < n; i++)
            if (i < rev[i]) swap(a[i], a[rev[i]]);
        for (int m = 1; m < n; m <<= 1) {
            int wn = pw(inv ? inv_mod(g, p) : g, (p - 1) / (m << 1), p);
            for (int j = 0; j < n; j += m << 1) {
                int w = 1;
                for (int k = 0; k < m; k++, w = 1ll * w * wn % p) {
                    int x = a[j + k], y = 1ll * w * a[j + k + m] % p;
                    a[j + k]     = (x + y) % p;
                    a[j + k + m] = (x - y + p) % p;
                }
            }
        }
        if (inv) {
            int iv = inv_mod(n, p);
            for (int& x : a) x = 1ll * x * iv % p;
        }
    }

    // 单模数卷积：A * B (mod p)，p 必须 NTT 友好
    static vector<int> mul_p(vector<int> A, vector<int> B, int p) {
        if (A.empty() || B.empty()) return {};
        int len = (int)A.size() + (int)B.size() - 1, n = next_pow2(len);
        ntt(A, n, false, p);
        ntt(B, n, false, p);
        for (int i = 0; i < n; i++) A[i] = 1ll * A[i] * B[i] % p;
        ntt(A, n, true, p);
        A.resize(len);
        return A;
    }
```

---

## 三模数 NTT（精确卷积 / 任意模数）

三个 NTT 友好质数，原根均为 $3$：

$$P_1 = 998244353 = 119 \times 2^{23} + 1$$
$$P_2 = 1004535809 = 479 \times 2^{21} + 1$$
$$P_3 = 469762049 = 7 \times 2^{26} + 1$$

$P_1 P_2 P_3 \approx 4.7 \times 10^{26}$。`mul_exact` 假设结果 $< 2^{63}$ 直接返回 `long long`；任意模数下走 `mul_any`，CRT 合并时直接对目标模数取模，避免中间溢出。

```cpp
    // === 三模数 NTT ===
    static const int P1 = 998244353;   // 119 * 2^23 + 1
    static const int P2 = 1004535809;  // 479 * 2^21 + 1
    static const int P3 = 469762049;   // 7  * 2^26 + 1

    // 三模数 CRT → 直接对 p 取模（任意模数 NTT 用，避免大整数中转）
    static int crt3_mod(int r1, int r2, int r3, int p) {
        static const int  iv12  = inv_mod(P1 % P2, P2);
        static const long long M12 = 1ll * P1 * P2;
        static const int  iv123 = inv_mod((int)(M12 % P3), P3);
        long long t12  = 1ll * ((r2 - r1) % P2 + P2) % P2 * iv12 % P2;
        long long x12  = r1 + P1 * t12;                       // < P1*P2 < 2^60
        long long t123 = ((r3 - x12 % P3) % P3 + P3) % P3 * iv123 % P3;
        return (int)((x12 % p + (M12 % p) * t123 % p) % p);
    }

    // 三模数 CRT → 精确整数（要求结果 < 2^63）
    static long long crt3(int r1, int r2, int r3) {
        static const int  iv12  = inv_mod(P1 % P2, P2);
        static const long long M12 = 1ll * P1 * P2;
        static const int  iv123 = inv_mod((int)(M12 % P3), P3);
        long long t12  = 1ll * ((r2 - r1) % P2 + P2) % P2 * iv12 % P2;
        long long x12  = r1 + P1 * t12;
        long long t123 = ((r3 - x12 % P3) % P3 + P3) % P3 * iv123 % P3;
        return x12 + M12 * t123;
    }

    // 任意模数卷积（p 可非 NTT 友好）
    static vector<int> mul_any(const vector<int>& a, const vector<int>& b, int p) {
        auto c1 = mul_p(a, b, P1);
        auto c2 = mul_p(a, b, P2);
        auto c3 = mul_p(a, b, P3);
        vector<int> r(c1.size());
        for (int i = 0; i < (int)r.size(); i++)
            r[i] = crt3_mod(c1[i], c2[i], c3[i], p);
        return r;
    }

    // 精确卷积，结果系数需 < 2^63
    static vector<long long> mul_exact(const vector<int>& a, const vector<int>& b) {
        auto c1 = mul_p(a, b, P1);
        auto c2 = mul_p(a, b, P2);
        auto c3 = mul_p(a, b, P3);
        vector<long long> r(c1.size());
        for (int i = 0; i < (int)r.size(); i++)
            r[i] = crt3(c1[i], c2[i], c3[i]);
        return r;
    }
```

---

## 三模数 NTT（精确卷积 / 任意模数）

三个 NTT 友好质数，原根均为 $3$：

$$P_1 = 998244353 = 119 \times 2^{23} + 1$$
$$P_2 = 1004535809 = 479 \times 2^{21} + 1$$
$$P_3 = 469762049 = 7 \times 2^{26} + 1$$

$P_1 P_2 P_3 \approx 4.7 \times 10^{26}$，结果系数不超此值时精确。分别在三模数下 NTT 卷积，再 CRT 合并得精确整数结果。

```cpp
    // === 三模数 NTT（追加在 Poly 类内）===

    static const int P1 = 998244353;   // 119 * 2^23 + 1
    static const int P2 = 1004535809;  // 479 * 2^21 + 1
    static const int P3 = 469762049;   // 7 * 2^26 + 1
    static const int G3 = 3;           // 公共原根

    // 任意模数快速幂
    static int qpow_mod(int a, int k, int p) {
        int ans = 1;
        for (; k; a = 1ll * a * a % p, k >>= 1)
            if (k & 1) ans = 1ll * ans * a % p;
        return ans;
    }

    // 单模数 NTT（模数作为参数传入）
    static void NTT_mod(vector<int>& a, int limit, bool type, int p) {
        int bit = 0;
        while ((1 << bit) < limit) bit++;
        vector<int> rev(limit);
        for (int i = 0; i < limit; i++)
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
        for (int i = 0; i < limit; i++)
            if (i < rev[i]) swap(a[i], a[rev[i]]);

        for (int mid = 1; mid < limit; mid <<= 1) {
            int wn = qpow_mod(
                type ? G3 : qpow_mod(G3, p - 2, p),
                (p - 1) / (mid << 1), p);
            for (int j = 0; j < limit; j += (mid << 1)) {
                int w = 1;
                for (int k = 0; k < mid; k++, w = 1ll * w * wn % p) {
                    int x = a[j + k], y = 1ll * w * a[j + k + mid] % p;
                    a[j + k] = (x + y >= p) ? (x + y - p) : (x + y);
                    a[j + k + mid] = (x - y < 0) ? (x - y + p) : (x - y);
                }
            }
        }

        if (!type) {
            int inv_limit = qpow_mod(limit, p - 2, p);
            for (int i = 0; i < limit; i++)
                a[i] = 1ll * a[i] * inv_limit % p;
        }
    }

    // 单模数卷积
    static vector<int> conv_mod(const vector<int>& a, const vector<int>& b, int p) {
        if (a.empty() || b.empty()) return {};
        int n = (int)a.size(), m = (int)b.size();
        int limit = 1;
        while (limit < n + m - 1) limit <<= 1;

        vector<int> A = a, B = b;
        NTT_mod(A, limit, true, p);
        NTT_mod(B, limit, true, p);
        for (int i = 0; i < limit; i++)
            A[i] = 1ll * A[i] * B[i] % p;
        NTT_mod(A, limit, false, p);

        A.resize(n + m - 1);
        return A;
    }

    // CRT 合并三个余数 → 精确整数 (< P1*P2*P3)
    static long long crt_merge(int r1, int r2, int r3) {
        // 合并前两个模数 P1, P2
        static const int inv_p1_p2 = qpow_mod(P1 % P2, P2 - 2, P2);
        long long x12 = r1 + 1ll * P1 *
            (1ll * (r2 - r1 + P2) % P2 * inv_p1_p2 % P2);

        // 合并第三个模数 P3
        static const long long M12 = 1ll * P1 * P2;
        static const int inv_m12_p3 = qpow_mod(M12 % P3, P3 - 2, P3);
        __int128 x = (__int128)x12 +
            (__int128)M12 * ((r3 - x12 % P3 + P3) % P3 * inv_m12_p3 % P3);

        return (long long)(x % (__int128)(M12 * P3));
    }

    // 三模数精确卷积 → 返回 long long 系数
    static vector<long long> conv_exact(const vector<int>& a, const vector<int>& b) {
        auto c1 = conv_mod(a, b, P1);
        auto c2 = conv_mod(a, b, P2);
        auto c3 = conv_mod(a, b, P3);
        vector<long long> res(c1.size());
        for (size_t i = 0; i < res.size(); i++)
            res[i] = crt_merge(c1[i], c2[i], c3[i]);
        return res;
    }

    // 任意模数卷积（结果对 arbitrary_mod 取模，用于非 NTT 友好模数）
    static vector<int> conv_arbitrary(const vector<int>& a, const vector<int>& b,
                                      int mod) {
        auto exact = conv_exact(a, b);
        vector<int> res(exact.size());
        for (size_t i = 0; i < res.size(); i++)
            res[i] = exact[i] % mod;
        return res;
    }
```

---

## 基础运算

### 加 / 减 / 乘

```cpp
    Poly operator+(const Poly& rhs) const {
        Poly r; r.a.assign(max(size(), rhs.size()), 0);
        for (int i = 0; i < r.size(); i++)
            r.a[i] = ((*this)[i] + rhs[i]) % mod;
        return r;
    }
    Poly operator-(const Poly& rhs) const {
        Poly r; r.a.assign(max(size(), rhs.size()), 0);
        for (int i = 0; i < r.size(); i++)
            r.a[i] = ((*this)[i] - rhs[i] + mod) % mod;
        return r;
    }
    Poly operator*(const Poly& rhs) const {
        return Poly(mul_p(a, rhs.a, mod));
    }
    Poly& operator*=(const Poly& rhs) { return *this = *this * rhs; }
```

### 求导 / 积分

$$
\frac{d}{dx} \sum a_i x^i = \sum (i+1) a_{i+1} x^i,\qquad
\int \sum a_i x^i = \sum \frac{a_i}{i+1} x^{i+1}
$$

积分的逆元用线性递推 $i^{-1} = -\lfloor p/i \rfloor \cdot (p \bmod i)^{-1}$ 而非快速幂，常数更小。

```cpp
    Poly deriv() const {
        if (size() <= 1) return Poly({0});
        vector<int> r(size() - 1);
        for (int i = 1; i < size(); i++) r[i - 1] = 1ll * a[i] * i % mod;
        return Poly(move(r));
    }

    Poly integ() const {
        if (size() == 0) return Poly();
        int n = size();
        vector<int> iv(n + 1), r(n + 1, 0);
        iv[1] = 1;
        for (int i = 2; i <= n; i++)
            iv[i] = 1ll * (mod - mod / i) * iv[mod % i] % mod;
        for (int i = 0; i < n; i++) r[i + 1] = 1ll * a[i] * iv[i + 1] % mod;
        return Poly(move(r));
    }
```

---

## 多项式求逆

求 $F \cdot G \equiv 1 \pmod{x^n}$。牛顿迭代 $G_{k+1} = G_k (2 - F G_k)$，要求 $a[0] \neq 0$。

```cpp
    Poly inv(int n) const {
        if (n <= 0) return Poly();
        Poly r({inv_mod(a[0])});
        for (int len = 2; len < (n << 1); len <<= 1) {
            int N = len << 1;
            vector<int> A(len);
            for (int i = 0; i < len; i++) A[i] = (*this)[i];
            vector<int> B = r.a;
            ntt(A, N, false);
            ntt(B, N, false);
            for (int i = 0; i < N; i++)
                B[i] = 1ll * B[i] * (2 - 1ll * A[i] * B[i] % mod + mod) % mod;
            ntt(B, N, true);
            B.resize(len);
            r = Poly(move(B));
        }
        r.resize(n);
        return r;
    }
```

---

## 多项式 ln / exp

$$
\ln F = \int \frac{F'}{F},\quad a[0] = 1
$$

$$
\exp F = G,\ \ln G \equiv F,\quad G_{k+1} = G_k(1 - \ln G_k + F),\quad a[0] = 0
$$

```cpp
    Poly ln(int n) const {
        Poly r = (deriv() * inv(n)).integ();
        r.resize(n);
        return r;
    }

    Poly exp(int n) const {
        if (n <= 0) return Poly();
        Poly r({1});
        for (int len = 2; len < (n << 1); len <<= 1) {
            Poly L = r.ln(len);
            vector<int> t(len);
            for (int i = 0; i < len; i++)
                t[i] = ((*this)[i] - L[i] + mod) % mod;
            t[0] = (t[0] + 1) % mod;
            r = r * Poly(move(t));
            r.resize(len);
        }
        r.resize(n);
        return r;
    }
```

---

## 多项式开根

求 $G^2 \equiv F \pmod{x^n}$。牛顿迭代 $G_{k+1} = (G_k + F G_k^{-1}) / 2$，要求 $a[0] = 1$（一般情形需二次剩余 / Cipolla 配合）。

```cpp
    Poly sqrt(int n) const {
        if (n <= 0) return Poly();
        const int i2 = (mod + 1) / 2;
        Poly r({1});
        for (int len = 2; len < (n << 1); len <<= 1) {
            vector<int> f(len);
            for (int i = 0; i < len; i++) f[i] = (*this)[i];
            Poly t = Poly(move(f)) * r.inv(len);
            t.resize(len);
            for (int i = 0; i < len; i++)
                t.a[i] = 1ll * (t[i] + r[i]) % mod * i2 % mod;
            r = move(t);
        }
        r.resize(n);
        return r;
    }
```

---

## 多项式快速幂

$$
F^k = a[0]^k \cdot \exp(k \cdot \ln(F / a[0]))
$$

`pow_core` 抽出 `(km, kp)` 共用主体，整数版与字符串版（用于 $k$ 极大）共用。$a[0] = 0$ 的情况通过提取 $x^{\text{shift}}$ 因子递归处理。

```cpp
private:
    // a[0] != 0：F^k = a0^kp * exp(km * ln(F / a0)),  km = k mod p, kp = k mod (p-1)
    Poly pow_core(int n, long long km, long long kp) const {
        int iv = inv_mod(a[0]);
        int a0k = pw(a[0], kp);
        Poly q; q.a.assign(size(), 0);
        for (int i = 0; i < size(); i++) q.a[i] = 1ll * a[i] * iv % mod;
        Poly L = q.ln(n);
        for (int& x : L.a) x = 1ll * x * km % mod;
        Poly r = L.exp(n);
        for (int& x : r.a) x = 1ll * x * a0k % mod;
        return r;
    }

public:
    // F^k, k 较小（fits long long）
    Poly pow(int n, long long k) const {
        if (n <= 0) return Poly();
        if (k == 0) { Poly r; r.a.assign(n, 0); r.a[0] = 1; return r; }
        int s = 0;
        while (s < size() && a[s] == 0) s++;
        if (s == size()) { Poly r; r.a.assign(n, 0); return r; }
        if (s > 0) {
            if ((__int128)s * k >= n) { Poly r; r.a.assign(n, 0); return r; }
            int sk = (int)(s * k);
            vector<int> q(size() - s);
            for (int i = 0; i + s < size(); i++) q[i] = a[i + s];
            Poly p = Poly(move(q)).pow(n - sk, k);
            Poly r; r.a.assign(n, 0);
            for (int i = 0; i < p.size() && i + sk < n; i++) r.a[i + sk] = p[i];
            return r;
        }
        return pow_core(n, k % mod, k % (mod - 1));
    }

    // F^k, k 极大（十进制字符串）
    Poly pow(int n, const string& s) const {
        if (n <= 0) return Poly();
        bool zero = true;
        for (char c : s) if (c != '0') { zero = false; break; }
        if (zero) { Poly r; r.a.assign(n, 0); r.a[0] = 1; return r; }
        long long km = 0, kp = 0;
        for (char c : s) {
            km = (km * 10 + (c - '0')) % mod;
            kp = (kp * 10 + (c - '0')) % (mod - 1);
        }
        int sh = 0;
        while (sh < size() && a[sh] == 0) sh++;
        if (sh == size() || sh > 0) { Poly r; r.a.assign(n, 0); return r; }
        return pow_core(n, km, kp);
    }
```

---

## 多项式除法 / 取模

设 $\deg F = n - 1,\ \deg G = m - 1$，求 $F = QG + R$，$\deg R < m - 1$。

记 $F^R(x) = x^{n-1} F(1/x)$（系数翻转），则

$$F^R(x) \equiv Q^R(x) \cdot G^R(x) \pmod{x^{n-m+1}}$$

故 $Q^R = F^R \cdot (G^R)^{-1} \bmod x^{n-m+1}$，反转得 $Q$，再 $R = F - QG$。

```cpp
    // 返回 (Q, R)：F = Q*G + R, deg R < deg G。要求 G 末项非零
    pair<Poly, Poly> divmod(const Poly& g) const {
        int n = size(), m = g.size();
        if (n < m) return {Poly(), *this};
        int qd = n - m + 1;
        vector<int> fr(a.rbegin(), a.rend()), gr(g.a.rbegin(), g.a.rend());
        fr.resize(qd); gr.resize(qd);
        Poly q = Poly(move(fr)) * Poly(move(gr)).inv(qd);
        q.resize(qd);
        reverse(q.a.begin(), q.a.end());
        Poly r = *this - q * g;
        r.resize(max(0, m - 1));
        return {q, r};
    }
    Poly operator/(const Poly& g) const { return divmod(g).first; }
    Poly operator%(const Poly& g) const { return divmod(g).second; }
```

---

## 多项式平移 f(x+c)

$$
f(x + c) = \sum_j x^j \cdot \frac{1}{j!} \sum_{i \ge j} (a_i \cdot i!) \cdot \frac{c^{i-j}}{(i-j)!}
$$

把 $A_i = a_i i!$ 翻转后与 $C_k = c^k / k!$ 卷积，结果在 $n - 1 - j$ 处除以 $j!$ 即得。整体 $O(n \log n)$。

```cpp
    // f(x + c)，长度保持 n
    Poly shift(int c) const {
        int n = size();
        if (n == 0) return Poly();
        c = ((c % mod) + mod) % mod;

        vector<int> fac(n), ifac(n);
        fac[0] = 1;
        for (int i = 1; i < n; i++) fac[i] = 1ll * fac[i - 1] * i % mod;
        ifac[n - 1] = inv_mod(fac[n - 1]);
        for (int i = n - 2; i >= 0; i--) ifac[i] = 1ll * ifac[i + 1] * (i + 1) % mod;

        vector<int> A(n), C(n);
        for (int i = 0; i < n; i++) A[n - 1 - i] = 1ll * a[i] * fac[i] % mod;
        int ck = 1;
        for (int k = 0; k < n; k++) {
            C[k] = 1ll * ck * ifac[k] % mod;
            ck   = 1ll * ck * c % mod;
        }
        auto h = mul_p(A, C, mod);
        vector<int> r(n);
        for (int j = 0; j < n; j++) r[j] = 1ll * h[n - 1 - j] * ifac[j] % mod;
        return Poly(move(r));
    }
```

---

## 多点求值

给定 $n$ 个 $x_i$，对度 $< n$ 的 $f$ 求 $f(x_i)$。在线段树上维护 $T_u(x) = \prod_{i \in u}(x - x_i)$，递归时令 $f \leftarrow f \bmod T_u$，叶子处即 $f(x_i)$。复杂度 $O(n \log^2 n)$。

```cpp
private:
    static void mpe_build(vector<Poly>& tr, int u, int l, int r,
                          const vector<int>& xs) {
        if (l == r) { tr[u] = Poly({(mod - xs[l]) % mod, 1}); return; }
        int m = (l + r) >> 1;
        mpe_build(tr, u << 1,     l,     m, xs);
        mpe_build(tr, u << 1 | 1, m + 1, r, xs);
        tr[u] = tr[u << 1] * tr[u << 1 | 1];
    }
    static void mpe_query(const vector<Poly>& tr, int u, int l, int r,
                          Poly f, vector<int>& res) {
        f = f % tr[u];
        if (l == r) { res[l] = f[0]; return; }
        int m = (l + r) >> 1;
        mpe_query(tr, u << 1,     l,     m, f, res);
        mpe_query(tr, u << 1 | 1, m + 1, r, f, res);
    }

public:
    vector<int> multi_eval(const vector<int>& xs) const {
        int n = xs.size();
        if (n == 0) return {};
        vector<Poly> tr(4 * n);
        mpe_build(tr, 1, 0, n - 1, xs);
        vector<int> res(n);
        mpe_query(tr, 1, 0, n - 1, *this, res);
        return res;
    }
```

---

## 快速插值

给定 $n$ 个 $(x_i, y_i)$（$x_i$ 互异），构造度 $< n$ 的 $f$。

$$
f(x) = \sum_i \frac{y_i}{M'(x_i)} \prod_{j \ne i} (x - x_j),\quad M(x) = \prod_i (x - x_i)
$$

用多点求值算 $d_i = M'(x_i)$，再分治合并 $\frac{y_i}{d_i}$ 与子树乘积。复杂度 $O(n \log^2 n)$。

```cpp
private:
    static Poly intp_solve(const vector<Poly>& tr, int u, int l, int r,
                           const vector<int>& ys, const vector<int>& d) {
        if (l == r) return Poly({(int)(1ll * ys[l] * inv_mod(d[l]) % mod)});
        int m = (l + r) >> 1;
        return intp_solve(tr, u << 1,     l,     m, ys, d) * tr[u << 1 | 1]
             + intp_solve(tr, u << 1 | 1, m + 1, r, ys, d) * tr[u << 1];
    }

public:
    static Poly fast_interp(const vector<int>& xs, const vector<int>& ys) {
        int n = xs.size();
        if (n == 0) return Poly();
        vector<Poly> tr(4 * n);
        mpe_build(tr, 1, 0, n - 1, xs);
        vector<int> d = tr[1].deriv().multi_eval(xs);
        return intp_solve(tr, 1, 0, n - 1, ys, d);
    }
```

---

## 常系数齐次线性递推 (Bostan–Mori)

设 $a_n = \sum_{i=1}^{k} c_i \cdot a_{n-i}$，给定初值 $a_0,\dots,a_{k-1}$，求 $a_n$（$n$ 可达 $10^{18}$）。

生成函数 $A(x) = P(x) / Q(x)$，$Q(x) = 1 - \sum_{i=1}^k c_i x^i$，$P = (A \cdot Q) \bmod x^k$。Bostan–Mori 利用

$$
\frac{P(x)}{Q(x)} = \frac{P(x) Q(-x)}{Q(x) Q(-x)},\quad Q(x) Q(-x) \text{ 仅含偶次项}
$$

每轮把 $n$ 减半，单轮 $O(k \log k)$，总 $O(k \log k \log n)$，远快于矩阵快速幂。

```cpp
    // a_n = sum_{i=1..k} c[i-1] * a_{n-i}, 给定 a[0..k-1], c[0..k-1]
    static int linear_recur(vector<int> a, vector<int> c, long long n) {
        int k = (int)c.size();
        if (n < (long long)a.size()) return a[(int)n];
        vector<int> Q(k + 1, 0); Q[0] = 1;
        for (int i = 0; i < k; i++) Q[i + 1] = (mod - c[i]) % mod;
        vector<int> P = mul_p(a, Q, mod);
        P.resize(k);
        while (n > 0) {
            vector<int> Qm = Q;
            for (int i = 1; i < (int)Qm.size(); i += 2)
                Qm[i] = (mod - Qm[i]) % mod;             // Q(-x)
            auto U = mul_p(P, Qm, mod);
            auto V = mul_p(Q, Qm, mod);                  // 仅偶次非零
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
        return 1ll * P[0] * inv_mod(Q[0]) % mod;
    }
```

---

## 单点求值

霍纳法 $O(n)$。

```cpp
    int eval(int x) const {
        int r = 0;
        for (int i = size() - 1; i >= 0; i--)
            r = (1ll * r * x + a[i]) % mod;
        return r;
    }
```

---

## 拉格朗日插值

### 一般插值

给定 $n$ 个点 $(x_i, y_i)$，构造度 $< n$ 的多项式。先构造 $M(x) = \prod (x - x_i)$，再合成除法求 $M(x) / (x - x_i)$，$O(n^2)$。$n$ 大时改用 [快速插值](#快速插值)。

```cpp
    static Poly lagrange(const vector<int>& xs, const vector<int>& ys) {
        int n = (int)xs.size();
        vector<int> M(n + 1, 0); M[0] = 1;
        for (int i = 0; i < n; i++) {
            int xi = xs[i];
            for (int k = i; k >= 0; k--) {
                M[k + 1] = (M[k + 1] + M[k]) % mod;
                M[k]     = (mod - 1ll * xi * M[k] % mod) % mod;
            }
        }
        Poly r; r.a.assign(n, 0);
        for (int i = 0; i < n; i++) {
            int xi = xs[i];
            vector<int> Q(n, 0);
            Q[n - 1] = M[n];
            for (int k = n - 1; k >= 1; k--)
                Q[k - 1] = (M[k] + 1ll * xi * Q[k]) % mod;
            int den = 1;
            for (int j = 0; j < n; j++) if (j != i)
                den = 1ll * den * (xi - xs[j] + mod) % mod;
            int co = 1ll * ys[i] * inv_mod(den) % mod;
            for (int k = 0; k < n; k++)
                r.a[k] = (r.a[k] + 1ll * co * Q[k]) % mod;
        }
        return r;
    }
```

### 连续点值插值

已知 $f(1), f(2), \dots, f(n)$，求 $f(k)$，$k$ 可达 $10^{18}$。

$$
f(k) = \sum_{i=1}^{n} y_{i-1} \cdot \frac{\text{pre}_{i-1} \cdot \text{suf}_{i+1}}{(i-1)! \cdot (-1)^{n-i} (n-i)!}
$$

线性预处理阶乘逆元和前后缀积，$O(n)$。

```cpp
    static int lagrange_continuous(const vector<int>& y, long long k) {
        int n = (int)y.size();
        if (k <= n) return y[(int)k - 1];
        vector<int> fac(n), ifac(n);
        fac[0] = 1;
        for (int i = 1; i < n; i++) fac[i] = 1ll * fac[i - 1] * i % mod;
        ifac[n - 1] = inv_mod(fac[n - 1]);
        for (int i = n - 2; i >= 0; i--) ifac[i] = 1ll * ifac[i + 1] * (i + 1) % mod;
        vector<int> pre(n + 1), suf(n + 2);
        pre[0] = suf[n + 1] = 1;
        for (int i = 1; i <= n; i++)
            pre[i] = 1ll * pre[i - 1] * ((k - i) % mod + mod) % mod;
        for (int i = n; i >= 1; i--)
            suf[i] = 1ll * suf[i + 1] * ((k - i) % mod + mod) % mod;
        int r = 0;
        for (int i = 1; i <= n; i++) {
            int t = 1ll * y[i - 1] * pre[i - 1] % mod * suf[i + 1] % mod
                          * ifac[i - 1] % mod * ifac[n - i] % mod;
            if ((n - i) & 1) t = (mod - t) % mod;
            r = (r + t) % mod;
        }
        return r;
    }
};
```

---

## 位运算卷积 FWT

定义集合上的卷积：

- OR  卷积：$h_S = \sum_{T \cup U = S} f_T g_U$
- AND 卷积：$h_S = \sum_{T \cap U = S} f_T g_U$
- XOR 卷积：$h_S = \sum_{T \oplus U = S} f_T g_U$

三者形状一致：变换 → 逐点乘 → 逆变换，单次变换 $O(n \log n)$。OR/AND 用前缀和/后缀和的子集和（Möbius）变换，XOR 用 Hadamard 矩阵。下标长度必须为 $2^n$。

```cpp
const int mod = 998244353;

// h[S] = sum_{T | U = S} f[T] g[U]
void fwt_or(vector<int>& a, bool inv) {
    int n = (int)a.size();
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int& x = a[i + j]; int& y = a[i + j + m];
                y = inv ? (y - x + mod) % mod : (y + x) % mod;
            }
}

// h[S] = sum_{T & U = S} f[T] g[U]
void fwt_and(vector<int>& a, bool inv) {
    int n = (int)a.size();
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int& x = a[i + j]; int& y = a[i + j + m];
                x = inv ? (x - y + mod) % mod : (x + y) % mod;
            }
}

// h[S] = sum_{T ^ U = S} f[T] g[U]
void fwt_xor(vector<int>& a, bool inv) {
    int n = (int)a.size();
    const int i2 = (mod + 1) / 2;
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int x = a[i + j], y = a[i + j + m];
                a[i + j]     = (x + y) % mod;
                a[i + j + m] = (x - y + mod) % mod;
                if (inv) {
                    a[i + j]     = 1ll * a[i + j]     * i2 % mod;
                    a[i + j + m] = 1ll * a[i + j + m] * i2 % mod;
                }
            }
}

template <class F>
vector<int> bit_conv(vector<int> A, vector<int> B, F fwt) {
    fwt(A, false); fwt(B, false);
    for (int i = 0; i < (int)A.size(); i++)
        A[i] = 1ll * A[i] * B[i] % mod;
    fwt(A, true);
    return A;
}
```

---

## 子集卷积

$$h_S = \sum_{T \subseteq S} f_T \cdot g_{S \setminus T}$$

OR 卷积只要求 $T \cup U = S$，不限制 $T \cap U = \emptyset$；用 popcount 分层即可强制不相交：仅当 $\text{popcount}(T) + \text{popcount}(U) = \text{popcount}(S)$ 时贡献到第 $\text{popcount}(S)$ 层。把 $f, g$ 拆成 $n+1$ 层做 OR 卷积再合并，复杂度 $O(2^n n^2)$。

```cpp
// 长度 = 2^n
vector<int> subset_conv(const vector<int>& f, const vector<int>& g) {
    int sz = (int)f.size(), n = __builtin_ctz(sz);
    vector<vector<int>> F(n + 1, vector<int>(sz, 0)),
                         G(n + 1, vector<int>(sz, 0));
    for (int s = 0; s < sz; s++) {
        int p = __builtin_popcount(s);
        F[p][s] = f[s];
        G[p][s] = g[s];
    }
    for (int i = 0; i <= n; i++) { fwt_or(F[i], false); fwt_or(G[i], false); }
    vector<vector<int>> H(n + 1, vector<int>(sz, 0));
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= i; j++)
            for (int s = 0; s < sz; s++)
                H[i][s] = (H[i][s] + 1ll * F[j][s] * G[i - j][s]) % mod;
    for (int i = 0; i <= n; i++) fwt_or(H[i], true);
    vector<int> r(sz);
    for (int s = 0; s < sz; s++) r[s] = H[__builtin_popcount(s)][s];
    return r;
}
```

---

## 复杂度总结

| 操作 | 复杂度 | 备注 |
|------|--------|------|
| `operator+` / `operator-` | $O(n)$ | 系数加减 |
| `operator*` | $O(n \log n)$ | NTT 卷积 |
| `mul_exact` | $O(n \log n)$ | 三模数精确卷积，结果 $< 2^{63}$ |
| `mul_any` | $O(n \log n)$ | 任意模数卷积，CRT 直接对 $p$ 取模 |
| `deriv` / `integ` | $O(n)$ | 求导 / 积分（积分用线性逆元） |
| `inv` | $O(n \log n)$ | 牛顿迭代，要求 $a[0] \ne 0$ |
| `ln` | $O(n \log n)$ | 要求 $a[0] = 1$ |
| `exp` | $O(n \log n)$ | 要求 $a[0] = 0$ |
| `sqrt` | $O(n \log n)$ | 要求 $a[0] = 1$ |
| `pow` | $O(n \log n)$ | 自动处理 $a[0] = 0$，支持 `string` 形式的极大 $k$ |
| `divmod` / `operator/` / `operator%` | $O(n \log n)$ | 依赖系数翻转 + 求逆 |
| `shift` | $O(n \log n)$ | $f(x + c)$，单次卷积 |
| `multi_eval` | $O(n \log^2 n)$ | 分治线段树 + 取模 |
| `fast_interp` | $O(n \log^2 n)$ | 多点求值 + 分治合并 |
| `linear_recur` | $O(k \log k \log n)$ | Bostan–Mori，$n \le 10^{18}$ |
| `eval` | $O(n)$ | 霍纳法 |
| `lagrange` | $O(n^2)$ | 一般插值，小规模 |
| `lagrange_continuous` | $O(n)$ | 连续点值，$k$ 任意 |
| `fwt_or` / `fwt_and` / `fwt_xor` | $O(n \log n)$ | 位运算卷积，$n = 2^k$ |
| `subset_conv` | $O(2^n n^2)$ | 不相交子集卷积 |
