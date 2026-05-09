- [NTT 与基础框架](#ntt-与基础框架)
- [三模数 NTT（精确卷积 / 任意模数）](#三模数-ntt精确卷积--任意模数)
- [基础运算](#基础运算)
  - [乘法](#乘法)
  - [求导 / 积分](#求导--积分)
- [多项式求逆](#多项式求逆)
- [多项式 ln / exp](#多项式-ln--exp)
- [多项式开根](#多项式开根)
- [多项式快速幂](#多项式快速幂)
- [单点求值](#单点求值)
- [拉格朗日插值](#拉格朗日插值)
  - [一般插值](#一般插值)
  - [连续点值插值](#连续点值插值)
- [复杂度总结](#复杂度总结)

---

## NTT 与基础框架

模数 $998244353 = 119 \times 2^{23} + 1$，原根 $g = 3$，支持长度 $\le 2^{23}$ 的 NTT。

```cpp
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Poly {
private:
    static const int mod = 998244353;
    static const int G = 3;

    static int qpow(int a, int k) {
        int res = 1;
        for (; k; a = 1ll * a * a % mod, k >>= 1) {
            if (k & 1) res = 1ll * res * a % mod;
        }
        return res;
    }

    static vector<int> get_rev(int limit) {
        vector<int> rev(limit);
        int bit = 0;
        while ((1 << bit) < limit) bit++;
        for (int i = 0; i < limit; i++) {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
        }
        return rev;
    }

    static void NTT(vector<int>& a, int limit, bool type) {
        a.resize(limit, 0);
        vector<int> rev = get_rev(limit);
        for (int i = 0; i < limit; i++) {
            if (i < rev[i]) swap(a[i], a[rev[i]]);
        }
        for (int mid = 1; mid < limit; mid <<= 1) {
            int wn = qpow(type ? G : qpow(G, mod - 2), (mod - 1) / (mid << 1));
            for (int j = 0; j < limit; j += (mid << 1)) {
                int w = 1;
                for (int k = 0; k < mid; k++, w = 1ll * w * wn % mod) {
                    int x = a[j + k], y = 1ll * w * a[j + k + mid] % mod;
                    a[j + k] = (x + y >= mod) ? (x + y - mod) : (x + y);
                    a[j + k + mid] = (x - y < 0) ? (x - y + mod) : (x - y);
                }
            }
        }
        if (!type) {
            int inv_limit = qpow(limit, mod - 2);
            for (int i = 0; i < limit; i++) {
                a[i] = 1ll * a[i] * inv_limit % mod;
            }
        }
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

```cpp
public:
    vector<int> a;

    Poly() {}
    Poly(const vector<int>& _a) : a(_a) {}
    Poly(int val) { a.push_back(val); }

    int size() const { return a.size(); }
    void resize(int n) { a.resize(n, 0); }
    int operator[](int idx) const { return idx < size() ? a[idx] : 0; }
```

### 乘法

卷积 $C = A \times B$，结果长度为 $n + m - 1$。

```cpp
    Poly operator*(const Poly& rhs) const {
        if (size() == 0 || rhs.size() == 0) return Poly();
        int n = size(), m = rhs.size();
        int limit = 1;
        while (limit < n + m - 1) limit <<= 1;

        vector<int> A = a, B = rhs.a;
        NTT(A, limit, true);
        NTT(B, limit, true);
        for (int i = 0; i < limit; i++) A[i] = 1ll * A[i] * B[i] % mod;
        NTT(A, limit, false);

        A.resize(n + m - 1);
        return Poly(A);
    }
```

### 求导 / 积分

$$
\frac{d}{dx} \sum a_i x^i = \sum (i+1)a_{i+1} x^i
$$

$$
\int \sum a_i x^i = \sum \frac{a_i}{i+1} x^{i+1}
$$

```cpp
    Poly deriv() const {
        if (size() <= 1) return Poly(vector<int>{0});
        vector<int> res(size() - 1);
        for (int i = 1; i < size(); i++) {
            res[i - 1] = 1ll * a[i] * i % mod;
        }
        return Poly(res);
    }

    Poly integ() const {
        if (size() == 0) return Poly();
        vector<int> res(size() + 1, 0);
        for (int i = 0; i < size(); i++) {
            res[i + 1] = 1ll * a[i] * qpow(i + 1, mod - 2) % mod;
        }
        return Poly(res);
    }
```

---

## 多项式求逆

求 $F \cdot G \equiv 1 \pmod{x^n}$。

牛顿迭代：$G_{k+1} = G_k (2 - F \cdot G_k)$，要求 $a[0] \neq 0$。

```cpp
    Poly inv(int deg) const {
        if (deg <= 0) return Poly();
        Poly res(vector<int>{qpow(a[0], mod - 2)});
        for (int len = 2; len < (deg << 1); len <<= 1) {
            int limit = len << 1;
            vector<int> A(len);
            for(int i = 0; i < len; i++) A[i] = (*this)[i];

            vector<int> B = res.a;
            NTT(A, limit, true);
            NTT(B, limit, true);
            for (int i = 0; i < limit; i++) {
                B[i] = 1ll * B[i] * (2 - 1ll * A[i] * B[i] % mod + mod) % mod;
            }
            NTT(B, limit, false);
            B.resize(len);
            res = Poly(B);
        }
        res.resize(deg);
        return res;
    }
```

---

## 多项式 ln / exp

$$
\ln F = \int \frac{F'}{F}, \quad a[0] = 1
$$

$$
\exp F = G,\ \ln G - F = 0,\quad G_{k+1} = G_k (1 - \ln G_k + F), \quad a[0] = 0
$$

```cpp
    Poly ln(int deg) const {
        // 要求 a[0] == 1
        Poly res = (this->deriv() * this->inv(deg)).integ();
        res.resize(deg);
        return res;
    }

    Poly exp(int deg) const {
        if (deg <= 0) return Poly();
        // 要求 a[0] == 0
        Poly res(vector<int>{1});
        for (int len = 2; len < (deg << 1); len <<= 1) {
            Poly lnx = res.ln(len);
            vector<int> A(len);
            for(int i = 0; i < len; i++) A[i] = (*this)[i];

            vector<int> tmp = lnx.a;
            tmp.resize(len, 0);
            for(int i = 0; i < len; i++) {
                tmp[i] = (A[i] - tmp[i] + mod) % mod;
            }
            tmp[0] = (tmp[0] + 1) % mod;

            res = res * Poly(tmp);
            res.resize(len);
        }
        res.resize(deg);
        return res;
    }
```

---

## 多项式开根

求 $G^2 \equiv F \pmod{x^n}$。

牛顿迭代：$G_{k+1} = \dfrac{G_k + F \cdot G_k^{-1}}{2}$，要求 $a[0] = 1$。

```cpp
    // 多项式开根，要求 a[0] == 1
    Poly sqrt(int deg) const {
        if (deg <= 0) return Poly();
        Poly res(vector<int>{1});
        for (int len = 2; len < (deg << 1); len <<= 1) {
            vector<int> f_vec(len, 0);
            for (int i = 0; i < len; i++) f_vec[i] = (*this)[i];
            Poly F(f_vec);

            Poly inv_res = res.inv(len);
            Poly tmp = F * inv_res;
            tmp.resize(len);
            for (int i = 0; i < len; i++) {
                tmp.a[i] = (res[i] + tmp.a[i]) % mod;
                tmp.a[i] = 1ll * tmp.a[i] * ((mod + 1) / 2) % mod;
            }
            res = tmp;
        }
        res.resize(deg);
        return res;
    }
```

---

## 多项式快速幂

$$
F^k = a[0]^k \cdot \exp(k \cdot \ln(F / a[0]))
$$

自动处理 $a[0] = 0$ 的情况（提取 $x^{\text{shift}}$ 因子）。

```cpp
    // 多项式快速幂 F^k (k 较小时)
    Poly pow(int deg, long long k) const {
        if (deg <= 0) return Poly();
        if (k == 0) {
            Poly res; res.a.resize(deg, 0); res.a[0] = 1;
            return res;
        }

        int shift = 0;
        while (shift < size() && a[shift] == 0) shift++;
        if (shift == size()) { Poly res; res.a.resize(deg, 0); return res; }

        if (shift > 0) {
            if (k >= (deg + shift - 1) / shift) {
                Poly res; res.a.resize(deg, 0);
                return res;
            }
            int new_deg = deg - shift * k;
            vector<int> shifted_vec(size() - shift);
            for (int i = 0; i < size() - shift; i++) shifted_vec[i] = a[i + shift];
            Poly shifted(shifted_vec);
            Poly powed = shifted.pow(new_deg, k);

            vector<int> res_vec(deg, 0);
            for (int i = 0; i < powed.size() && i + shift * k < deg; i++) {
                res_vec[i + shift * k] = powed.a[i];
            }
            return Poly(res_vec);
        }

        int inv_a0 = qpow(a[0], mod - 2);
        int a0_k = qpow(a[0], k % (mod - 1));

        Poly norm;
        norm.a.resize(size());
        for (int i = 0; i < size(); i++) norm.a[i] = 1ll * a[i] * inv_a0 % mod;

        Poly lnx = norm.ln(deg);
        for (int i = 0; i < lnx.size() && i < deg; i++) {
            lnx.a[i] = 1ll * lnx.a[i] * (k % mod) % mod;
        }
        Poly res = lnx.exp(deg);
        for (int i = 0; i < res.size(); i++) {
            res.a[i] = 1ll * res.a[i] * a0_k % mod;
        }
        return res;
    }

    // 多项式快速幂 F^k (k 极大, 字符串表示)
    Poly pow(int deg, const string& s) const {
        if (s.size() <= 18 || (s.size() == 19 && s <= "9223372036854775807")) {
            return pow(deg, stoll(s));
        }
        int km = 0, kp = 0;
        for (char c : s) {
            km = (10ll * km + (c - '0')) % mod;
            kp = (10ll * kp + (c - '0')) % (mod - 1);
        }
        if (s == "0") {
            Poly res; res.a.resize(deg, 0); res.a[0] = 1;
            return res;
        }

        int shift = 0;
        while (shift < size() && a[shift] == 0) shift++;
        if (shift == size()) { Poly res; res.a.resize(deg, 0); return res; }
        if (shift > 0) { Poly res; res.a.resize(deg, 0); return res; }

        int inv_a0 = qpow(a[0], mod - 2);
        int a0_k = qpow(a[0], kp);

        Poly norm;
        norm.a.resize(size());
        for (int i = 0; i < size(); i++) norm.a[i] = 1ll * a[i] * inv_a0 % mod;

        Poly lnx = norm.ln(deg);
        for (int i = 0; i < lnx.size() && i < deg; i++) {
            lnx.a[i] = 1ll * lnx.a[i] * km % mod;
        }
        Poly res = lnx.exp(deg);
        for (int i = 0; i < res.size(); i++) {
            res.a[i] = 1ll * res.a[i] * a0_k % mod;
        }
        return res;
    }
```

---

## 单点求值

霍纳法（Horner）求 $f(x)$，$O(n)$。

```cpp
    // 单点求值 (Horner)
    int eval(int x) const {
        int res = 0;
        for (int i = size() - 1; i >= 0; i--) {
            res = (1ll * res * x + a[i]) % mod;
        }
        return res;
    }
```

---

## 拉格朗日插值

### 一般插值

给定 $n$ 个点 $(x_i, y_i)$，构造 $n-1$ 次多项式。

$$
P(x) = \sum y_i \cdot \prod_{j \ne i} \frac{x - x_j}{x_i - x_j}
$$

实现：先构造 $M(x) = \prod (x - x_i)$，再对每个 $i$ 用合成除法求 $M(x) / (x - x_i)$。

```cpp
    static Poly lagrange(const vector<int>& xs, const vector<int>& ys) {
        int n = xs.size();
        vector<int> M(n + 1, 0);
        M[0] = 1;
        for (int i = 0; i < n; i++) {
            int xi = xs[i];
            for (int k = i; k >= 0; k--) {
                M[k+1] = (M[k+1] + M[k]) % mod;
                M[k] = (mod - 1ll * xi * M[k] % mod) % mod;
            }
        }

        Poly ans(vector<int>(n, 0));
        for (int i = 0; i < n; i++) {
            int xi = xs[i];
            vector<int> Q(n, 0);
            Q[n-1] = M[n];
            for (int k = n-1; k >= 1; k--) {
                Q[k-1] = (M[k] + 1ll * xi * Q[k]) % mod;
            }
            int den = 1;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                den = 1ll * den * (xi - xs[j] + mod) % mod;
            }
            int coeff = 1ll * ys[i] * qpow(den, mod - 2) % mod;
            for (int k = 0; k < n; k++) {
                ans.a[k] = (ans.a[k] + 1ll * coeff * Q[k]) % mod;
            }
        }
        return ans;
    }
```

### 连续点值插值

已知 $f(1), f(2), \dots, f(n)$，求 $f(k)$。

$$
f(k) = \sum_{i=1}^{n} y_{i-1} \cdot \frac{\text{pref}_{i-1} \cdot \text{suff}_{i+1}}{(i-1)! \cdot (-1)^{n-i} \cdot (n-i)!}
$$

$O(n)$ 预处理阶乘和前后缀积，$O(n)$ 计算。支持 $k$ 到 $10^{18}$ 取模。

```cpp
    static int lagrange_continuous(const vector<int>& y, long long k) {
        int n = y.size();
        if (k <= n) return y[k - 1];

        vector<int> fact(n), inv_fact(n);
        fact[0] = 1;
        for (int i = 1; i < n; i++) fact[i] = 1ll * fact[i-1] * i % mod;
        inv_fact[n-1] = qpow(fact[n-1], mod - 2);
        for (int i = n-2; i >= 0; i--) inv_fact[i] = 1ll * inv_fact[i+1] * (i+1) % mod;

        vector<int> pref(n + 1);
        pref[0] = 1;
        for (int i = 1; i <= n; i++)
            pref[i] = 1ll * pref[i-1] * ((k - i) % mod + mod) % mod;

        vector<int> suff(n + 2);
        suff[n+1] = 1;
        for (int i = n; i >= 1; i--)
            suff[i] = 1ll * suff[i+1] * ((k - i) % mod + mod) % mod;

        int ans = 0;
        for (int i = 1; i <= n; i++) {
            int num = 1ll * pref[i-1] * suff[i+1] % mod;
            int den = 1ll * inv_fact[i-1] * inv_fact[n-i] % mod;
            int term = 1ll * y[i-1] * num % mod * den % mod;
            if ((n - i) & 1) term = mod - term;
            ans = (ans + term) % mod;
        }
        return ans;
    }
};
```

---

## 复杂度总结

| 操作 | 复杂度 | 备注 |
|------|--------|------|
| `operator*` | $O(n \log n)$ | NTT 卷积 |
| `conv_exact` | $O(n \log n)$ | 三模数精确卷积 |
| `conv_arbitrary` | $O(n \log n)$ | 任意模数卷积 |
| `deriv` / `integ` | $O(n)$ | 求导 / 积分 |
| `inv` | $O(n \log n)$ | 牛顿迭代 |
| `ln` | $O(n \log n)$ | $a[0] = 1$ |
| `exp` | $O(n \log n)$ | $a[0] = 0$ |
| `sqrt` | $O(n \log n)$ | $a[0] = 1$ |
| `pow` | $O(n \log n)$ | 自动处理 $a[0]=0$ |
| `eval` | $O(n)$ | 霍纳法 |
| `lagrange` | $O(n^2)$ | 一般插值 |
| `lagrange_continuous` | $O(n)$ | 连续点值 |
