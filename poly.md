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

- `pw` 接受任意素数模 `p`，默认 `mod`。
- `ntt` 用 `bool inv` 表示逆变换（true=逆），`p` 为模数参数，三模数共用。
- `rev` 用 `static` 按长度缓存，避免每次卷积重建。

```cpp
#include <bits/stdc++.h>
using namespace std;

struct Poly {
    static const int mod = 998244353, G = 3, G_inv = 332748118;

    vector<int> a;

    Poly() = default;
    Poly(vector<int> v) : a(move(v)) {}

    int size() const { return (int)a.size(); }
    void resize(int n) { a.resize(n, 0); }
    int& operator[](int i) { return a[i]; }
    int operator[](int i) const { return i < size() ? a[i] : 0; }

    static int pw(int x, long long k, int p = mod) {
        int r = 1; x = (x % p + p) % p;
        for (; k; x = 1ll * x * x % p, k >>= 1)
            if (k & 1) r = 1ll * r * x % p;
        return r;
    }
    static int inv_mod(int x, int p = mod) { return pw(x, p - 2, p); }

    static void ntt(vector<int> &a, int n, bool inv, int p = mod) {
        a.resize(n, 0);
        static vector<int> rev; static int last = -1;
        if (last != n) {
            rev.resize(n); last = n;
            for (int i = 1; i < n; i++)
                rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? n >> 1 : 0);
        }
        for (int i = 0; i < n; i++)
            if (i < rev[i]) swap(a[i], a[rev[i]]);
        for (int m = 1; m < n; m <<= 1) {
            int wn = pw(inv ? G_inv : G, (p - 1) / (m << 1), p);
            for (int i = 0; i < n; i += m << 1) {
                int w = 1;
                for (int j = 0; j < m; j++, w = 1ll * w * wn % p) {
                    int x = a[i + j], y = 1ll * w * a[i + j + m] % p;
                    a[i + j] = (x + y) >= p ? x + y - p : x + y;
                    a[i + j + m] = (x - y + p) >= p ? x - y : x - y + p;
                }
            }
        }
        if (inv) {
            int iv = inv_mod(n, p);
            for (int &x : a) x = 1ll * x * iv % p;
        }
    }

    static vector<int> mul_p(vector<int> A, vector<int> B, int p = mod) {
        if (A.empty() || B.empty()) return {};
        int len = (int)A.size() + (int)B.size() - 1, n = 1;
        while (n < len) n <<= 1;
        ntt(A, n, false, p); ntt(B, n, false, p);
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
    static const int P1 = 998244353, P2 = 1004535809, P3 = 469762049;

    // CRT 合并三个模数 → 对 p 取模
    static int crt_mod(int r1, int r2, int r3, int p) {
        static const int i12 = inv_mod(P1 % P2, P2);
        static const long long M12 = 1ll * P1 * P2;
        static const int i123 = inv_mod((int)(M12 % P3), P3);
        long long t1 = 1ll * ((r2 - r1) % P2 + P2) % P2 * i12 % P2;
        long long x = r1 + P1 * t1;
        long long t2 = ((r3 - x % P3) % P3 + P3) % P3 * i123 % P3;
        return (int)((x % p + (M12 % p) * t2 % p) % p);
    }

    // CRT 合并三个模数 → 精确整数（结果 < 2^63）
    static long long crt(ll r1, ll r2, ll r3) {
        static const int i12 = inv_mod(P1 % P2, P2);
        static const long long M12 = 1ll * P1 * P2;
        static const int i123 = inv_mod((int)(M12 % P3), P3);
        long long t1 = 1ll * ((r2 - r1) % P2 + P2) % P2 * i12 % P2;
        long long x = r1 + P1 * t1;
        long long t2 = ((r3 - x % P3) % P3 + P3) % P3 * i123 % P3;
        return x + M12 * t2;
    }

    static vector<int> mul_any(const vector<int> &a, const vector<int> &b, int p) {
        auto c1 = mul_p(a, b, P1), c2 = mul_p(a, b, P2), c3 = mul_p(a, b, P3);
        vector<int> r(c1.size());
        for (int i = 0; i < (int)r.size(); i++)
            r[i] = crt_mod(c1[i], c2[i], c3[i], p);
        return r;
    }

    static vector<long long> mul_exact(const vector<int> &a, const vector<int> &b) {
        auto c1 = mul_p(a, b, P1), c2 = mul_p(a, b, P2), c3 = mul_p(a, b, P3);
        vector<long long> r(c1.size());
        for (int i = 0; i < (int)r.size(); i++)
            r[i] = crt(c1[i], c2[i], c3[i]);
        return r;
    }
```


---

## 基础运算

### 加 / 减 / 乘

```cpp
    Poly operator+(const Poly &b) const {
        Poly r; r.a.resize(max(size(), b.size()));
        for (int i = 0; i < r.size(); i++) r[i] = ((*this)[i] + b[i]) % mod;
        return r;
    }
    Poly operator-(const Poly &b) const {
        Poly r; r.a.resize(max(size(), b.size()));
        for (int i = 0; i < r.size(); i++) r[i] = ((*this)[i] - b[i] + mod) % mod;
        return r;
    }
    Poly operator*(const Poly &b) const { return Poly(mul_p(a, b.a)); }
    Poly &operator*=(const Poly &b) { return *this = *this * b; }
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
    Poly pow_core(int n, long long km, long long kp) const {
        int iv = inv_mod(a[0]), a0k = pw(a[0], kp);
        Poly q; q.a.resize(size());
        for (int i = 0; i < size(); i++) q[i] = 1ll * a[i] * iv % mod;
        q = q.ln(n);
        for (int &x : q.a) x = 1ll * x * km % mod;
        q = q.exp(n);
        for (int &x : q.a) x = 1ll * x * a0k % mod;
        return q;
    }

    // k 较小时
    Poly pow(int n, long long k) const {
        if (n <= 0) return Poly();
        if (k == 0) { Poly r; r.a.assign(n, 0); r[0] = 1; return r; }
        int s = 0;
        while (s < size() && a[s] == 0) s++;
        if (s == size()) { Poly r; r.a.assign(n, 0); return r; }
        if (s > 0) {
            if ((__int128)s * k >= n) { Poly r; r.a.assign(n, 0); return r; }
            int sk = (int)(s * k);
            vector<int> q(a.begin() + s, a.end());
            Poly p = Poly(move(q)).pow(n - sk, k);
            Poly r; r.a.assign(n, 0);
            for (int i = 0; i < p.size() && i + sk < n; i++) r[i + sk] = p[i];
            return r;
        }
        return pow_core(n, k % mod, k % (mod - 1));
    }

    // k 极大（字符串）
    Poly pow(int n, const string &s) const {
        if (n <= 0) return Poly();
        bool zero = true;
        for (char c : s) if (c != '0') { zero = false; break; }
        if (zero) { Poly r; r.a.assign(n, 0); r[0] = 1; return r; }
        long long km = 0, kp = 0;
        for (char c : s) {
            km = (km * 10 + c - '0') % mod;
            kp = (kp * 10 + c - '0') % (mod - 1);
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
    pair<Poly, Poly> divmod(const Poly &g) const {
        int n = size(), m = g.size();
        if (n < m) return {Poly(), *this};
        int qd = n - m + 1;
        vector<int> fr(a.rbegin(), a.rend()), gr(g.a.rbegin(), g.a.rend());
        fr.resize(qd); gr.resize(qd);
        Poly q = Poly(move(fr)) * Poly(move(gr)).inv(qd);
        q.resize(qd);
        reverse(q.a.begin(), q.a.end());
        Poly r = *this - q * g;
        r.resize(m - 1);
        return {q, r};
    }
    Poly operator/(const Poly &g) const { return divmod(g).first; }
    Poly operator%(const Poly &g) const { return divmod(g).second; }
```

---

## 多项式平移 f(x+c)

$$
f(x + c) = \sum_j x^j \cdot \frac{1}{j!} \sum_{i \ge j} (a_i \cdot i!) \cdot \frac{c^{i-j}}{(i-j)!}
$$

把 $A_i = a_i i!$ 翻转后与 $C_k = c^k / k!$ 卷积，结果在 $n - 1 - j$ 处除以 $j!$ 即得。整体 $O(n \log n)$。

```cpp
    // f(x + c)
    Poly shift(int c) const {
        int n = size();
        if (!n) return Poly();
        c = (c % mod + mod) % mod;
        vector<int> fac(n), ifac(n), A(n), C(n), r(n);
        fac[0] = 1;
        for (int i = 1; i < n; i++) fac[i] = 1ll * fac[i - 1] * i % mod;
        ifac[n - 1] = inv_mod(fac[n - 1]);
        for (int i = n - 2; i >= 0; i--) ifac[i] = 1ll * ifac[i + 1] * (i + 1) % mod;
        for (int i = 0; i < n; i++) A[n - 1 - i] = 1ll * a[i] * fac[i] % mod;
        for (int i = 0, ck = 1; i < n; i++, ck = 1ll * ck * c % mod)
            C[i] = 1ll * ck * ifac[i] % mod;
        auto h = mul_p(A, C);
        for (int j = 0; j < n; j++) r[j] = 1ll * h[n - 1 - j] * ifac[j] % mod;
        return Poly(move(r));
    }
```

---

## 多点求值

给定 $n$ 个 $x_i$，对度 $< n$ 的 $f$ 求 $f(x_i)$。在线段树上维护 $T_u(x) = \prod_{i \in u}(x - x_i)$，递归时令 $f \leftarrow f \bmod T_u$，叶子处即 $f(x_i)$。复杂度 $O(n \log^2 n)$。

```cpp
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

    vector<int> multi_eval(const vector<int> &xs) const {
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
    static Poly intp_solve(const vector<Poly>& tr, int u, int l, int r,
                           const vector<int>& ys, const vector<int>& d) {
        if (l == r) return Poly({(int)(1ll * ys[l] * inv_mod(d[l]) % mod)});
        int m = (l + r) >> 1;
        return intp_solve(tr, u << 1,     l,     m, ys, d) * tr[u << 1 | 1]
             + intp_solve(tr, u << 1 | 1, m + 1, r, ys, d) * tr[u << 1];
    }

    static Poly fast_interp(const vector<int> &xs, const vector<int> &ys) {
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
