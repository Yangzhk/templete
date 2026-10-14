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
vector<int> stirling1_unsigned(int n) {
    assert(0 <= n && n < mod);
    if (n == 0) return {1};
    if (n == 1) return {0, 1};
    int m = n / 2;
    auto a = stirling1_unsigned(m);
    auto b = poly_shift(a, m);
    auto c = mul_p(a, b);
    if (n & 1) c = mul_p(c, {n - 1, 1});
    c.resize(n + 1);
    return c;
}
vector<int> stirling1_signed(int n) {
    auto a = stirling1_unsigned(n);
    for (int k = 0; k <= n; k++)
        if ((n - k) & 1) a[k] = a[k] ? mod - a[k] : 0;
    return a;
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

// FWT：数组长度必须是 2 的幂，全部在 mod 下运算
void fwt_or(vector<int> &a, bool inv) {
    int n = a.size();
    assert(n > 0 && (n & (n - 1)) == 0);
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int &x = a[i + j], &y = a[i + j + m];
                y = inv ? (y - x + mod) % mod : (y + x) % mod;
            }
}
void fwt_and(vector<int> &a, bool inv) {
    int n = a.size();
    assert(n > 0 && (n & (n - 1)) == 0);
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int &x = a[i + j], &y = a[i + j + m];
                x = inv ? (x - y + mod) % mod : (x + y) % mod;
            }
}
void fwt_xor(vector<int> &a, bool inv) {
    int n = a.size();
    assert(n > 0 && (n & (n - 1)) == 0);
    int inv2 = (mod + 1) / 2;
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int x = a[i + j], y = a[i + j + m];
                a[i + j] = (x + y) % mod;
                a[i + j + m] = (x - y + mod) % mod;
                if (inv) {
                    a[i + j] = 1LL * a[i + j] * inv2 % mod;
                    a[i + j + m] = 1LL * a[i + j + m] * inv2 % mod;
                }
            }
}
template<class F>
vector<int> bit_conv(vector<int> a, vector<int> b, F fwt) {
    assert(a.size() == b.size());
    fwt(a, false); fwt(b, false);
    for (int i = 0; i < (int)a.size(); i++) a[i] = 1LL * a[i] * b[i] % mod;
    fwt(a, true);
    return a;
}

// 子集卷积：h[S] = sum_{T subseteq S} f[T]*g[S\T]，O(n^2*2^n)
vector<int> subset_conv(const vector<int> &f, const vector<int> &g) {
    int sz = f.size();
    assert(sz > 0 && sz == (int)g.size() && (sz & (sz - 1)) == 0);
    int n = __builtin_ctz((unsigned)sz);
    vector<vector<int>> F(n + 1, vector<int>(sz)), GG(n + 1, vector<int>(sz));
    for (int s = 0; s < sz; s++) {
        int p = __builtin_popcount((unsigned)s);
        F[p][s] = f[s]; GG[p][s] = g[s];
    }
    for (int i = 0; i <= n; i++) { fwt_or(F[i], false); fwt_or(GG[i], false); }
    vector<vector<int>> H(n + 1, vector<int>(sz));
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= i; j++)
            for (int s = 0; s < sz; s++)
                H[i][s] = (H[i][s] + 1LL * F[j][s] * GG[i - j][s]) % mod;
    for (int i = 0; i <= n; i++) fwt_or(H[i], true);
    vector<int> r(sz);
    for (int s = 0; s < sz; s++) r[s] = H[__builtin_popcount((unsigned)s)][s];
    return r;
}
```
