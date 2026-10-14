```
#include <bits/stdc++.h>
using namespace std;
#define int long long

const int mod = 998244353, G = 3;
const int P1 = 998244353, P2 = 1004535809, P3 = 469762049;

int norm(int x, int p = mod) {
    x %= p;
    return x < 0 ? x + p : x;
}

int pw(int x, int k, int p = mod) {
    int r = 1;
    x = norm(x, p);
    for (; k; k >>= 1, x = x * x % p)
        if (k & 1) r = r * x % p;
    return r;
}

int inv_mod(int x, int p = mod) { return pw(x, p - 2, p); }

void ntt(vector<int> &a, bool inv, int p = mod) {
    int n = a.size();
    assert(n > 0 && (n & (n - 1)) == 0);
    assert((p - 1) % n == 0);
    vector<int> rev(n);
    for (int i = 1; i < n; i++)
        rev[i] = (rev[i >> 1] >> 1) | ((i & 1) ? n >> 1 : 0);
    for (int i = 0; i < n; i++)
        if (i < rev[i]) swap(a[i], a[rev[i]]);
    for (int m = 1; m < n; m <<= 1) {
        int w0 = pw(G, (p - 1) / (m << 1), p);
        if (inv) w0 = inv_mod(w0, p);
        for (int i = 0; i < n; i += m << 1) {
            int w = 1;
            for (int j = 0; j < m; j++, w = w * w0 % p) {
                int x = a[i + j], y = w * a[i + j + m] % p;
                a[i + j] = x + y >= p ? x + y - p : x + y;
                a[i + j + m] = x - y < 0 ? x - y + p : x - y;
            }
        }
    }
    if (inv) {
        int t = inv_mod(n, p);
        for (int &x : a) x = x * t % p;
    }
}

vector<int> mul_p(vector<int> a, vector<int> b, int p = mod) {
    if (a.empty() || b.empty()) return {};
    int len = a.size() + b.size() - 1, n = 1;
    while (n < len) n <<= 1;
    assert(p == P1 || p == P2 || p == P3);
    assert((p - 1) % n == 0);
    for (int &x : a) x = norm(x, p);
    for (int &x : b) x = norm(x, p);
    a.resize(n); b.resize(n);
    ntt(a, false, p); ntt(b, false, p);
    for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % p;
    ntt(a, true, p);
    a.resize(len);
    return a;
}

const int M12 = P1 * P2;
const int I12 = inv_mod(P1 % P2, P2);
const int I123 = inv_mod(M12 % P3, P3);

int crt_mod(int r1, int r2, int r3, int p) {
    int t1 = norm(r2 - r1, P2) * I12 % P2;
    int x = r1 + P1 * t1;
    int t2 = norm(r3 - x % P3, P3) * I123 % P3;
    return ((x % p) + (__int128)(M12 % p) * t2) % p;
}

int crt(int r1, int r2, int r3) {
    int t1 = norm(r2 - r1, P2) * I12 % P2;
    int x = r1 + P1 * t1;
    int t2 = norm(r3 - x % P3, P3) * I123 % P3;
    __int128 ans = (__int128)x + (__int128)M12 * t2;
    assert(ans <= LLONG_MAX);
    return (int)ans;
}

vector<int> mul_any(vector<int> a, vector<int> b, int p) {
    if (a.empty() || b.empty()) return {};
    assert(p > 0);
    for (int &x : a) x = norm(x, p);
    for (int &x : b) x = norm(x, p);
    auto c1 = mul_p(a, b, P1), c2 = mul_p(a, b, P2), c3 = mul_p(a, b, P3);
    vector<int> res(c1.size());
    for (int i = 0; i < (int)res.size(); i++)
        res[i] = crt_mod(c1[i], c2[i], c3[i], p);
    return res;
}

vector<int> mul_exact(const vector<int> &a, const vector<int> &b) {
    if (a.empty() || b.empty()) return {};
    auto c1 = mul_p(a, b, P1), c2 = mul_p(a, b, P2), c3 = mul_p(a, b, P3);
    vector<int> res(c1.size());
    for (int i = 0; i < (int)res.size(); i++)
        res[i] = crt(c1[i], c2[i], c3[i]);
    return res;
}

struct Poly {
    vector<int> a;
    Poly(vector<int> v = {}) : a(move(v)) {
        for (int &x : a) x = norm(x);
    }
    int size() const { return a.size(); }
    void resize(int n) { a.resize(n); }
    int &operator[](int i) { return a[i]; }
    int operator[](int i) const { return i >= 0 && i < size() ? a[i] : 0; }
};

Poly operator+(const Poly &f, const Poly &g) {
    Poly h;
    h.resize(max(f.size(), g.size()));
    for (int i = 0; i < h.size(); i++) h[i] = (f[i] + g[i]) % mod;
    return h;
}
Poly operator-(const Poly &f, const Poly &g) {
    Poly h;
    h.resize(max(f.size(), g.size()));
    for (int i = 0; i < h.size(); i++) h[i] = norm(f[i] - g[i]);
    return h;
}
Poly operator*(const Poly &f, const Poly &g) { return Poly(mul_p(f.a, g.a)); }
Poly &operator*=(Poly &f, const Poly &g) { return f = f * g; }

Poly deriv(const Poly &f) {
    if (f.size() <= 1) return Poly({0});
    vector<int> a(f.size() - 1);
    for (int i = 1; i < f.size(); i++) a[i - 1] = f[i] * i % mod;
    return Poly(move(a));
}

Poly integ(const Poly &f) {
    if (!f.size()) return Poly();
    int n = f.size();
    assert(n < mod);
    vector<int> iv(n + 1), a(n + 1);
    iv[1] = 1;
    for (int i = 2; i <= n; i++)
        iv[i] = (mod - mod / i) * iv[mod % i] % mod;
    for (int i = 0; i < n; i++) a[i + 1] = f[i] * iv[i + 1] % mod;
    return Poly(move(a));
}

Poly poly_inv(const Poly &f, int n) {
    if (n <= 0) return Poly();
    assert(f.size() && f[0] != 0);
    Poly g({inv_mod(f[0])});
    for (int m = 2; m < 2 * n; m <<= 1) {
        int N = m << 1;
        vector<int> a(m);
        for (int i = 0; i < m; i++) a[i] = f[i];
        vector<int> b = g.a;
        a.resize(N); b.resize(N);
        ntt(a, false); ntt(b, false);
        for (int i = 0; i < N; i++)
            b[i] = b[i] * norm(2 - a[i] * b[i] % mod) % mod;
        ntt(b, true);
        b.resize(m);
        g = Poly(move(b));
    }
    g.resize(n);
    return g;
}

Poly poly_ln(const Poly &f, int n) {
    if (n <= 0) return Poly();
    assert(f.size() && f[0] == 1);
    Poly g = integ(deriv(f) * poly_inv(f, n));
    g.resize(n);
    return g;
}

Poly poly_exp(const Poly &f, int n) {
    if (n <= 0) return Poly();
    assert(!f.size() || f[0] == 0);
    Poly g({1});
    for (int m = 2; m < 2 * n; m <<= 1) {
        Poly L = poly_ln(g, m);
        vector<int> t(m);
        for (int i = 0; i < m; i++) t[i] = norm(f[i] - L[i]);
        t[0] = (t[0] + 1) % mod;
        g = g * Poly(move(t));
        g.resize(m);
    }
    g.resize(n);
    return g;
}

Poly poly_sqrt(const Poly &f, int n) {
    if (n <= 0) return Poly();
    assert(f.size() && f[0] == 1);
    int i2 = (mod + 1) / 2;
    Poly g({1});
    for (int m = 2; m < 2 * n; m <<= 1) {
        Poly t = f * poly_inv(g, m);
        t.resize(m);
        for (int i = 0; i < m; i++)
            t[i] = (t[i] + (i < g.size() ? g.a[i] : 0)) * i2 % mod;
        g = move(t);
    }
    g.resize(n);
    return g;
}

Poly poly_pow_bin(Poly f, int k, int n) {
    assert(k >= 0);
    if (n <= 0) return Poly();
    Poly g({1});
    while (k) {
        if (k & 1) {
            g = g * f;
            if (g.size() > n) g.resize(n);
        }
        k >>= 1;
        if (k) {
            f = f * f;
            if (f.size() > n) f.resize(n);
        }
    }
    g.resize(n);
    return g;
}

Poly poly_pow_core(const Poly &f, int n, int km, int kp) {
    int c = f[0], invc = inv_mod(c);
    vector<int> v(min(n, f.size()));
    for (int i = 0; i < (int)v.size(); i++) v[i] = f[i] * invc % mod;
    Poly g = poly_ln(Poly(move(v)), n);
    for (int &x : g.a) x = x * km % mod;
    g = poly_exp(g, n);
    int ck = pw(c, kp);
    for (int &x : g.a) x = x * ck % mod;
    return g;
}

Poly poly_pow(const Poly &f, int n, int k) {
    assert(k >= 0);
    if (n <= 0) return Poly();
    if (k == 0) { Poly g; g.resize(n); g[0] = 1; return g; }
    int s = 0;
    while (s < f.size() && f[s] == 0) s++;
    if (s == f.size() || (__int128)s * k >= n) return Poly(vector<int>(n));
    if (s) {
        int sk = s * k;
        Poly g = poly_pow(Poly(vector<int>(f.a.begin() + s, f.a.end())), n - sk, k);
        Poly h; h.resize(n);
        for (int i = 0; i < g.size(); i++) h[i + sk] = g[i];
        return h;
    }
    return poly_pow_core(f, n, k % mod, k % (mod - 1));
}

Poly poly_pow(const Poly &f, int n, const string &s) {
    if (n <= 0) return Poly();
    bool zero = true;
    for (char c : s) if (c != '0') zero = false;
    if (zero) { Poly g; g.resize(n); g[0] = 1; return g; }
    int km = 0, kp = 0;
    for (char c : s) {
        km = (km * 10 + c - '0') % mod;
        kp = (kp * 10 + c - '0') % (mod - 1);
    }
    int sh = 0;
    while (sh < f.size() && f[sh] == 0) sh++;
    if (sh == f.size()) return Poly(vector<int>(n));
    if (!sh) return poly_pow_core(f, n, km, kp);
    int limit = (n - 1) / sh, k = 0;
    for (char c : s) {
        int d = c - '0';
        if (k > limit / 10 || (k == limit / 10 && d > limit % 10))
            return Poly(vector<int>(n));
        k = k * 10 + d;
    }
    int sk = sh * k;
    Poly g = poly_pow_core(Poly(vector<int>(f.a.begin() + sh, f.a.end())), n - sk, km, kp);
    Poly h; h.resize(n);
    for (int i = 0; i < g.size(); i++) h[i + sk] = g[i];
    return h;
}

pair<Poly, Poly> poly_divmod(const Poly &f, const Poly &g) {
    int n = f.size(), m = g.size();
    while (n && f[n - 1] == 0) n--;
    while (m && g[m - 1] == 0) m--;
    assert(m > 0);
    if (n < m) return {Poly(), f};
    int d = n - m + 1;
    vector<int> a(f.a.begin(), f.a.begin() + n), b(g.a.begin(), g.a.begin() + m);
    reverse(a.begin(), a.end()); reverse(b.begin(), b.end());
    a.resize(d); b.resize(d);
    Poly q = Poly(move(a)) * poly_inv(Poly(move(b)), d);
    q.resize(d);
    reverse(q.a.begin(), q.a.end());
    Poly r = f - q * g;
    r.resize(m - 1);
    return {q, r};
}
Poly operator/(const Poly &f, const Poly &g) { return poly_divmod(f, g).first; }
Poly operator%(const Poly &f, const Poly &g) { return poly_divmod(f, g).second; }

Poly poly_shift(const Poly &f, int c) {
    int n = f.size();
    if (!n) return Poly();
    assert(n < mod);
    c = norm(c);
    vector<int> fac(n), ifac(n), A(n), B(n), r(n);
    fac[0] = 1;
    for (int i = 1; i < n; i++) fac[i] = fac[i - 1] * i % mod;
    ifac[n - 1] = inv_mod(fac[n - 1]);
    for (int i = n - 2; i >= 0; i--) ifac[i] = ifac[i + 1] * (i + 1) % mod;
    for (int i = 0; i < n; i++) A[n - 1 - i] = f[i] * fac[i] % mod;
    for (int i = 0, ck = 1; i < n; i++, ck = ck * c % mod)
        B[i] = ck * ifac[i] % mod;
    auto h = mul_p(A, B);
    for (int i = 0; i < n; i++) r[i] = h[n - 1 - i] * ifac[i] % mod;
    return Poly(move(r));
}

vector<int> stirling1_unsigned(int n) {
    assert(n >= 0 && n < mod);
    if (n == 0) return {1};
    if (n == 1) return {0, 1};
    int m = n / 2;
    Poly f(stirling1_unsigned(m));
    Poly g = f * poly_shift(f, m);
    if (n & 1) g = g * Poly({n - 1, 1});
    g.resize(n + 1);
    return g.a;
}
vector<int> stirling1_signed(int n) {
    auto a = stirling1_unsigned(n);
    for (int k = 0; k <= n; k++) if ((n - k) & 1) a[k] = norm(-a[k]);
    return a;
}
vector<int> stirling2(int n) {
    assert(n >= 0 && n < mod);
    vector<int> fac(n + 1), ifac(n + 1), a(n + 1), b(n + 1);
    fac[0] = 1;
    for (int i = 1; i <= n; i++) fac[i] = fac[i - 1] * i % mod;
    ifac[n] = inv_mod(fac[n]);
    for (int i = n; i >= 1; i--) ifac[i - 1] = ifac[i] * i % mod;
    for (int i = 0; i <= n; i++) {
        a[i] = (i & 1) ? norm(-ifac[i]) : ifac[i];
        b[i] = pw(i, n) * ifac[i] % mod;
    }
    auto c = mul_p(a, b);
    c.resize(n + 1);
    return c;
}

void mpe_build(vector<Poly> &tr, int u, int l, int r, const vector<int> &xs) {
    if (l == r) { tr[u] = Poly({norm(-xs[l]), 1}); return; }
    int m = (l + r) / 2;
    mpe_build(tr, u * 2, l, m, xs);
    mpe_build(tr, u * 2 + 1, m + 1, r, xs);
    tr[u] = tr[u * 2] * tr[u * 2 + 1];
}
void mpe_query(const vector<Poly> &tr, int u, int l, int r, Poly f, vector<int> &res) {
    f = f % tr[u];
    if (l == r) { res[l] = f.size() ? f[0] : 0; return; }
    int m = (l + r) / 2;
    mpe_query(tr, u * 2, l, m, f, res);
    mpe_query(tr, u * 2 + 1, m + 1, r, f, res);
}
vector<int> multi_eval(const Poly &f, const vector<int> &xs) {
    int n = xs.size();
    if (!n) return {};
    vector<Poly> tr(4 * n);
    mpe_build(tr, 1, 0, n - 1, xs);
    vector<int> res(n);
    mpe_query(tr, 1, 0, n - 1, f, res);
    return res;
}

Poly interp_solve(const vector<Poly> &tr, int u, int l, int r,
                  const vector<int> &ys, const vector<int> &d) {
    if (l == r) return Poly({norm(ys[l]) * inv_mod(d[l]) % mod});
    int m = (l + r) / 2;
    return interp_solve(tr, u * 2, l, m, ys, d) * tr[u * 2 + 1]
         + interp_solve(tr, u * 2 + 1, m + 1, r, ys, d) * tr[u * 2];
}
Poly fast_interp(const vector<int> &xs, const vector<int> &ys) {
    int n = xs.size();
    assert((int)ys.size() == n);
    if (!n) return Poly();
    vector<Poly> tr(4 * n);
    mpe_build(tr, 1, 0, n - 1, xs);
    auto d = multi_eval(deriv(tr[1]), xs);
    return interp_solve(tr, 1, 0, n - 1, ys, d);
}

int linear_recur(vector<int> a, vector<int> c, int n) {
    int k = c.size();
    assert(k > 0 && (int)a.size() == k && n >= 0);
    if (n < k) return norm(a[n]);
    vector<int> Q(k + 1); Q[0] = 1;
    for (int i = 0; i < k; i++) Q[i + 1] = norm(-c[i]);
    vector<int> P = mul_p(a, Q);
    P.resize(k);
    while (n) {
        vector<int> Qm = Q;
        for (int i = 1; i < (int)Qm.size(); i += 2) Qm[i] = norm(-Qm[i]);
        auto U = mul_p(P, Qm), V = mul_p(Q, Qm);
        for (int i = 0; i < (int)P.size(); i++) {
            int j = 2 * i + (n & 1);
            P[i] = j < (int)U.size() ? U[j] : 0;
        }
        for (int i = 0; i < (int)Q.size(); i++) {
            int j = 2 * i;
            Q[i] = j < (int)V.size() ? V[j] : 0;
        }
        n >>= 1;
    }
    return P[0] * inv_mod(Q[0]) % mod;
}

int poly_eval(const Poly &f, int x) {
    x = norm(x);
    int r = 0;
    for (int i = f.size() - 1; i >= 0; i--) r = (r * x + f[i]) % mod;
    return r;
}

int lagrange_general(const vector<int> &xs, const vector<int> &ys, int k, int p = mod) {
    int n = xs.size();
    assert((int)ys.size() == n);
    int ans = 0;
    for (int i = 0; i < n; i++) {
        int a = 1, b = 1;
        for (int j = 0; j < n; j++) if (j != i) {
            a = a * norm(k - xs[j], p) % p;
            b = b * norm(xs[i] - xs[j], p) % p;
        }
        ans = (ans + norm(ys[i], p) * a % p * inv_mod(b, p)) % p;
    }
    return ans;
}

int lagrange_continuous(const vector<int> &y, int k, int p = mod) {
    int n = y.size();
    if (!n) return 0;
    assert(n < p);
    if (1 <= k && k <= n) return norm(y[k - 1], p);
    vector<int> fac(n), ifac(n), pre(n + 1), suf(n + 2);
    fac[0] = 1;
    for (int i = 1; i < n; i++) fac[i] = fac[i - 1] * i % p;
    ifac[n - 1] = inv_mod(fac[n - 1], p);
    for (int i = n - 2; i >= 0; i--) ifac[i] = ifac[i + 1] * (i + 1) % p;
    pre[0] = suf[n + 1] = 1;
    for (int i = 1; i <= n; i++) pre[i] = pre[i - 1] * norm(k - i, p) % p;
    for (int i = n; i >= 1; i--) suf[i] = suf[i + 1] * norm(k - i, p) % p;
    int ans = 0;
    for (int i = 1; i <= n; i++) {
        int t = norm(y[i - 1], p) * pre[i - 1] % p * suf[i + 1] % p
              * ifac[i - 1] % p * ifac[n - i] % p;
        if ((n - i) & 1) t = norm(-t, p);
        ans = (ans + t) % p;
    }
    return ans;
}

void fwt_or(vector<int> &a, bool inv) {
    int n = a.size();
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int &x = a[i + j], &y = a[i + j + m];
                y = inv ? norm(y - x) : (y + x) % mod;
            }
}
void fwt_and(vector<int> &a, bool inv) {
    int n = a.size();
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int &x = a[i + j], &y = a[i + j + m];
                x = inv ? norm(x - y) : (x + y) % mod;
            }
}
void fwt_xor(vector<int> &a, bool inv) {
    int n = a.size(), i2 = (mod + 1) / 2;
    for (int m = 1; m < n; m <<= 1)
        for (int i = 0; i < n; i += m << 1)
            for (int j = 0; j < m; j++) {
                int x = a[i + j], y = a[i + j + m];
                a[i + j] = (x + y) % mod;
                a[i + j + m] = norm(x - y);
                if (inv) {
                    a[i + j] = a[i + j] * i2 % mod;
                    a[i + j + m] = a[i + j + m] * i2 % mod;
                }
            }
}
template <class F>
vector<int> bit_conv(vector<int> a, vector<int> b, F fwt) {
    assert(a.size() == b.size() && a.size() && (a.size() & (a.size() - 1)) == 0);
    fwt(a, false); fwt(b, false);
    for (int i = 0; i < (int)a.size(); i++) a[i] = a[i] * b[i] % mod;
    fwt(a, true);
    return a;
}

vector<int> subset_conv(const vector<int> &f, const vector<int> &g) {
    int sz = f.size();
    assert(sz && (sz & (sz - 1)) == 0 && g.size() == f.size());
    int n = __builtin_ctzll(sz);
    vector<vector<int>> F(n + 1, vector<int>(sz)), H0(n + 1, vector<int>(sz));
    for (int s = 0; s < sz; s++) {
        int p = __builtin_popcountll(s);
        F[p][s] = norm(f[s]); H0[p][s] = norm(g[s]);
    }
    for (int i = 0; i <= n; i++) { fwt_or(F[i], false); fwt_or(H0[i], false); }
    vector<vector<int>> H(n + 1, vector<int>(sz));
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= i; j++)
            for (int s = 0; s < sz; s++)
                H[i][s] = (H[i][s] + F[j][s] * H0[i - j][s]) % mod;
    for (int i = 0; i <= n; i++) fwt_or(H[i], true);
    vector<int> res(sz);
    for (int s = 0; s < sz; s++) res[s] = H[__builtin_popcountll(s)][s];
    return res;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    return 0;
}

```
