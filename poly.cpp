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

public:
    vector<int> a;

    Poly() {}
    Poly(const vector<int>& _a) : a(_a) {}
    Poly(int val) { a.push_back(val); }

    int size() const { return a.size(); }
    void resize(int n) { a.resize(n, 0); }
    int operator[](int idx) const { return idx < size() ? a[idx] : 0; }

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
        
        A.resize(n + m - 1); // 去除多余的 0
        return Poly(A);
    }

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
            B.resize(len); // 截断到当前精度
            res = Poly(B);
        }
        res.resize(deg);
        return res;
    }

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
            tmp[0] = (tmp[0] + 1) % mod; // (1 - ln(G) + F)

            res = res * Poly(tmp);
            res.resize(len); // 截断
        }
        res.resize(deg);
        return res;
    }

    // 单点求值 (Horner)
    int eval(int x) const {
        int res = 0;
        for (int i = size() - 1; i >= 0; i--) {
            res = (1ll * res * x + a[i]) % mod;
        }
        return res;
    }

    // 多项式开根，要求 a[0] == 1
    Poly sqrt(int deg) const {
        if (deg <= 0) return Poly();
        Poly res(vector<int>{1}); // sqrt(1) = 1
        for (int len = 2; len < (deg << 1); len <<= 1) {
            vector<int> f_vec(len, 0);
            for (int i = 0; i < len; i++) f_vec[i] = (*this)[i];
            Poly F(f_vec);

            Poly inv_res = res.inv(len);
            // tmp = F * inv(res)  =>  res = (res + tmp) / 2
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

    // 多项式快速幂 F^k (k 较小时)
    Poly pow(int deg, long long k) const {
        if (deg <= 0) return Poly();
        if (k == 0) {
            Poly res; res.a.resize(deg, 0); res.a[0] = 1;
            return res;
        }

        // 找到第一个非零系数
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

        // a[0] != 0: F^k = a[0]^k * exp(k * ln(F / a[0]))
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

    // 连续点值插值: y[0..n-1] = f(1), f(2), ..., f(n), 求 f(k)
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
            if ((n - i) & 1) term = mod - term;  // (-1)^{n-i}
            ans = (ans + term) % mod;
        }
        return ans;
    }
};
