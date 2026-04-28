#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Poly {
private:
    static const int mod = 998244353;
    static const int G = 3;
    
    // 快速幂
    static int qpow(int a, int k) {
        int res = 1;
        for (; k; a = 1ll * a * a % mod, k >>= 1) {
            if (k & 1) res = 1ll * res * a % mod;
        }
        return res;
    }

    // 动态计算蝴蝶变换的 rev 数组
    static vector<int> get_rev(int limit) {
        vector<int> rev(limit);
        int bit = 0;
        while ((1 << bit) < limit) bit++;
        for (int i = 0; i < limit; i++) {
            rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
        }
        return rev;
    }

    // 底层 NTT，就地变换
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

    // 构造函数
    Poly() {}
    Poly(const vector<int>& _a) : a(_a) {}
    Poly(int val) { a.push_back(val); }

    int size() const { return a.size(); }
    void resize(int n) { a.resize(n, 0); }
    int operator[](int idx) const { return idx < size() ? a[idx] : 0; }

    // 多项式乘法
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

    // 求导
    Poly deriv() const {
        if (size() <= 1) return Poly(vector<int>{0});
        vector<int> res(size() - 1);
        for (int i = 1; i < size(); i++) {
            res[i - 1] = 1ll * a[i] * i % mod;
        }
        return Poly(res);
    }

    // 积分
    Poly integ() const {
        if (size() == 0) return Poly();
        vector<int> res(size() + 1, 0);
        for (int i = 0; i < size(); i++) {
            res[i + 1] = 1ll * a[i] * qpow(i + 1, mod - 2) % mod;
        }
        return Poly(res);
    }

    // 多项式求逆 (mod x^deg)
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

    // 多项式 Ln (mod x^deg)
    Poly ln(int deg) const {
        // 要求 a[0] == 1
        Poly res = (this->deriv() * this->inv(deg)).integ();
        res.resize(deg);
        return res;
    }

    // 多项式 Exp (mod x^deg)
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
};
