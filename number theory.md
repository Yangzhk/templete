- [类欧几里得](#类欧几里得)
- [莫比乌斯函数](#莫比乌斯函数)
- [数论分块](#数论分块)
- [杜教筛](#杜教筛)
- [中国剩余定理](#中国剩余定理)
- [Lucas 定理](#lucas-定理)
- [BSGS（离散对数）](#bsgs离散对数)
- [二次剩余](#二次剩余)
- [Miller-Rabin \& Pollard Rho](#miller-rabin--pollard-rho)
### 类欧几里得
```
//(ai + b) / c, 0 <= i <= n
#define int long long

int floor_sum(int n, int a, int b, int c) {
    n++;
    int ans = 0;
    while (true) {
        if (a >= c) {
            ans += n * (n - 1) / 2 * (a / c);
            a %= c;
        }
        if (b >= c) {
            ans += n * (b / c);
            b %= c;
        }
        int y = a * n + b;
        if (y < c) break;
        n = y / c;
        b = y % c;
        swap(a, c);
    }
    return ans;
}
```

```
//f : 1, g : i, h : ()^2

#define int long long

struct Node {
    int f, g, h;
};

Node euclid(int n, int a, int b, int c) {
    int qa = a / c, qb = b / c;
    a %= c, b %= c;

    int s1 = n * (n + 1) / 2;
    int s2 = n * (n + 1) * (2 * n + 1) / 6;

    int f = 0, g = 0, h = 0;
    int m = (a * n + b) / c;

    if (m) {
        auto [F, G, H] = euclid(m - 1, c, c - b - 1, a);
        f = n * m - F;
        g = m * s1 - (H + F) / 2;
        h = n * m * m - 2 * G - F;
    }

    return {
        f + qa * s1 + qb * (n + 1),
        g + qa * s2 + qb * s1,
        h + 2 * qa * g + 2 * qb * f
          + qa * qa * s2 + 2 * qa * qb * s1
          + qb * qb * (n + 1)
    };
}

```

```
//带模数

#define int long long
const int mod = 998244353;
const int inv2 = 499122177, inv6 = 166374059;

struct Z {
    int x;
    Z(int v = 0) : x((v % mod + mod) % mod) {}
    Z operator + (Z b) const { return x + b.x; }
    Z operator - (Z b) const { return x - b.x; }
    Z operator * (Z b) const { return x * b.x; }
};

struct Node {
    Z f, g, h;
};

Node euclid(int n, int a, int b, int c) {
    Z A = a / c, B = b / c;
    a %= c, b %= c;

    Z N = n;
    Z s1 = N * (N + 1) * inv2;
    Z s2 = N * (N + 1) * (N * 2 + 1) * inv6;

    Z f = 0, g = 0, h = 0;
    int m = ((__int128)a * n + b) / c;

    if (m) {
        auto [F, G, H] = euclid(m - 1, c, c - b - 1, a);
        f = N * m - F;
        g = s1 * m - (H + F) * inv2;
        h = N * m * m - G * 2 - F;
    }

    return {
        f + A * s1 + B * (N + 1),
        g + A * s2 + B * s1,
        h + A * g * 2 + B * f * 2
          + A * A * s2 + A * B * s1 * 2 + B * B * (N + 1)
    };
}


```

### 莫比乌斯函数

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

### 数论分块

```
// sum_{i=1..n} floor(n/i), O(sqrt(n))
long long sumDiv(long long n) {
    long long ans = 0;
    for (long long l = 1, r; l <= n; l = r + 1) {
        r = n / (n / l);
        ans += (r - l + 1) * (n / l);
    }
    return ans;
}

// 二维分块: floor(n/i) 和 floor(m/i) 同时分块
// sf[i] = f(1) + f(2) + ... + f(i) (例如 mu 的前缀和)
long long sumDiv2(long long n, long long m, vector<long long>& sf) {
    long long ans = 0, lim = min(n, m);
    for (long long l = 1, r; l <= lim; l = r + 1) {
        r = min(n / (n / l), m / (m / l));
        ans += (sf[r] - sf[l - 1]) * (n / l) * (m / l);
    }
    return ans;
}
```

### 杜教筛

设 $F(n)=\sum_{i=1}^{n}f(i)$。利用狄利克雷卷积 $f*1=g$，可得

$$
F(n)=\sum_{i=1}^{n}g(i)-\sum_{i=2}^{n}F\left(\left\lfloor\frac ni\right\rfloor\right).
$$

后面的求和使用**整除分块 + 记忆化搜索**。预处理到 $B=\lfloor n^{2/3}\rfloor$；以下两份模板**独立使用**。

#### 1. 莫比乌斯函数前缀和

由 $\mu*1=\varepsilon$，其中 $\varepsilon(1)=1$、其余为 $0$，得到

$$
\boxed{M(n)=\sum_{i=1}^{n}\mu(i)=1-\sum_{i=2}^{n}M\left(\left\lfloor\frac ni\right\rfloor\right)}.
$$

```cpp
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define rep(i, a, b) for(int i = a; i <= b; i++)

int B;
vector<int> mu, p, vis, s;
map<int, int> mp;

void init(int n) {
    B = n;
    p.clear(); mp.clear();
    mu.assign(n + 1, 0);
    vis.assign(n + 1, 0);
    s.assign(n + 1, 0);
    mu[1] = 1;

    rep(i, 2, n) {
        if(!vis[i]) {
            p.push_back(i);
            mu[i] = -1;
        }
        for(auto j : p) {
            if(i * j > n) break;
            vis[i * j] = 1;
            if(i % j == 0) {
                mu[i * j] = 0;
                break;
            }
            mu[i * j] = -mu[i];
        }
    }
    rep(i, 1, n) s[i] = s[i - 1] + mu[i];
}

int get(int n) {
    if(n <= B) return s[n];
    if(mp.count(n)) return mp[n];

    int ans = 1;
    for(int l = 2, r; l <= n; l = r + 1) {
        r = n / (n / l);
        ans -= (r - l + 1) * get(n / l);
    }
    return mp[n] = ans;
}
```

#### 2. 欧拉函数前缀和

由 $\varphi*1=\operatorname{id}$（即 $\sum_{d\mid n}\varphi(d)=n$），得到

$$
\boxed{\Phi(n)=\sum_{i=1}^{n}\varphi(i)=\frac{n(n+1)}2-\sum_{i=2}^{n}\Phi\left(\left\lfloor\frac ni\right\rfloor\right)}.
$$

```cpp
#include <bits/stdc++.h>
using namespace std;
#define int long long
#define rep(i, a, b) for(int i = a; i <= b; i++)

int B;
vector<int> phi, p, vis, s;
map<int, int> mp;

void init(int n) {
    B = n;
    p.clear(); mp.clear();
    phi.assign(n + 1, 0);
    vis.assign(n + 1, 0);
    s.assign(n + 1, 0);
    phi[1] = 1;

    rep(i, 2, n) {
        if(!vis[i]) {
            p.push_back(i);
            phi[i] = i - 1;
        }
        for(auto j : p) {
            if(i * j > n) break;
            vis[i * j] = 1;
            if(i % j == 0) {
                phi[i * j] = phi[i] * j;
                break;
            }
            phi[i * j] = phi[i] * (j - 1);
        }
    }
    rep(i, 1, n) s[i] = s[i - 1] + phi[i];
}

int get(int n) {
    if(n <= B) return s[n];
    if(mp.count(n)) return mp[n];

    int ans = n * (n + 1) / 2;
    for(int l = 2, r; l <= n; l = r + 1) {
        r = n / (n / l);
        ans -= (r - l + 1) * get(n / l);
    }
    return mp[n] = ans;
}
```

#### 调用与复杂度

```cpp
init(max(1LL, (int)pow((long double)n, 2.0L / 3)));
cout << get(n) << '\n';
```

- 整除分块右端点：$r=\left\lfloor n/\lfloor n/l\rfloor\right\rfloor$。
- `map` 版本：时间 $O(n^{2/3}\log n)$，空间 $O(n^{2/3})$。
- 不取模；$\varphi$ 前缀和注意 `long long` 溢出（例如建议 $n\le 10^9$）。


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

### Lucas 定理

对**素数** $p$：

$$
\binom{n}{m} \equiv \binom{n \bmod p}{m \bmod p} \cdot \binom{\lfloor n/p \rfloor}{\lfloor m/p \rfloor} \pmod p
$$

适用于 $n, m$ 极大但 $p$ 较小（$p \le 10^5$ 量级）的二项式系数取模。时间复杂度 $O(p + \log_p n)$。

```
long long lucas_pow(long long a, long long b, long long p) {
    long long r = 1; a %= p;
    for (; b; b >>= 1, a = a * a % p) if (b & 1) r = r * a % p;
    return r;
}

// 直接组合: 要求 0 <= m, n < p, p 为素数
long long C_small(long long n, long long m, long long p) {
    if (m > n || m < 0) return 0;
    long long a = 1, b = 1;
    for (long long i = 0; i < m; i++) {
        a = a * ((n - i) % p) % p;
        b = b * (i + 1) % p;
    }
    return a * lucas_pow(b, p - 2, p) % p;
}

long long Lucas(long long n, long long m, long long p) {
    if (m == 0) return 1;
    return C_small(n % p, m % p, p) * Lucas(n / p, m / p, p) % p;
}
```

### BSGS（离散对数）

求解 $a^x \equiv b \pmod p$，先考虑 $\gcd(a, p) = 1$ 的情形。

设 $m = \lceil \sqrt p \rceil$，写 $x = im - j$，$0 \le j < m,\ 1 \le i \le m$，则
$$
a^{im} \equiv b \cdot a^j \pmod p
$$

枚举 $j$ 把 $b a^j$ 入哈希表，再枚举 $i$ 查表。时间复杂度 $O(\sqrt p)$。

```
long long bsgs_pow(long long a, long long b, long long p) {
    long long r = 1; a %= p;
    for (; b; b >>= 1, a = a * a % p) if (b & 1) r = r * a % p;
    return r;
}

// a^x ≡ b (mod p), gcd(a, p) = 1; 无解返回 -1
long long BSGS(long long a, long long b, long long p) {
    a %= p; b %= p;
    if (b == 1 || p == 1) return 0;
    long long m = (long long)ceil(sqrt((double)p));
    unordered_map<long long, long long> mp;
    long long cur = b;
    for (long long j = 0; j < m; j++) {
        mp[cur] = j;
        cur = cur * a % p;
    }
    long long am = bsgs_pow(a, m, p);
    cur = 1;
    for (long long i = 1; i <= m; i++) {
        cur = cur * am % p;
        if (mp.count(cur)) return i * m - mp[cur];
    }
    return -1;
}

// 不要求 gcd(a, p) = 1: 不断把公因子从两边消掉, 直到互素再调用 BSGS
long long exBSGS(long long a, long long b, long long p) {
    a = (a % p + p) % p; b = (b % p + p) % p;
    if (b == 1 || p == 1) return 0;
    long long g, k = 0, ak = 1;
    while ((g = __gcd(a, p)) > 1) {
        if (b % g) return -1;
        b /= g; p /= g;
        ak = ak * (a / g) % p;
        k++;
        if (ak == b) return k;
    }
    // 此时 gcd(ak, p) = 1, 用 exgcd 求 ak^{-1} (mod p)
    long long x = 1, y = 0, x1 = 0, y1 = 1, A = ak, P = p;
    while (P) {
        long long q = A / P;
        tie(x, x1) = make_pair(x1, x - q * x1);
        tie(y, y1) = make_pair(y1, y - q * y1);
        tie(A, P) = make_pair(P, A - q * P);
    }
    long long inv = (x % p + p) % p;
    long long r = BSGS(a, b * inv % p, p);
    return r < 0 ? -1 : r + k;
}
```

### 二次剩余

```
#include <iostream>
#include <random>

long long p, v;  // 分别是模数和 r^2 - a 的值

struct Poly {
  long long a, b;

  Poly(long long _a = 0, long long _b = 0) : a(_a), b(_b) {}
};

Poly operator*(const Poly& x, const Poly& y) {
  // 重载乘法，可以参考上面有关运算性质的说明
  return Poly((x.a * y.a + v * (x.b * y.b % p)) % p,
              (x.a * y.b + x.b * y.a) % p);
}

// 多项式的快速幂，用于计算答案
Poly modpow(Poly a, long long b) {
  Poly res(1, 0);
  while (b) {
    if (b & 1) res = res * a;
    a = a * a;
    b >>= 1;
  }
  return res;
}

// 普通的快速幂，用于判断二次非剩余
long long modpow(long long a, long long b) {
  long long res = 1;
  while (b) {
    if (b & 1) res = res * a % p;
    a = a * a % p;
    b >>= 1;
  }
  return res;
}

// 用于生成随机数
std::mt19937 rng(std::random_device{}());

long long cipolla(long long a, long long _p) {
  p = _p;
  if (a == 0)
    return 0;  // 特判一下 0 的情况
  else if (modpow(a, (p - 1) / 2) == p - 1)
    return -1;  // 判断二次非剩余，此时无解
  else {
    // 随机 r，使得 r^2 - a 是一个二次非剩余
    long long r;
    for (r = rng() % p;; r = rng() % p) {
      if (modpow((r * r - a + p) % p, (p - 1) / 2) == p - 1) break;
    }
    v = (r * r - a + p) % p;
    return modpow(Poly(r, 1), (p + 1) / 2).a;  // 根据结论式计算结果
  }
}

int main() {
  int t, a, p;
  std::cin >> t;
  while (t--) {
    std::cin >> a >> p;
    int ans = cipolla(a, p);
    if (ans == -1)
      std::cout << "Hola!" << std::endl;
    else if (ans == 0)
      std::cout << 0 << std::endl;
    else {
      // 相反数是另一个解
      int ans2 = (p - ans) % p;
      if (ans2 < ans) std::swap(ans, ans2);
      std::cout << ans << " " << ans2 << std::endl;
    }
  }
  return 0;
}
```

### Miller-Rabin & Pollard Rho

**Miller-Rabin**：写 $n - 1 = 2^s \cdot d$（$d$ 奇）。若 $n$ 为奇素数，则对任意 $\gcd(a, n) = 1$，要么 $a^d \equiv 1$，要么存在 $0 \le r < s$ 使 $a^{2^r d} \equiv -1 \pmod n$。取见证集 $\{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37\}$ 即可对 $n < 2^{64}$ 给出**确定性**判定。

**Pollard Rho**：迭代 $x_{k+1} = x_k^2 + c \pmod n$ 产生 $\rho$ 形循环，期望 $O(n^{1/4})$ 步通过 $\gcd(|x - y|, n)$ 找到非平凡因子。

```
using u64 = unsigned long long;
using u128 = __uint128_t;

u64 mulmod(u64 a, u64 b, u64 m) { return (u128)a * b % m; }

u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    for (; e; e >>= 1, a = mulmod(a, a, m)) if (e & 1) r = mulmod(r, a, m);
    return r;
}

bool millerRabin(u64 n) {
    if (n < 2) return false;
    for (u64 p : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) d >>= 1, s++;
    for (u64 a : {2ULL, 3ULL, 5ULL, 7ULL, 11ULL, 13ULL, 17ULL, 19ULL, 23ULL, 29ULL, 31ULL, 37ULL}) {
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 0; r < s - 1; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = false; break; }
        }
        if (comp) return false;
    }
    return true;
}

u64 pollardRho(u64 n) {
    if (n % 2 == 0) return 2;
    static mt19937_64 rng(20231103);
    while (true) {
        u64 c = rng() % (n - 1) + 1;
        auto f = [&](u64 x) { return (mulmod(x, x, n) + c) % n; };
        u64 x = rng() % n, y = x, d = 1;
        do {
            x = f(x);
            y = f(f(y));
            u64 diff = x > y ? x - y : y - x;
            d = __gcd(diff, n);
        } while (d == 1);
        if (d != n) return d;
    }
}

// 分解 n 的所有素因子（含重复）, 放入 res
void factorize(u64 n, vector<u64>& res) {
    if (n == 1) return;
    if (millerRabin(n)) { res.push_back(n); return; }
    u64 d = n;
    while (d == n) d = pollardRho(n);
    factorize(d, res);
    factorize(n / d, res);
}

// n 的最大质因子 (典型应用)
u64 maxPrimeFactor(u64 n) {
    vector<u64> p; factorize(n, p);
    return *max_element(p.begin(), p.end());
}
```
