
## 目录

- [模运算](#模运算)
- [矩阵运算类](#矩阵运算类)
- [莫队](#莫队)
- [xor 线性基](#xor-线性基)
- [mobius 函数 (反演)](#mobius-函数-反演)
- [数论分块（整除分块）](#数论分块整除分块)
- [exgcd 求 线性同余方程](#exgcd-求-线性同余方程)
- [中国剩余定理](#中国剩余定理)
- [Lucas 定理](#lucas-定理)
- [BSGS（离散对数）](#bsgs离散对数)
- [Miller-Rabin & Pollard Rho](#miller-rabin--pollard-rho)
- [最大流](#最大流)
- [zkw 费用流](#zkw-费用流)
- [计算几何 (dls version)](#计算几何-dls-version)
- [三维点类](#三维点类)

---

### 模运算

```
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


##### 5. 图论问题

- 在图上维护路径异或（如树上路径 xor、带环图）
- 求：
  - 两点路径 xor 最大值
  - 环空间的 xor 基

---

##### 6. 区间 / 前缀问题

- 前缀异或 + 线性基：
  - 区间最大 xor
  - 子数组 xor 相关问题

```
// 线性基 (60位, 适用于 1e18)
struct XORBasis {
    long long p[60] = {};

    bool insert(long long x) {
        for (int i = 60; ~i; i--) {
            if (!(x >> i)) continue;
            if (!p[i]) { p[i] = x; return true; }
            x ^= p[i];
        }
        return false;
    }

    long long maxXor(long long init = 0) {
        long long res = init;
        for (int i = 60; ~i; i--)
            if ((res ^ p[i]) > res) res ^= p[i];
        return res;
    }

    void rebuild() { // 行最简形 (用于 kth)
        for (int i = 60; ~i; i--) {
            if (!p[i]) continue;
            for (int j = i-1; ~j; j--)
                if (p[i] >> j & 1) p[i] ^= p[j];
        }
    }

    long long kth(long long k) { // 需先 rebuild, k 从 1 开始
        long long res = 0;
        for (int i = 60; ~i; i--) {
            if (!p[i]) continue;
            if (k & 1) res ^= p[i];
            k >>= 1;
        }
        return res;
    }
};

// 图论异或路径: 任意路径 xor = dis[u] ^ dis[v] ^ (环的任意组合)
// dfs 找生成树, 环插入线性基
XORBasis cycleBase;
long long dis[N];
bool vis[N];

void dfs(int u, int fa) {
    vis[u] = true;
    for (auto [v, w] : g[u]) {
        if (v == fa) continue;
        if (!vis[v]) {
            dis[v] = dis[u] ^ w;
            dfs(v, u);
        } else {
            cycleBase.insert(dis[u] ^ dis[v] ^ w);
        }
    }
}
// 查询: cycleBase.maxXor(dis[u] ^ dis[v])

// 前缀线性基: 每个前缀维护一个基, insert 时尽量保留靠右的元素
struct PrefixBasis {
    static const int LOG = 60;
    array<long long, LOG> b{};
    array<int, LOG> pos{};

    PrefixBasis() { pos.fill(-1); }

    void insert(long long x, int id) {
        for (int j = LOG-1; ~j; j--) {
            if (!(x >> j)) continue;
            if (!b[j]) { b[j] = x; pos[j] = id; break; }
            if (id > pos[j]) { swap(x, b[j]); swap(id, pos[j]); }
            x ^= b[j];
        }
    }

    long long maxXor(int l, long long x = 0) { // 查询 [l, n-1]
        long long res = x;
        for (int j = LOG-1; ~j; j--)
            if (b[j] && pos[j] >= l && (res ^ b[j]) > res)
                res ^= b[j];
        return res;
    }
};
// 使用: for (int i = 0; i < n; i++) pref.insert(a[i], i);
//       查询 a[l..r] 最大 xor: 对每个 r 维护一个版本, 查询 pref[r].maxXor(l)
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

---

#### 关键恒等式

$$
\sum_{d \mid n} \mu(d) = [n = 1]
$$

这是反演公式成立的根本原因（即 $\mu * 1 = \varepsilon$，狄利克雷卷积下 $\mu$ 是常函数 $1$ 的逆元，$\varepsilon(n) = [n=1]$ 是单位元）。

---

#### 互素判别（最常用技巧）

$$
[\gcd(i, j) = 1] = \sum_{d \mid \gcd(i, j)} \mu(d)
$$

可把含 $\gcd$ 的求和限制改写为对约数的求和，进而交换求和顺序、提取出 $\lfloor n/d \rfloor$ 形式。

---

#### 典型应用：互素对计数

求 $\displaystyle\sum_{i=1}^{n} \sum_{j=1}^{m} [\gcd(i,j) = 1]$：

$$
\sum_{i=1}^{n}\sum_{j=1}^{m} \sum_{d \mid \gcd(i,j)} \mu(d)
= \sum_{d=1}^{\min(n,m)} \mu(d) \left\lfloor \frac{n}{d} \right\rfloor \left\lfloor \frac{m}{d} \right\rfloor
$$

配合**数论分块**可在 $O(\sqrt{n})$ 时间求值。若题目要求 $\gcd(i,j) = k$，先令 $i = ki', j = kj'$ 化为 $\gcd = 1$。

---

#### 与积性函数的卷积

记 $\text{id}(n) = n, \mathbf{1}(n) = 1$，常用关系：

- $\varphi * \mathbf{1} = \text{id}$，即 $\sum_{d \mid n} \varphi(d) = n$
- $\varphi = \mu * \text{id}$，即 $\varphi(n) = \sum_{d \mid n} \mu(d) \cdot \dfrac{n}{d}$
- $d(n) = \mathbf{1} * \mathbf{1}$，$\sigma(n) = \mathbf{1} * \text{id}$

---

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

### 数论分块（整除分块）

利用 $\lfloor n/d \rfloor$ 至多有 $O(\sqrt n)$ 个不同取值。固定值 $v = \lfloor n/d \rfloor$ 时，最大的 $d$ 等于 $\lfloor n / v \rfloor$，因此可以把这些相同取值的 $d$ 合并为一个区间 $[l, r]$ 一起处理。

适用场景：与莫比乌斯反演组合（$\sum \mu(d) \lfloor n/d \rfloor \lfloor m/d \rfloor$）、约数和、欧拉函数前缀和、$\sum_{i=1}^n \lfloor n/i \rfloor$ 等。

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

### 计算几何 (dls version)

```
#include <bits/stdc++.h>

#define rep(i, a, b) for(int i = a; i <= b; i++)
#define db double

using namespace std;

const db EPS = 1e-9;

inline int sign(db a) { return a < -EPS ? -1 : a > EPS; }

inline int cmp(db a, db b) { return sign(a - b); }

struct P {
    db x, y;

    P() {}
    P(db _x, db _y) : x(_x), y(_y) {}
    
    P operator + (P p) { return {x + p.x, y + p.y}; }
    P operator - (P p) { return {x - p.x, y - p.y}; }
    P operator * (db d) { return {x * d, y * d}; }
    P operator / (db d) { return {x / d, y / d}; }

    bool operator < (P p) const {
        int c = cmp(x, p.x);
        if(c) return c == -1;
        return cmp(y, p.y) == -1;
    }

    /*极角排序
    bool operator < (P p) const {
        if(quad() != b.quad()) return quad() < p.quad();
        return sign(det(p)) == -1;
    } 
    */

    bool operator == (P o) const {
        return cmp(x, o.x) == 0 && cmp(y, o.y) == 0;
    }

    db distTo(P p) { return (*this - p).abs(); }
    db alpha() { return atan2(y, x); }
    void read() { cin >> x >> y; }
    void write() { cout << "(" << x << "," << y << ")" << "\n"; }
    db abs() { return sqrt(abs2()); }
    db abs2() { return x * x + y * y; }
    P rot90() { return P(-y, x); }
    P unit() { return *this / abs(); }
    int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }

    db dot(P p) { return x * p.x + y * p.y; }
    db det(P p) { return x * p.y - y * p.x; }

    P rot(db an) { return {x * cos(an) - y * sin(an), x * sin(an) + y * cos(an)}; }

};

db cross(P p1, P p2, P p3) { return ((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)); }
db crossOp(P p1, P p2, P p3) { return sign(cross(p1, p2, p3)); }

bool chkLL(P p1, P p2, P q1, P q2) {//直线平行
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return sign(a1 + a2) != 0;
}

P isLL(P p1, P p2, P q1, P q2) { // 直线交点
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return (p1 * a2 + p2 * a1) / (a1 + a2);
}

bool intersect(db l1, db r1, db l2, db r2) { // 跨立实验
    if(l1 > r1) swap(l1, r1); if(l2 > r2) swap(l2, r2);
    return !( cmp(r1, l2) == -1 || cmp(r2, l1) == -1 );
}

bool isSS(P p1, P p2, P q1, P q2) { // 线段相交
    return intersect(p1.x, p2.x, q1.x, q2.x) && intersect(p1.y, p2.y, q1.y, q2.y) &&
    crossOp(p1, p2, q1) * crossOp(p1, p2, q2) <= 0 && crossOp(q1, q2, p1)
            * crossOp(q1, q2, p2) <= 0;
}

bool isSS_strict(P p1, P p2, P q1, P q2) { // 是否严格相交
    return crossOp(p1, p2, q1) * crossOp(p1, p2, q2) < 0 && crossOp(q1, q2, p1)
            * crossOp(q1, q2, p2) < 0;
}

bool isMiddle(db a, db m, db b) {
    return sign(a - m) == 0 || sign(b - m) == 0 || (a < m != b < m);
}

bool isMiddle(P a, P m, P b) {
    return isMiddle(a.x, m.x, b.x) && isMiddle(a.y, m.y, b.y);
}

bool onSeg(P p1, P p2, P q) { // 点在线段上
    return crossOp(p1, p2, q) == 0 && isMiddle(p1, q, p2);
}

bool onSeg_strict(P p1, P p2, P q) {
    return crossOp(p1, p2, q) == 0 && sign((q - p1).dot(p1 - p2) * sign((q - p2).dot(p1 - p2)) < 0);
}

P proj(P p1, P p2, P q) { // 投影
    P dir = p2 - p1;
    return p1 + dir * (dir.dot(q - p1) / dir.abs2());
}

P reflect(P p1, P p2, P q) { // 对称点 
    return proj(p1, p2, q) * 2 - q;
}

db nearest(P p1, P p2, P q) { // 点到线段的最短距离
    P h = proj(p1, p2, q);
    if(isMiddle(p1, h, p2))
        return q.distTo(h);
    return min(p1.distTo(q), p2.distTo(q));
}

db disSS(P p1, P p2, P q1, P q2) { // 线段到线段的最短距离
    if(isSS(p1, p2, q1, q2)) return 0;
    return min(min(nearest(p1, p2, q1), nearest(p1, p2, q2)), min(nearest(q1, q2, p1), nearest(q1, q2, p2)));
}

db rad(P p1, P p2) { // 向量夹角
    return atan2l(p1.det(p2), p1.dot(p2));
}

db area(vector<P> ps) { // 多边形面积
    db ret = 0; rep(i, 0, ps.size() - 1) ret += ps[i].det(ps[(i + 1) % ps.size()]);

    return ret / 2;
}

int contain(vector<P> ps, P p) { // 2: inside, 1: on_seg, 0: outside
    int n = ps.size(), ret = 0;
    rep(i, 0, n-1) {
        P u = ps[i], v = ps[(i + 1) % n];
        if(onSeg(u, v, p)) return 1;
        if(cmp(u.y, v.y) <= 0) swap(u, v);
        if(cmp(p.y, u.y) > 0 || cmp(p.y, v.y) <= 0) continue;
        ret ^= crossOp(p, u, v) > 0;
    }
    return ret * 2;
}

vector<P> convexHull(vector<P> ps) {
    int n = ps.size(); if(n <= 1) return ps;
    sort(ps.begin(), ps.end());
    vector<P> qs(n * 2); int k = 0;
    for(int i = 0; i < n; qs[k++] = ps[i++])
        while(k > 1 && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0) k--;
    for(int i = n - 2, t = k; i >= 0; qs[k++] = ps[i--])
        while(k > t && crossOp(qs[k - 2], qs[k - 1], ps[i]) <= 0) k--;
    qs.resize(k - 1);
    return qs;
}

db dist2(P a, P b) { return (a - b).abs2(); }

db convexDiameter(vector<P> ps) { // 旋转卡尺
    int n = ps.size();
    if (n <= 1) return 0;
    if (n == 2) return ps[0].distTo(ps[1]);
    
    db maxd2 = 0;
    // 旋转卡尺寻找最远点对
    for (int i = 0, j = 1; i < n; i++) {
        // 当 (i, i+1, j+1) 的面积大于 (i, i+1, j) 时，j 向后移
        // 这里使用 det (叉积) 计算平行四边形面积
        while (sign((ps[(i + 1) % n] - ps[i]).det(ps[(j + 1) % n] - ps[i]) - 
                   (ps[(i + 1) % n] - ps[i]).det(ps[j] - ps[i])) > 0) {
            j = (j + 1) % n;
        }
        maxd2 = max({maxd2, dist2(ps[i], ps[j]), dist2(ps[(i + 1) % n], ps[j])});
    }
    return sqrt(maxd2);
}

// 辅助函数：将凸包的顶点调整为从左下角开始，且为逆时针
void reorder(vector<P>& ps) {
    int pos = 0;
    for(int i = 1; i < ps.size(); i++) {
        if(ps[i].y < ps[pos].y || (ps[i].y == ps[pos].y && ps[i].x < ps[pos].x)) pos = i;
    }
    rotate(ps.begin(), ps.begin() + pos, ps.end());
}

vector<P> minkowskiSum(vector<P> A, vector<P> B) { // 闵可夫斯基和
    reorder(A); reorder(B);
    int n = A.size(), m = B.size();
    
    // 构造差分向量
    vector<P> v1(n), v2(m);
    for(int i = 0; i < n; i++) v1[i] = A[(i + 1) % n] - A[i];
    for(int i = 0; i < m; i++) v2[i] = B[(i + 1) % m] - B[i];
    
    vector<P> res;
    res.push_back(A[0] + B[0]);
    
    // 归并排序边向量（按极角）
    int i = 0, j = 0;
    while(i < n && j < m) {
        if(sign(v1[i].det(v2[j])) > 0) res.push_back(res.back() + v1[i++]);
        else if(sign(v1[i].det(v2[j])) < 0) res.push_back(res.back() + v2[j++]);
        else res.push_back(res.back() + v1[i++] + v2[j++]); // 平行向量合并
    }
    while(i < n) res.push_back(res.back() + v1[i++]);
    while(j < m) res.push_back(res.back() + v2[j++]);
    
    // 重新跑一遍凸包以去重或处理共线（可选，视具体题目精度要求而定）
    return convexHull(res); 
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
