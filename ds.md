# 数据结构模板

## 目录

- [并查集](#并查集)
  - [路径压缩 + 按秩合并](#路径压缩--按秩合并)
  - [带权并查集](#带权并查集)
  - [可撤销并查集](#可撤销并查集)
- [树状数组](#树状数组)
  - [单点改 + 区间查](#单点改--区间查)
  - [区间改 + 单点查（差分）](#区间改--单点查差分)
  - [区间改 + 区间查（双 BIT）](#区间改--区间查双-bit)
  - [二维树状数组（区间改区间查）](#二维树状数组区间改区间查)
  - [BIT 上二分](#bit-上二分)
- [ST 表](#st-表)
  - [一维 ST 表](#一维-st-表)
  - [二维 ST 表](#二维-st-表)
- [单调栈与单调队列](#单调栈与单调队列)
- [线段树](#线段树)
  - [综合 lazy 线段树（加 / 乘 / 赋值）](#综合-lazy-线段树加--乘--赋值)
  - [Segment Tree Beats（吉司机）](#segment-tree-beats吉司机)
  - [历史最值线段树](#历史最值线段树)
  - [动态开点线段树](#动态开点线段树)
  - [线段树合并 / 分裂](#线段树合并--分裂)
  - [可持久化线段树（主席树）](#可持久化线段树主席树)
  - [李超线段树](#李超线段树)
- [平衡树](#平衡树)
  - [FHQ Treap（无旋）](#fhq-treap无旋)
  - [Splay](#splay)
  - [pb\_ds Tree](#pb_ds-tree)
- [分块与莫队](#分块与莫队)
  - [分块模板](#分块模板)
  - [普通莫队](#普通莫队)
  - [带修莫队](#带修莫队)
  - [回滚莫队](#回滚莫队)
  - [树上莫队](#树上莫队)
- [高阶数据结构](#高阶数据结构)
  - [K-D Tree](#k-d-tree)
  - [笛卡尔树](#笛卡尔树)
  - [Link-Cut Tree](#link-cut-tree)
  - [线段树分治 + 可撤销并查集](#线段树分治--可撤销并查集)
- [复杂度速查](#复杂度速查)

---

## 并查集

### 路径压缩 + 按秩合并

```cpp
struct DSU {
    vector<int> fa, sz;
    int comp;
    DSU(int n) : fa(n), sz(n, 1), comp(n) { iota(fa.begin(), fa.end(), 0); }
    int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
    bool same(int x, int y) { return find(x) == find(y); }
    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (sz[x] < sz[y]) swap(x, y);
        fa[y] = x; sz[x] += sz[y]; --comp;
        return true;
    }
    int size(int x) { return sz[find(x)]; }
};
```

复杂度近似 $O(\alpha(n))$。

### 带权并查集

维护 $w[x]$ 表示 $x$ 到 `find(x)` 的相对权值（如距离、奇偶性、模 $k$ 余数）。

```cpp
struct WDSU {
    vector<int> fa; vector<long long> w; // w[x] = dist(x, fa[x])
    WDSU(int n) : fa(n), w(n, 0) { iota(fa.begin(), fa.end(), 0); }
    int find(int x) {
        if (fa[x] == x) return x;
        int r = find(fa[x]);
        w[x] += w[fa[x]];
        return fa[x] = r;
    }
    // 加约束 a - b = d (a 比 b 大 d)
    bool merge(int a, int b, long long d) {
        int ra = find(a), rb = find(b);
        if (ra == rb) return w[b] - w[a] == d; // 校验
        fa[ra] = rb; w[ra] = w[b] + d - w[a];
        return true;
    }
    long long diff(int a, int b) { // 仅在 same(a,b) 时有效
        find(a); find(b);
        return w[a] - w[b];
    }
};
```

经典题：食物链、奇偶性约束、二分图判定（边权 0/1，merge(u,v,1) 表示异色）。

### 可撤销并查集

只能用按秩合并（**禁止路径压缩**），栈记录每次合并修改，撤销时回退。线段树分治、动态图离线判二分图必备。

```cpp
struct RollbackDSU {
    vector<int> fa, sz;
    vector<tuple<int,int,int>> stk; // (y, fa[y], sz[x])
    int comp;
    RollbackDSU(int n) : fa(n), sz(n, 1), comp(n) { iota(fa.begin(), fa.end(), 0); }
    int find(int x) { while (fa[x] != x) x = fa[x]; return x; }
    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) { stk.emplace_back(-1, 0, 0); return false; }
        if (sz[x] < sz[y]) swap(x, y);
        stk.emplace_back(y, fa[y], sz[x]);
        fa[y] = x; sz[x] += sz[y]; --comp;
        return true;
    }
    void rollback() {
        auto [y, fy, sx] = stk.back(); stk.pop_back();
        if (y == -1) return;
        sz[fa[y]] = sx; fa[y] = fy; ++comp;
    }
};
```

复杂度 $O(\log n)$ 每次操作。

---

## 树状数组

### 单点改 + 区间查

```cpp
template<class T>
struct BIT {
    int n; vector<T> c;
    BIT(int n = 0) : n(n), c(n + 1, T{}) {}
    void add(int x, T v) { for (; x <= n; x += x & -x) c[x] += v; }
    T sum(int x) { T s{}; for (; x > 0; x -= x & -x) s += c[x]; return s; }
    T query(int l, int r) { return sum(r) - sum(l - 1); }
};
```

### 区间改 + 单点查（差分）

```cpp
// 维护差分数组 d[i] = a[i] - a[i-1]
// 区间 [l,r] 加 v：d[l] += v, d[r+1] -= v；单点查 a[x] = sum(d[1..x])
BIT<ll> d;
void update(int l, int r, ll v) { d.add(l, v); d.add(r + 1, -v); }
ll point(int x) { return d.sum(x); }
```

### 区间改 + 区间查（双 BIT）

$\sum_{i=1}^{p} a_i = (p+1)\sum d_i - \sum i \cdot d_i$。

```cpp
BIT<ll> b1, b2; // d[i], i*d[i]
void update(int l, int r, ll v) {
    b1.add(l, v); b1.add(r + 1, -v);
    b2.add(l, (ll)l * v); b2.add(r + 1, -(ll)(r + 1) * v);
}
ll prefix(int p) { return (ll)(p + 1) * b1.sum(p) - b2.sum(p); }
ll query(int l, int r) { return prefix(r) - prefix(l - 1); }
```

### 二维树状数组（区间改区间查）

四个 BIT 维护 $d, i\cdot d, j\cdot d, ij\cdot d$。

```cpp
struct BIT2D {
    int n, m;
    vector<vector<ll>> t[4];
    BIT2D(int n, int m) : n(n), m(m) {
        for (int k = 0; k < 4; k++) t[k].assign(n + 2, vector<ll>(m + 2, 0));
    }
    void add(int x, int y, ll v) {
        ll v1 = v, v2 = v * x, v3 = v * y, v4 = v * x * y;
        for (int i = x; i <= n; i += i & -i)
            for (int j = y; j <= m; j += j & -j)
                t[0][i][j] += v1, t[1][i][j] += v2, t[2][i][j] += v3, t[3][i][j] += v4;
    }
    void update(int x1, int y1, int x2, int y2, ll v) {
        add(x1, y1, v); add(x1, y2 + 1, -v); add(x2 + 1, y1, -v); add(x2 + 1, y2 + 1, v);
    }
    ll prefix(int x, int y) {
        ll s = 0;
        for (int i = x; i > 0; i -= i & -i)
            for (int j = y; j > 0; j -= j & -j)
                s += (ll)(x + 1) * (y + 1) * t[0][i][j]
                   - (y + 1) * t[1][i][j] - (x + 1) * t[2][i][j] + t[3][i][j];
        return s;
    }
    ll query(int x1, int y1, int x2, int y2) {
        return prefix(x2, y2) - prefix(x1 - 1, y2) - prefix(x2, y1 - 1) + prefix(x1 - 1, y1 - 1);
    }
};
```

### BIT 上二分

求最小 $p$ 满足 $\sum_{i\le p} a_i \ge k$（要求 $a_i \ge 0$）。$O(\log n)$。

```cpp
int kth(BIT<ll>& b, ll k) {
    int p = 0, lg = __lg(b.n);
    for (int i = lg; i >= 0; i--) {
        int nx = p + (1 << i);
        if (nx <= b.n && b.c[nx] < k) p = nx, k -= b.c[nx];
    }
    return p + 1;
}
```

---

## ST 表

### 一维 ST 表

支持可重复贡献的运算（min/max/gcd/按位与或）。

```cpp
template<class T, class Op = function<T(T,T)>>
struct SparseTable {
    int n; vector<vector<T>> st; Op op;
    SparseTable(const vector<T>& a, Op op) : n(a.size()), op(op) {
        int LOG = __lg(n) + 1;
        st.assign(LOG, vector<T>(n));
        st[0] = a;
        for (int k = 1; k < LOG; k++)
            for (int i = 0; i + (1 << k) <= n; i++)
                st[k][i] = op(st[k-1][i], st[k-1][i + (1 << (k-1))]);
    }
    T query(int l, int r) { // [l, r]
        int k = __lg(r - l + 1);
        return op(st[k][l], st[k][r - (1 << k) + 1]);
    }
};
// 用法: SparseTable<int> st(a, [](int x,int y){return min(x,y);});
```

预处理 $O(n\log n)$，查询 $O(1)$。

### 二维 ST 表

```cpp
struct ST2D {
    int n, m, LN, LM;
    vector<vector<vector<vector<int>>>> st;
    ST2D(const vector<vector<int>>& a) : n(a.size()), m(a[0].size()) {
        LN = __lg(n) + 1; LM = __lg(m) + 1;
        st.assign(LN, vector<vector<vector<int>>>(LM,
            vector<vector<int>>(n, vector<int>(m))));
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) st[0][0][i][j] = a[i][j];
        for (int kj = 1; kj < LM; kj++)
            for (int i = 0; i < n; i++)
                for (int j = 0; j + (1 << kj) <= m; j++)
                    st[0][kj][i][j] = max(st[0][kj-1][i][j], st[0][kj-1][i][j + (1 << (kj-1))]);
        for (int ki = 1; ki < LN; ki++)
            for (int kj = 0; kj < LM; kj++)
                for (int i = 0; i + (1 << ki) <= n; i++)
                    for (int j = 0; j < m; j++)
                        st[ki][kj][i][j] = max(st[ki-1][kj][i][j], st[ki-1][kj][i + (1 << (ki-1))][j]);
    }
    int query(int x1, int y1, int x2, int y2) {
        int ki = __lg(x2 - x1 + 1), kj = __lg(y2 - y1 + 1);
        return max({
            st[ki][kj][x1][y1],
            st[ki][kj][x2 - (1 << ki) + 1][y1],
            st[ki][kj][x1][y2 - (1 << kj) + 1],
            st[ki][kj][x2 - (1 << ki) + 1][y2 - (1 << kj) + 1]
        });
    }
};
```

预处理 $O(nm\log n\log m)$，查询 $O(1)$。

---

## 单调栈与单调队列

```cpp
// 单调栈：求每个位置左/右最近比它小/大的位置
vector<int> prev_less(const vector<int>& a) {
    int n = a.size(); vector<int> l(n, -1); stack<int> st;
    for (int i = 0; i < n; i++) {
        while (!st.empty() && a[st.top()] >= a[i]) st.pop();
        if (!st.empty()) l[i] = st.top();
        st.push(i);
    }
    return l;
}

// 单调队列：滑动窗口最小值，窗口大小 k
vector<int> sliding_min(const vector<int>& a, int k) {
    deque<int> dq; vector<int> res;
    for (int i = 0; i < (int)a.size(); i++) {
        while (!dq.empty() && a[dq.back()] >= a[i]) dq.pop_back();
        dq.push_back(i);
        if (dq.front() <= i - k) dq.pop_front();
        if (i >= k - 1) res.push_back(a[dq.front()]);
    }
    return res;
}
```

---

## 线段树

### 综合 lazy 线段树（加 / 乘 / 赋值）

支持 区间加、区间乘、区间赋值、区间求和。lazy 优先级：赋值 > 乘 > 加。

```cpp
const int MOD = 1'000'000'007;

struct SegTree {
    struct Node { ll sum, add, mul, set_; bool has_set; int len; };
    int n; vector<Node> t;
    SegTree(const vector<ll>& a) : n(a.size()), t(4 * n) { build(1, 0, n - 1, a); }

    void apply_set(int p, ll v) {
        t[p].sum = v % MOD * t[p].len % MOD;
        t[p].set_ = v; t[p].has_set = true;
        t[p].add = 0; t[p].mul = 1;
    }
    void apply_mul(int p, ll v) {
        t[p].sum = t[p].sum * v % MOD;
        t[p].mul = t[p].mul * v % MOD;
        t[p].add = t[p].add * v % MOD;
    }
    void apply_add(int p, ll v) {
        t[p].sum = (t[p].sum + v % MOD * t[p].len) % MOD;
        t[p].add = (t[p].add + v) % MOD;
    }
    void push(int p) {
        if (t[p].has_set) {
            apply_set(p << 1, t[p].set_);
            apply_set(p << 1 | 1, t[p].set_);
            t[p].has_set = false;
        }
        if (t[p].mul != 1) {
            apply_mul(p << 1, t[p].mul);
            apply_mul(p << 1 | 1, t[p].mul);
            t[p].mul = 1;
        }
        if (t[p].add) {
            apply_add(p << 1, t[p].add);
            apply_add(p << 1 | 1, t[p].add);
            t[p].add = 0;
        }
    }
    void pull(int p) { t[p].sum = (t[p << 1].sum + t[p << 1 | 1].sum) % MOD; }

    void build(int p, int l, int r, const vector<ll>& a) {
        t[p] = {0, 0, 1, 0, false, r - l + 1};
        if (l == r) { t[p].sum = a[l] % MOD; return; }
        int m = (l + r) >> 1;
        build(p << 1, l, m, a); build(p << 1 | 1, m + 1, r, a);
        pull(p);
    }

    void update(int p, int l, int r, int ql, int qr, int op, ll v) {
        if (ql <= l && r <= qr) {
            if (op == 0) apply_add(p, v);
            else if (op == 1) apply_mul(p, v);
            else apply_set(p, v);
            return;
        }
        push(p);
        int m = (l + r) >> 1;
        if (ql <= m) update(p << 1, l, m, ql, qr, op, v);
        if (qr > m) update(p << 1 | 1, m + 1, r, ql, qr, op, v);
        pull(p);
    }
    ll query(int p, int l, int r, int ql, int qr) {
        if (ql <= l && r <= qr) return t[p].sum;
        push(p);
        int m = (l + r) >> 1; ll s = 0;
        if (ql <= m) s = (s + query(p << 1, l, m, ql, qr)) % MOD;
        if (qr > m) s = (s + query(p << 1 | 1, m + 1, r, ql, qr)) % MOD;
        return s;
    }
};
```

### Segment Tree Beats（吉司机）

支持 **区间对 $v$ 取 min**、区间加、区间求和、区间最大值。
节点维护：最大值 `mx`、严格次大值 `se`、最大值出现次数 `cnt`、和。
区间取 min 时：若 $v \ge mx$ 跳过；若 $se < v < mx$ 直接更新；否则递归。

```cpp
struct SegBeats {
    struct Node { ll sum, mx, se, lz_add; int cnt; };
    int n; vector<Node> t;
    static const ll INF = (ll)4e18;

    SegBeats(const vector<ll>& a) : n(a.size()), t(4 * n) { build(1, 0, n - 1, a); }

    void pull(int p) {
        auto &L = t[p << 1], &R = t[p << 1 | 1];
        t[p].sum = L.sum + R.sum;
        if (L.mx == R.mx) {
            t[p].mx = L.mx; t[p].cnt = L.cnt + R.cnt;
            t[p].se = max(L.se, R.se);
        } else if (L.mx > R.mx) {
            t[p].mx = L.mx; t[p].cnt = L.cnt;
            t[p].se = max(L.se, R.mx);
        } else {
            t[p].mx = R.mx; t[p].cnt = R.cnt;
            t[p].se = max(L.mx, R.se);
        }
    }
    void apply_add(int p, int len, ll v) {
        t[p].sum += (ll)len * v; t[p].mx += v;
        if (t[p].se != -INF) t[p].se += v;
        t[p].lz_add += v;
    }
    void apply_chmin(int p, ll v) { // 仅当 v < mx 时调用
        if (v >= t[p].mx) return;
        t[p].sum -= (ll)(t[p].mx - v) * t[p].cnt;
        t[p].mx = v;
    }
    void push(int p, int l, int r) {
        int m = (l + r) >> 1;
        if (t[p].lz_add) {
            apply_add(p << 1, m - l + 1, t[p].lz_add);
            apply_add(p << 1 | 1, r - m, t[p].lz_add);
            t[p].lz_add = 0;
        }
        if (t[p].mx < t[p << 1].mx) apply_chmin(p << 1, t[p].mx);
        if (t[p].mx < t[p << 1 | 1].mx) apply_chmin(p << 1 | 1, t[p].mx);
    }
    void build(int p, int l, int r, const vector<ll>& a) {
        t[p] = {0, 0, -INF, 0, 1};
        if (l == r) { t[p].sum = t[p].mx = a[l]; return; }
        int m = (l + r) >> 1;
        build(p << 1, l, m, a); build(p << 1 | 1, m + 1, r, a);
        pull(p);
    }
    void chmin(int p, int l, int r, int ql, int qr, ll v) {
        if (qr < l || r < ql || v >= t[p].mx) return;
        if (ql <= l && r <= qr && v > t[p].se) { apply_chmin(p, v); return; }
        push(p, l, r);
        int m = (l + r) >> 1;
        chmin(p << 1, l, m, ql, qr, v); chmin(p << 1 | 1, m + 1, r, ql, qr, v);
        pull(p);
    }
    void add(int p, int l, int r, int ql, int qr, ll v) {
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr) { apply_add(p, r - l + 1, v); return; }
        push(p, l, r);
        int m = (l + r) >> 1;
        add(p << 1, l, m, ql, qr, v); add(p << 1 | 1, m + 1, r, ql, qr, v);
        pull(p);
    }
    ll qsum(int p, int l, int r, int ql, int qr) {
        if (qr < l || r < ql) return 0;
        if (ql <= l && r <= qr) return t[p].sum;
        push(p, l, r);
        int m = (l + r) >> 1;
        return qsum(p << 1, l, m, ql, qr) + qsum(p << 1 | 1, m + 1, r, ql, qr);
    }
};
```

均摊复杂度 $O(n \log^2 n)$。同时支持 chmax 时还需对称维护最小值与严格次小值。

### 历史最值线段树

维护历史最大值 $mxh$ 与历史最大加法标记 $lz_h$。每次 push 时：$mxh_{son} = \max(mxh_{son}, mx_{son} + lz_h)$，再合并 lazy。

```cpp
struct HistSeg {
    struct Node { ll mx, mxh, lz, lzh; }; // lz 当前加, lzh 历史加最大值
    int n; vector<Node> t;

    void apply(int p, ll a, ll ah) {
        t[p].mxh = max(t[p].mxh, t[p].mx + ah);
        t[p].lzh = max(t[p].lzh, t[p].lz + ah);
        t[p].mx += a; t[p].lz += a;
    }
    void push(int p) {
        apply(p << 1, t[p].lz, t[p].lzh);
        apply(p << 1 | 1, t[p].lz, t[p].lzh);
        t[p].lz = t[p].lzh = 0;
    }
    void pull(int p) {
        t[p].mx = max(t[p << 1].mx, t[p << 1 | 1].mx);
        t[p].mxh = max(t[p << 1].mxh, t[p << 1 | 1].mxh);
    }
    // build / update / query 同标准线段树, lazy 用 apply(p, v, max(0, v))
};
```

### 动态开点线段树

值域大或区间稀疏时使用，节点按需分配。

```cpp
struct DynSeg {
    struct Node { int ls, rs; ll sum, lz; };
    vector<Node> t; int rt = 0;
    DynSeg() { t.push_back({0, 0, 0, 0}); } // 0 号为 null
    int newNode() { t.push_back({0, 0, 0, 0}); return t.size() - 1; }
    void apply(int p, ll len, ll v) { t[p].sum += len * v; t[p].lz += v; }
    void push(int p, ll l, ll r) {
        if (!t[p].lz) return;
        ll m = (l + r) >> 1;
        if (!t[p].ls) t[p].ls = newNode();
        if (!t[p].rs) t[p].rs = newNode();
        apply(t[p].ls, m - l + 1, t[p].lz);
        apply(t[p].rs, r - m, t[p].lz);
        t[p].lz = 0;
    }
    void update(int& p, ll l, ll r, ll ql, ll qr, ll v) {
        if (!p) p = newNode();
        if (ql <= l && r <= qr) { apply(p, r - l + 1, v); return; }
        push(p, l, r);
        ll m = (l + r) >> 1;
        if (ql <= m) update(t[p].ls, l, m, ql, qr, v);
        if (qr > m) update(t[p].rs, m + 1, r, ql, qr, v);
        t[p].sum = t[t[p].ls].sum + t[t[p].rs].sum;
    }
    ll query(int p, ll l, ll r, ll ql, ll qr) {
        if (!p || qr < l || r < ql) return 0;
        if (ql <= l && r <= qr) return t[p].sum;
        push(p, l, r);
        ll m = (l + r) >> 1;
        return query(t[p].ls, l, m, ql, qr) + query(t[p].rs, m + 1, r, ql, qr);
    }
};
```

### 线段树合并 / 分裂

权值线段树合并：把两棵权值线段树合成一棵。常用于树上启发式合并（如「Lomsat gelral」）。

```cpp
// 合并：返回合并后的根
int merge(int x, int y, int l, int r) {
    if (!x || !y) return x | y;
    if (l == r) { t[x].sum += t[y].sum; return x; }
    int m = (l + r) >> 1;
    t[x].ls = merge(t[x].ls, t[y].ls, l, m);
    t[x].rs = merge(t[x].rs, t[y].rs, m + 1, r);
    t[x].sum = t[t[x].ls].sum + t[t[x].rs].sum;
    return x;
}

// 分裂：把 x 中 [k+1, r] 的部分分给 y
void split(int x, int& y, ll k) { // k = 左侧应保留的元素数
    if (!x) return;
    y = newNode();
    ll s = t[t[x].ls].sum;
    if (k > s) split(t[x].rs, t[y].rs, k - s);
    else swap(t[x].rs, t[y].rs);
    if (k < s) split(t[x].ls, t[y].ls, k);
    t[y].sum = t[x].sum - k;
    t[x].sum = k;
}
```

合并复杂度 $O(\text{合并节点数})$，整体 $O((n+m)\log n)$。

### 可持久化线段树（主席树）

按时间维护多版本权值线段树，每次修改新建 $O(\log n)$ 节点。区间第 K 小：在 $\text{root}_r - \text{root}_{l-1}$ 上二分。

```cpp
struct PersistSeg {
    struct Node { int ls, rs, cnt; };
    vector<Node> t; vector<int> root;
    int n;
    PersistSeg(int n) : n(n) { t.push_back({0, 0, 0}); root.push_back(0); }

    int update(int pre, int l, int r, int x) {
        t.push_back(t[pre]);
        int p = t.size() - 1;
        ++t[p].cnt;
        if (l == r) return p;
        int m = (l + r) >> 1;
        if (x <= m) t[p].ls = update(t[pre].ls, l, m, x);
        else t[p].rs = update(t[pre].rs, m + 1, r, x);
        return p;
    }
    void insert(int x) { root.push_back(update(root.back(), 1, n, x)); }
    // 查询 [l, r]（版本号区间）的第 k 小
    int kth(int u, int v, int l, int r, int k) {
        if (l == r) return l;
        int m = (l + r) >> 1;
        int lc = t[t[v].ls].cnt - t[t[u].ls].cnt;
        if (k <= lc) return kth(t[u].ls, t[v].ls, l, m, k);
        return kth(t[u].rs, t[v].rs, m + 1, r, k - lc);
    }
    int kth(int l, int r, int k) { return kth(root[l - 1], root[r], 1, n, k); }
};
```

值需先离散化，$n$ 为离散化后大小。

### 李超线段树

支持插入 直线 / 线段，单点查 $\max$（或 $\min$）。每个区间存当前**最优**直线，递归下放被覆盖的次优直线。

```cpp
struct Line {
    double k, b;
    double eval(int x) const { return k * x + b; }
};

struct LiChao {
    int n; vector<int> id; vector<Line> seg;
    LiChao(int n) : n(n), id(4 * n, 0) { seg.push_back({0, -1e18}); } // 0 号为 -inf 直线

    void modify(int p, int l, int r, int ql, int qr, int u) {
        if (ql > r || qr < l) return;
        int m = (l + r) >> 1;
        if (ql <= l && r <= qr) {
            if (seg[u].eval(m) > seg[id[p]].eval(m)) swap(u, id[p]);
            if (l == r) return;
            if (seg[u].eval(l) > seg[id[p]].eval(l)) modify(p << 1, l, m, ql, qr, u);
            else if (seg[u].eval(r) > seg[id[p]].eval(r)) modify(p << 1 | 1, m + 1, r, ql, qr, u);
            return;
        }
        modify(p << 1, l, m, ql, qr, u);
        modify(p << 1 | 1, m + 1, r, ql, qr, u);
    }
    void add_line(double k, double b) {
        seg.push_back({k, b});
        modify(1, 1, n, 1, n, seg.size() - 1);
    }
    void add_segment(double k, double b, int l, int r) {
        seg.push_back({k, b});
        modify(1, 1, n, l, r, seg.size() - 1);
    }
    double query(int p, int l, int r, int x) {
        double res = seg[id[p]].eval(x);
        if (l == r) return res;
        int m = (l + r) >> 1;
        if (x <= m) res = max(res, query(p << 1, l, m, x));
        else res = max(res, query(p << 1 | 1, m + 1, r, x));
        return res;
    }
    double query(int x) { return query(1, 1, n, x); }
};
```

直线插入 $O(\log n)$，线段插入 $O(\log^2 n)$，查询 $O(\log n)$。整数版本注意比较时开 `long long` 防溢出。

---

## 平衡树

### FHQ Treap（无旋）

通用平衡树，按值分裂或按大小分裂；支持区间翻转、合并、可持久化。

```cpp
struct FHQ {
    struct Node {
        int ls, rs, val, pri, sz;
        int rev; // 区间翻转
    };
    vector<Node> t; mt19937 rng{(unsigned)chrono::steady_clock::now().time_since_epoch().count()};
    int root = 0;
    FHQ() { t.push_back({0, 0, 0, 0, 0, 0}); } // 哨兵

    int newNode(int v) { t.push_back({0, 0, v, (int)rng(), 1, 0}); return t.size() - 1; }
    void pull(int p) { t[p].sz = t[t[p].ls].sz + t[t[p].rs].sz + 1; }
    void apply_rev(int p) { swap(t[p].ls, t[p].rs); t[p].rev ^= 1; }
    void push(int p) {
        if (t[p].rev) {
            if (t[p].ls) apply_rev(t[p].ls);
            if (t[p].rs) apply_rev(t[p].rs);
            t[p].rev = 0;
        }
    }
    // 按值分裂：≤ v 入 x，> v 入 y
    void splitVal(int p, int v, int& x, int& y) {
        if (!p) { x = y = 0; return; }
        push(p);
        if (t[p].val <= v) { x = p; splitVal(t[p].rs, v, t[p].rs, y); }
        else { y = p; splitVal(t[p].ls, v, x, t[p].ls); }
        pull(p);
    }
    // 按大小分裂：前 k 个入 x
    void splitSz(int p, int k, int& x, int& y) {
        if (!p) { x = y = 0; return; }
        push(p);
        if (t[t[p].ls].sz + 1 <= k) { x = p; splitSz(t[p].rs, k - t[t[p].ls].sz - 1, t[p].rs, y); }
        else { y = p; splitSz(t[p].ls, k, x, t[p].ls); }
        pull(p);
    }
    int merge(int a, int b) {
        if (!a || !b) return a | b;
        if (t[a].pri < t[b].pri) { push(a); t[a].rs = merge(t[a].rs, b); pull(a); return a; }
        else { push(b); t[b].ls = merge(a, t[b].ls); pull(b); return b; }
    }
    void insert(int v) {
        int x, y; splitVal(root, v - 1, x, y);
        root = merge(merge(x, newNode(v)), y);
    }
    void erase(int v) {
        int x, y, z; splitVal(root, v - 1, x, y); splitVal(y, v, y, z);
        y = merge(t[y].ls, t[y].rs);
        root = merge(merge(x, y), z);
    }
    int rank(int v) { // 比 v 小的元素数 + 1
        int x, y; splitVal(root, v - 1, x, y);
        int r = t[x].sz + 1;
        root = merge(x, y); return r;
    }
    int kth(int p, int k) {
        push(p);
        if (t[t[p].ls].sz + 1 == k) return t[p].val;
        if (t[t[p].ls].sz >= k) return kth(t[p].ls, k);
        return kth(t[p].rs, k - t[t[p].ls].sz - 1);
    }
    // 区间 [l, r] 翻转（用作序列平衡树，按下标 splitSz）
    void reverse(int l, int r) {
        int a, b, c; splitSz(root, l - 1, a, b); splitSz(b, r - l + 1, b, c);
        apply_rev(b);
        root = merge(merge(a, b), c);
    }
};
```

可持久化：`splitVal` / `merge` 中每个被修改的节点先 `t.push_back(t[p])` 复制，记录新版本根。

### Splay

经典区间平衡树，常数大但实现紧凑；序列翻转、LCT 实现核心。

```cpp
struct Splay {
    struct Node { int fa, ch[2], val, sz, cnt, rev; };
    vector<Node> t; int root = 0;
    Splay() { t.push_back({0, {0, 0}, 0, 0, 0, 0}); }
    int newNode(int v, int fa) { t.push_back({fa, {0, 0}, v, 1, 1, 0}); return t.size() - 1; }
    void pull(int p) { t[p].sz = t[t[p].ch[0]].sz + t[t[p].ch[1]].sz + t[p].cnt; }
    void apply_rev(int p) { swap(t[p].ch[0], t[p].ch[1]); t[p].rev ^= 1; }
    void push(int p) {
        if (t[p].rev) {
            if (t[p].ch[0]) apply_rev(t[p].ch[0]);
            if (t[p].ch[1]) apply_rev(t[p].ch[1]);
            t[p].rev = 0;
        }
    }
    bool dir(int x) { return t[t[x].fa].ch[1] == x; }
    void rotate(int x) {
        int y = t[x].fa, z = t[y].fa, d = dir(x), w = t[x].ch[d ^ 1];
        if (z) t[z].ch[dir(y)] = x;
        t[x].ch[d ^ 1] = y; t[y].ch[d] = w;
        if (w) t[w].fa = y;
        t[y].fa = x; t[x].fa = z;
        pull(y); pull(x);
    }
    void splay(int x, int goal = 0) {
        while (t[x].fa != goal) {
            int y = t[x].fa, z = t[y].fa;
            if (z != goal) rotate(dir(x) == dir(y) ? y : x);
            rotate(x);
        }
        if (!goal) root = x;
    }
    // 其余 insert/erase/kth/find 略
};
```

### pb\_ds Tree

GNU 扩展，提供 `order_of_key` / `find_by_order`，比赛允许时省事。

```cpp
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag,
                         tree_order_statistics_node_update>;

ordered_set<int> s;
s.insert(5);
auto it = s.find_by_order(0); // 第 1 小
int rk = s.order_of_key(5);   // 严格小于 5 的个数
```

---

## 分块与莫队

### 分块模板

块大小 $B = \sqrt n$，每块维护可加 lazy 与块内元素副本（区间排序时用）。

```cpp
struct Block {
    int n, B; vector<ll> a, tag;
    vector<vector<ll>> blk; // 每块排序副本
    int bid(int i) { return i / B; }
    Block(vector<ll> a) : a(a) {
        n = a.size(); B = max(1, (int)sqrt((double)n));
        int cnt = (n + B - 1) / B;
        tag.assign(cnt, 0); blk.assign(cnt, {});
        for (int i = 0; i < n; i++) blk[bid(i)].push_back(a[i]);
        for (auto& v : blk) sort(v.begin(), v.end());
    }
    void rebuild(int b) {
        blk[b].clear();
        int l = b * B, r = min(n, l + B);
        for (int i = l; i < r; i++) blk[b].push_back(a[i]);
        sort(blk[b].begin(), blk[b].end());
    }
    void update(int l, int r, ll v) {
        int bl = bid(l), br = bid(r);
        if (bl == br) {
            for (int i = l; i <= r; i++) a[i] += v;
            rebuild(bl); return;
        }
        for (int i = l; i < (bl + 1) * B; i++) a[i] += v;
        rebuild(bl);
        for (int b = bl + 1; b < br; b++) tag[b] += v;
        for (int i = br * B; i <= r; i++) a[i] += v;
        rebuild(br);
    }
    // 查询「区间内 < x 的个数」示例
    int countLess(int l, int r, ll x) {
        int bl = bid(l), br = bid(r), res = 0;
        if (bl == br) {
            for (int i = l; i <= r; i++) if (a[i] + tag[bl] < x) ++res;
            return res;
        }
        for (int i = l; i < (bl + 1) * B; i++) if (a[i] + tag[bl] < x) ++res;
        for (int b = bl + 1; b < br; b++)
            res += lower_bound(blk[b].begin(), blk[b].end(), x - tag[b]) - blk[b].begin();
        for (int i = br * B; i <= r; i++) if (a[i] + tag[br] < x) ++res;
        return res;
    }
};
```

### 普通莫队

排序后双指针离线处理区间询问。块大小取 $\sqrt n$，奇偶分块（奇数块 $r$ 升序 / 偶数块 $r$ 降序）减少常数。

```cpp
struct Query { int l, r, id; };
int B;
bool cmp(const Query& a, const Query& b) {
    if (a.l / B != b.l / B) return a.l < b.l;
    return (a.l / B) & 1 ? a.r < b.r : a.r > b.r;
}

void mo(vector<int>& a, vector<Query>& qs, vector<long long>& ans) {
    int n = a.size(); B = max(1, (int)sqrt((double)n));
    sort(qs.begin(), qs.end(), cmp);
    int l = 1, r = 0; long long cur = 0;
    auto add = [&](int x) { /* ... 更新 cur */ };
    auto del = [&](int x) { /* ... 更新 cur */ };
    for (auto& q : qs) {
        while (r < q.r) add(a[++r]);
        while (l > q.l) add(a[--l]);
        while (r > q.r) del(a[r--]);
        while (l < q.l) del(a[l++]);
        ans[q.id] = cur;
    }
}
```

复杂度 $O((n + q)\sqrt n)$。

### 带修莫队

加入「时间」维：每次修改是一次时间推进。块大小取 $n^{2/3}$，三关键字排序。

```cpp
struct Q { int l, r, t, id; };
struct U { int p, v, old; };
int B;
bool cmp(const Q& a, const Q& b) {
    if (a.l / B != b.l / B) return a.l / B < b.l / B;
    if (a.r / B != b.r / B) return a.r / B < b.r / B;
    return a.t < b.t;
}
// 维护当前时间 t，apply(t)/revoke(t) 应用/撤销修改
```

复杂度 $O(n^{5/3})$。

### 回滚莫队

适用于「只能加不能删」型问题（如区间最大值、区间众数）。同块的询问左指针在块右边界处归位重做，右指针单调右移。

```cpp
// 排序：左端点同块时按右端点升序
sort(qs.begin(), qs.end(), [&](const Query& a, const Query& b) {
    if (a.l / B != b.l / B) return a.l < b.l;
    return a.r < b.r;
});
// 同块内：l 指针每次回到块右边界 + 1，从该位置重新延伸
// 块内单独询问 (l, r 同块)：暴力扫描
```

复杂度 $O((n+q)\sqrt n)$。

### 树上莫队

将树拍成欧拉序（每个点出现两次：进入 in[u] 与离开 out[u]）。
区间 $[in[u], in[v]]$（$u$ 是 $v$ 祖先时）或 $[out[u], in[v]]$ + LCA 单独处理（一般情况）。

```cpp
// 入栈 in[u] 时若已访问则删除，未访问则添加；出栈 out[u] 同理 ⇒ 出现奇数次的点正是路径上的点。
auto toggle = [&](int u) { if (vis[u]) del(a[u]); else add(a[u]); vis[u] ^= 1; };
// 处理 LCA 不在 [in[u], in[v]] 中的情形：toggle(lca) 前后各一次
```

---

## 高阶数据结构

### K-D Tree

二维点集近邻查询、矩形询问、带替换/删除（懒重构）。每隔 $\sqrt{\text{节点数}}$ 重构保持平衡。

```cpp
struct KD {
    struct Node { int ls, rs; int p[2], mn[2], mx[2]; ll v, sum; bool del; };
    vector<Node> t; int rt = 0;
    int cmp_dim;

    int build(vector<int>& idx, int l, int r, int d) {
        if (l > r) return 0;
        int m = (l + r) >> 1;
        cmp_dim = d;
        nth_element(idx.begin() + l, idx.begin() + m, idx.begin() + r + 1,
            [&](int a, int b){ return t[a].p[cmp_dim] < t[b].p[cmp_dim]; });
        int p = idx[m];
        t[p].ls = build(idx, l, m - 1, d ^ 1);
        t[p].rs = build(idx, m + 1, r, d ^ 1);
        for (int k = 0; k < 2; k++) t[p].mn[k] = t[p].mx[k] = t[p].p[k];
        t[p].sum = t[p].v;
        for (int c : {t[p].ls, t[p].rs}) if (c) {
            for (int k = 0; k < 2; k++) {
                t[p].mn[k] = min(t[p].mn[k], t[c].mn[k]);
                t[p].mx[k] = max(t[p].mx[k], t[c].mx[k]);
            }
            t[p].sum += t[c].sum;
        }
        return p;
    }

    // 矩形 [x1,x2] x [y1,y2] 内权值和
    ll query(int p, int x1, int x2, int y1, int y2) {
        if (!p || t[p].mx[0] < x1 || t[p].mn[0] > x2
              || t[p].mx[1] < y1 || t[p].mn[1] > y2) return 0;
        if (x1 <= t[p].mn[0] && t[p].mx[0] <= x2
         && y1 <= t[p].mn[1] && t[p].mx[1] <= y2) return t[p].sum;
        ll s = 0;
        if (!t[p].del && x1 <= t[p].p[0] && t[p].p[0] <= x2
                     && y1 <= t[p].p[1] && t[p].p[1] <= y2) s += t[p].v;
        s += query(t[p].ls, x1, x2, y1, y2);
        s += query(t[p].rs, x1, x2, y1, y2);
        return s;
    }
};
```

矩形询问 $O(\sqrt n)$，最近点 $O(\log n)$ 期望。

### 笛卡尔树

按下标作 BST、按权值作堆。$O(n)$ 单调栈构造，常用于 RMQ ↔ LCA 转化、最大子矩形。

```cpp
// 小根笛卡尔树：a[i] 为权值
vector<int> ls(n, -1), rs(n, -1); stack<int> stk;
int root = 0;
for (int i = 0; i < n; i++) {
    int last = -1;
    while (!stk.empty() && a[stk.top()] > a[i]) { last = stk.top(); stk.pop(); }
    if (!stk.empty()) rs[stk.top()] = i;
    else root = i;
    ls[i] = last;
    stk.push(i);
}
```

### Link-Cut Tree

动态树问题（树边插删 + 路径查询）。底层是若干棵 Splay，按 preferred path 维护 prefer-child 关系。

```cpp
struct LCT {
    struct Node { int fa, ch[2], rev; ll val, sum; };
    vector<Node> t;
    LCT(int n) : t(n + 1) {}

    bool isRoot(int x) { return t[t[x].fa].ch[0] != x && t[t[x].fa].ch[1] != x; }
    bool dir(int x) { return t[t[x].fa].ch[1] == x; }
    void pull(int x) { t[x].sum = t[t[x].ch[0]].sum ^ t[x].val ^ t[t[x].ch[1]].sum; }
    void apply_rev(int x) { swap(t[x].ch[0], t[x].ch[1]); t[x].rev ^= 1; }
    void push(int x) {
        if (t[x].rev) {
            if (t[x].ch[0]) apply_rev(t[x].ch[0]);
            if (t[x].ch[1]) apply_rev(t[x].ch[1]);
            t[x].rev = 0;
        }
    }
    void pushPath(int x) {
        if (!isRoot(x)) pushPath(t[x].fa);
        push(x);
    }
    void rotate(int x) {
        int y = t[x].fa, z = t[y].fa, d = dir(x), w = t[x].ch[d ^ 1];
        if (!isRoot(y)) t[z].ch[dir(y)] = x;
        t[x].ch[d ^ 1] = y; t[y].ch[d] = w;
        if (w) t[w].fa = y;
        t[y].fa = x; t[x].fa = z;
        pull(y); pull(x);
    }
    void splay(int x) {
        pushPath(x);
        while (!isRoot(x)) {
            int y = t[x].fa;
            if (!isRoot(y)) rotate(dir(x) == dir(y) ? y : x);
            rotate(x);
        }
    }
    void access(int x) {
        for (int y = 0; x; y = x, x = t[x].fa) {
            splay(x); t[x].ch[1] = y; pull(x);
        }
    }
    void makeRoot(int x) { access(x); splay(x); apply_rev(x); }
    int findRoot(int x) {
        access(x); splay(x);
        push(x);
        while (t[x].ch[0]) {
            x = t[x].ch[0];
            push(x);
        }

        splay(x); return x;
    }

    bool connected(int x, int y) {
        if (x == y) return true;
        return findRoot(x) == findRoot(y);
    }

    void link(int x, int y) {
        makeRoot(x);
        if (findRoot(y) == x) return; // 已连通
        t[x].fa = y;
    }
    void cut(int x, int y) {
        makeRoot(x); access(y); splay(y);
        if (t[y].ch[0] != x || t[x].ch[1]) return; // 不直接相邻
        t[y].ch[0] = t[x].fa = 0; pull(y);
    }
    ll pathSum(int x, int y) { makeRoot(x); access(y); splay(y); return t[y].sum; }
    void modify(int x, ll v) { access(x); splay(x); t[x].val = v; pull(x); }
};
```

```
// 有根树的写法
bool link(int x, int y) {
    if (findRoot(x) != x) return false;
    if (findRoot(y) == x) return false;
    t[x].fa = y;
    return true;
}
void cut(int x) {
    access(x);
    splay(x);

    int left = t[x].ch[0];
    t[x].ch[0] = 0;
    t[left].fa = 0;
    pull(x);
}
```

```
//有根树的完整模板
struct LCT {
    struct Node { int fa, ch[2]; ll val, sum; };
    vector<Node> t;

    LCT(int n) : t(n + 1) {}

    bool isRoot(int x) {
        int f = t[x].fa;
        return t[f].ch[0] != x && t[f].ch[1] != x;
    }

    bool dir(int x) { return t[t[x].fa].ch[1] == x; }

    void pull(int x) {
        t[x].sum = t[t[x].ch[0]].sum ^ t[x].val ^ t[t[x].ch[1]].sum;
    }

    void rotate(int x) {
        int y = t[x].fa, z = t[y].fa, d = dir(x), w = t[x].ch[d ^ 1];
        if (!isRoot(y)) t[z].ch[dir(y)] = x;
        t[x].ch[d ^ 1] = y; t[y].ch[d] = w;
        if (w) t[w].fa = y;
        t[y].fa = x; t[x].fa = z;
        pull(y); pull(x);
    }

    void splay(int x) {
        while (!isRoot(x)) {
            int y = t[x].fa;
            if (!isRoot(y)) rotate(dir(x) == dir(y) ? y : x);
            rotate(x);
        }
    }

    int access(int x) {
        int y = 0;
        for (; x; x = t[y = x].fa) {
            splay(x);
            t[x].ch[1] = y;
            pull(x);
        }
        return y;
    }

    int findRoot(int x) {
        access(x); splay(x);
        while (t[x].ch[0]) x = t[x].ch[0];
        splay(x);
        return x;
    }

    bool connected(int x, int y) {
        return x == y || findRoot(x) == findRoot(y);
    }

    int parent(int x) {
        access(x); splay(x);
        if (!t[x].ch[0]) return 0;
        int y = t[x].ch[0];
        while (t[y].ch[1]) y = t[y].ch[1];
        splay(y);
        return y;
    }

    bool link(int x, int y) {
        if (findRoot(x) != x || findRoot(y) == x) return false;
        t[x].fa = y;
        return true;
    }

    bool cut(int x) {
        access(x); splay(x);
        if (!t[x].ch[0]) return false;
        int y = t[x].ch[0];
        t[x].ch[0] = 0;
        t[y].fa = 0;
        pull(x);
        return true;
    }

    bool cut(int x, int y) {
        if (parent(x) != y) return false;
        return cut(x);
    }

    int lca(int x, int y) {
        if (!connected(x, y)) return 0;
        access(x);
        return access(y);
    }

    ll rootXor(int x) {
        access(x); splay(x);
        return t[x].sum;
    }

    ll pathXor(int x, int y) {
        int z = lca(x, y);
        if (!z) return 0;
        return rootXor(x) ^ rootXor(y) ^ t[z].val;
    }

    void modify(int x, ll v) {
        access(x); splay(x);
        t[x].val = v;
        pull(x);
    }

    void setVal(int x, ll v) {
        t[x].val = t[x].sum = v;
    }

    ll rootSum(int x) {
        access(x); splay(x);
        return t[x].sum;
    }
    
    ll pathSum(int x, int y) {
        int z = lca(x, y);
        if (!z) return 0;
        return rootSum(x) + rootSum(y) - 2 * rootSum(z) + t[z].val;
    }
};
```

均摊 $O(\log n)$ 每次操作。

### 线段树分治 + 可撤销并查集

适用于「边在某段时间存在」的动态连通性 / 二分图判定 / MST 类问题：把每条边的存在区间挂到时间线段树对应节点上，DFS 该树时 merge 进入、回溯时 rollback。

```cpp
struct SegDiv {
    int T; vector<vector<pair<int,int>>> seg; // 每个时间节点上的边集
    SegDiv(int T) : T(T), seg(4 * T) {}
    void add(int p, int l, int r, int ql, int qr, pair<int,int> e) {
        if (qr < l || r < ql) return;
        if (ql <= l && r <= qr) { seg[p].push_back(e); return; }
        int m = (l + r) >> 1;
        add(p << 1, l, m, ql, qr, e); add(p << 1 | 1, m + 1, r, ql, qr, e);
    }
    void add(int l, int r, pair<int,int> e) { add(1, 1, T, l, r, e); }

    RollbackDSU dsu{0}; // 外部初始化
    void solve(int p, int l, int r) {
        int snap = dsu.stk.size();
        for (auto [u, v] : seg[p]) dsu.merge(u, v);
        if (l == r) {
            // 在时刻 l 回答询问，例如输出 dsu.comp
        } else {
            int m = (l + r) >> 1;
            solve(p << 1, l, m); solve(p << 1 | 1, m + 1, r);
        }
        while ((int)dsu.stk.size() > snap) dsu.rollback();
    }
};
```

总复杂度 $O(m \log T \log n)$，$T$ 为时间长度。

---

## 复杂度速查

| 数据结构 | 单次操作 | 备注 |
|---|---|---|
| 并查集（路径压缩+按秩） | $O(\alpha(n))$ | 不可撤销 |
| 可撤销并查集 | $O(\log n)$ | 仅按秩合并 |
| BIT | $O(\log n)$ | 区改区查需双 BIT |
| ST 表 | 预处理 $O(n\log n)$，查 $O(1)$ | 仅幂等运算 |
| 线段树 | $O(\log n)$ | 多 lazy 注意优先级 |
| Segment Tree Beats | 均摊 $O(\log^2 n)$ | chmin/chmax + 求和 |
| 动态开点 | $O(\log V)$ 空间随操作数 | 值域大时使用 |
| 主席树 | $O(\log n)$ 空间 $O(n\log n)$ | 区间第 K 小 |
| 李超树 | 直线 $O(\log n)$ / 线段 $O(\log^2 n)$ | DP 斜率优化 |
| FHQ Treap | 期望 $O(\log n)$ | 易写，可持久化 |
| Splay | 均摊 $O(\log n)$ | 区间翻转 / LCT 内核 |
| 普通莫队 | $O((n+q)\sqrt n)$ | 离线 |
| 带修莫队 | $O(n^{5/3})$ | 块 $n^{2/3}$ |
| 回滚莫队 | $O((n+q)\sqrt n)$ | 不可删元素 |
| K-D Tree | 矩形 $O(\sqrt n)$ | 周期重构防偏 |
| LCT | 均摊 $O(\log n)$ | 动态树 |
| 线段树分治 | $O(m\log T\log n)$ | 离线动态连通 |
