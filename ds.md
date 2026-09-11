# 数据结构

## 目录
- [01 trie](#01-trie)
- [Segment Tree Beats（吉司机）](#segment-tree-beats吉司机)
- [李超线段树](#李超线段树)
- [Splay](#splay)
- [普通莫队](#普通莫队)
- [带修莫队](#带修莫队)
- [回滚莫队](#回滚莫队)
- [树上莫队](#树上莫队)
- [K-D Tree](#k-d-tree)
- [笛卡尔树](#笛卡尔树)
- [Link-Cut Tree](#link-cut-tree)

---

## 01 trie

```
struct Trie {
    static const int B = 60;
    int tot = 1;
    vector<array<int, 2>> ch;
    vector<int> cnt;

    Trie(int n) : ch(n * (B + 2) + 5), cnt(n * (B + 2) + 5) {}

    void insert(long long x) {
        int p = 1;
        cnt[p]++;
        for(int i = B; i >= 0; i--) {
            int c = (x >> i) & 1;
            if(!ch[p][c]) ch[p][c] = ++tot;
            p = ch[p][c];
            cnt[p]++;
        }
    }

    void erase(long long x) {
        int p = 1;
        cnt[p]--;
        for(int i = B; i >= 0; i--) {
            int c = (x >> i) & 1;
            p = ch[p][c];
            cnt[p]--;
        }
    }

    long long queryMax(long long x) {
        int p = 1;
        long long ans = 0;
        for(int i = B; i >= 0; i--) {
            int c = (x >> i) & 1;
            if(ch[p][c ^ 1] && cnt[ch[p][c ^ 1]]) {
                ans |= 1LL << i;
                p = ch[p][c ^ 1];
            } else {
                p = ch[p][c];
            }
        }
        return ans;
    }

    long long queryMin(long long x) {
        int p = 1;
        long long ans = 0;
        for(int i = B; i >= 0; i--) {
            int c = (x >> i) & 1;
            if(ch[p][c] && cnt[ch[p][c]]) {
                p = ch[p][c];
            } else {
                ans |= 1LL << i;
                p = ch[p][c ^ 1];
            }
        }
        return ans;
    }
};
```

---

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

---

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
//无根树模板
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
// 有根树的写法(在无根树基础上改)
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
