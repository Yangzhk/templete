
## 目录

- [目录](#目录)
  - [模运算重载](#模运算重载)
  - [模运算](#模运算)
  - [矩阵运算类](#矩阵运算类)
  - [普通高斯消元](#普通高斯消元)
  - [行列式](#行列式)
  - [莫队](#莫队)
  - [xor 线性基](#xor-线性基)
  - [异或高斯消元](#异或高斯消元)
  - [最大流](#最大流)
  - [zkw 费用流](#zkw-费用流)
  - [dij 费用流](#dij-费用流)
  - [上下界网络流](#上下界网络流)
  - [计算几何 (dls version)](#计算几何-dls-version)
  - [三维点类](#三维点类)
  - [集合操作（位运算）](#集合操作位运算)
  - [三元环计数](#三元环计数)
  - [四元环计数](#四元环计数)
  - [弱周期引理](#弱周期引理)
  - [fail 链路径压缩](#fail-链路径压缩)
  - [ntt](#ntt)
---

### 模运算重载
```
struct Z {
    ll x;
    Z(ll v = 0) : x((v % mod + mod) % mod) {}

    Z operator + (Z b) const { return x + b.x; }
    Z operator - (Z b) const { return x - b.x; }
    Z operator * (Z b) const { return x * b.x; }

    Z qpow(ll b) const {
        Z a = *this, res = 1;
        for(; b; b >>= 1, a = a * a)
            if(b & 1) res = res * a;
        return res;
    }

    Z operator / (Z b) const {
        return *this * b.qpow(mod - 2);
    }
};
```

### 模运算

```
int fac[N], inv[N];

int qpow(int a, int v = mod - 2){
    int res = 1;
    for(; v; a = a * a % mod, v >>= 1){
        if(v & 1) res = res * a % mod;
    }
    return res;
}
 
void init(){
    fac[0] = 1;
    for(int i = 1; i <= N - 10; i++){
        fac[i] = fac[i-1] * i % mod;
    }
    inv[N - 10] = qpow(fac[N - 10], mod - 2);
 
    for(int i = N - 11; i >= 0; i--) {
        inv[i] = inv[i+1] * (i+1) % mod;
    }

}

int C(int x, int y){
    if(x < y || x < 0 || y < 0) return 0;
    return (fac[x] * inv[y]) % mod * inv[x - y] % mod;
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

### 普通高斯消元
```
#define int long long
const int mod = 998244353;

int qpow(int a, int b) {
    int res = 1;
    for(; b; b >>= 1, a = a * a % mod)
        if(b & 1) res = res * a % mod;
    return res;
}

int gauss(vector<vector<int>> a, vector<int> &x) {
    int n = a.size(), m = a[0].size() - 1;
    vector<int> pos(m, -1);
    int rk = 0;

    for(auto &v : a)
        for(auto &t : v) t = (t % mod + mod) % mod;

    for(int j = 0; j < m && rk < n; j++) {
        int p = rk;
        while(p < n && !a[p][j]) p++;
        if(p == n) continue;

        swap(a[p], a[rk]);
        pos[j] = rk;

        int inv = qpow(a[rk][j], mod - 2);
        for(int k = j; k <= m; k++)
            a[rk][k] = a[rk][k] * inv % mod;

        for(int i = 0; i < n; i++) {
            if(i == rk || !a[i][j]) continue;
            int t = a[i][j];
            for(int k = j; k <= m; k++)
                a[i][k] = (a[i][k] - t * a[rk][k] % mod + mod) % mod;
        }
        rk++;
    }

    for(int i = rk; i < n; i++)
        if(a[i][m]) return -1;

    x.assign(m, 0);
    for(int j = 0; j < m; j++)
        if(pos[j] != -1) x[j] = a[pos[j]][m];

    return m - rk;
}
```

### 行列式
```
int det(vector<vector<int>> a) {
    int n = a.size();
    int ans = 1;

    for(auto &v : a)
        for(auto &t : v) t = (t % mod + mod) % mod;

    for(int j = 0; j < n; j++) {
        int p = j;
        while(p < n && !a[p][j]) p++;
        if(p == n) return 0;

        if(p != j) {
            swap(a[p], a[j]);
            ans = (mod - ans) % mod;
        }

        ans = ans * a[j][j] % mod;
        int inv = qpow(a[j][j], mod - 2);

        for(int i = j + 1; i < n; i++) {
            int t = a[i][j] * inv % mod;
            for(int k = j; k < n; k++)
                a[i][k] = (a[i][k] - t * a[j][k] % mod + mod) % mod;
        }
    }
    return ans;
}
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

### 最大流

```cpp
#define For(i, u) for (int i = fir[u]; ~i; i = e[i].nxt)
#define rep(i, a, b) for (int i = (a); i <= (b); i++)
constexpr int inf = 1e9;

const int N = 1e5 + 10;
const int M = 1e6 + 10;

struct Edge {
    int v, c, nxt;
} e[M];

int fir[N], cur[N], d[N];
int cnt = 0;
int n, a[N], b[N];

inline void addedge(int u, int v, int c) {
    e[cnt] = {v, c, fir[u]};
    fir[u] = cnt++;
    
    e[cnt] = {u, 0, fir[v]};
    fir[v] = cnt++;
}

inline bool bfs(int S, int T, int n) {
    rep(i, 0, n) { 
        cur[i] = fir[i];
        d[i] = inf;
    }
    std::queue<int> q;
    q.push(S);
    d[S] = 0;
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        For(i, u) {
            int v = e[i].v;
            if (e[i].c && d[u] + 1 < d[v]) {
                d[v] = d[u] + 1;
                q.push(v);
            }
        }
    }
    return d[T] != inf;
}

inline int dfs(int u, int flow, int T) {
    if (u == T || !flow) return flow;
    int sum = 0;
    
    for (int &i = cur[u]; ~i; i = e[i].nxt) {
        int v = e[i].v;
        if (e[i].c && d[u] + 1 == d[v]) {
            int f = dfs(v, std::min(flow - sum, e[i].c), T);
            e[i].c -= f;
            e[i ^ 1].c += f;
            sum += f;
            if (sum == flow) return sum;
        }
    }
    return sum;
}

inline int dinic(int S, int T, int n) {
    int maxflow = 0;
    while (bfs(S, T, n)) {
        maxflow += dfs(S, inf, T);
    }
    return maxflow;
}

void init(int n) {
    cnt = 0;
    std::fill(fir, fir + n + 1, -1);
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

### dij 费用流(单路增广)

```cpp
#define For(i, u) for(int i = fir[u]; ~i; i = e[i].nx)
#define rep(i, a, b) for(int i = a; i <= b; i++)
#define int long long
#define pii pair<int, int>
constexpr int Inf = 1e18;
int n, m, s, t;

struct edge {
    int v, nx, f, c;
} e[M];

struct node {
    int v, e;
} p[N];

int fir[N], dis[N], h[N], vis[N], cnt;

void spfa()
{
    rep(i, 1, n) h[i] = Inf, vis[i] = false;
    queue<int> q;
    q.push(s);
    h[s] = 0, vis[s] = 1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        vis[u] = 0;
        For(i, u) {
            int v = e[i].v;
            if (e[i].f && h[v] > h[u] + e[i].c) {
                h[v] = h[u] + e[i].c;
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }
    }
}

bool dijkstra()
{
    rep(i, 1, n) dis[i] = Inf, vis[i] = false;
    dis[s] = 0;
    priority_queue<pii, vector<pii>, greater<pii>> q;
    q.push({0, s});

    while (!q.empty()) {
        auto [d, u] = q.top();
        q.pop();
        if (vis[u]) continue;
        vis[u] = true;

        For(i, u) {
            int v = e[i].v;
            if (!e[i].f || h[u] == Inf || h[v] == Inf) continue;

            int nw = e[i].c + h[u] - h[v];
            if (dis[v] > dis[u] + nw) {
                dis[v] = dis[u] + nw;
                p[v] = {u, i};
                q.push({dis[v], v});
            }
        }
    }

    if (dis[t] == Inf) return false;

    rep(i, 1, n)
        if (dis[i] < Inf)
            h[i] += dis[i];

    return true;
}

void addedge(int u, int v, int f, int c) {
    e[cnt] = {v, fir[u], f, c};
    fir[u] = cnt++;
    e[cnt] = {u, fir[v], 0, -c};
    fir[v] = cnt++;
}

pair<int, int> MCMF() {
    spfa();

    int maxf = 0, cost = 0;

    while (dijkstra()) {
        int minf = Inf;

        for (int u = t; u != s; u = p[u].v)
            minf = min(minf, e[p[u].e].f);

        for (int u = t; u != s; u = p[u].v) {
            e[p[u].e].f -= minf;
            e[p[u].e ^ 1].f += minf;
        }

        maxf += minf;
        cost += minf * h[t];
    }

    return {maxf, cost};
}

void init(int n) {
    cnt = 0;
    rep(i, 1, n) fir[i] = -1;
}
```

### dij 费用流(多路增广)

```cpp
#define For(i, u) for(int i = fir[u]; ~i; i = e[i].nx)
#define rep(i, a, b) for(int i = a; i <= b; i++)
#define int long long
#define pii pair<int, int>
constexpr int Inf = 1e18;
int n, m, s, t;

struct edge {
    int v, nx, f, c;
} e[M];

int fir[N], dis[N], h[N], dep[N], cur[N];
int vis[N], cnt;

void spfa()
{
    rep(i, 1, n) h[i] = Inf, vis[i] = false;
    queue<int> q;
    q.push(s);
    h[s] = 0, vis[s] = 1;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        vis[u] = 0;
        For(i, u) {
            int v = e[i].v;
            if (e[i].f && h[v] > h[u] + e[i].c) {
                h[v] = h[u] + e[i].c;
                if (!vis[v]) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }
    }
}

bool dijkstra()
{
    rep(i, 1, n) dis[i] = Inf, vis[i] = false;
    dis[s] = 0;
    priority_queue<pii, vector<pii>, greater<pii>> q;
    q.push({0, s});

    while (!q.empty()) {
        auto [c, u] = q.top();
        q.pop();
        if (vis[u]) continue;
        vis[u] = true;

        For(i, u) {
            int v = e[i].v;
            if (!e[i].f || h[v] == Inf) continue;
            int nw = e[i].c + h[u] - h[v];
            if (dis[v] > dis[u] + nw) {
                dis[v] = dis[u] + nw;
                if (!vis[v])
                    q.push({dis[v], v});
            }
        }
    }

    if (dis[t] == Inf) return false;

    rep(i, 1, n)
        if (dis[i] < Inf)
            h[i] += dis[i];

    return true;
}

bool bfs_mcmf()
{
    rep(i, 1, n) dep[i] = Inf, cur[i] = fir[i];
    queue<int> q;
    q.push(s);
    dep[s] = 0;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        For(i, u) {
            int v = e[i].v;
            if (e[i].f && dep[v] == Inf &&
                h[v] < Inf && e[i].c + h[u] - h[v] == 0) {
                dep[v] = dep[u] + 1;
                q.push(v);
            }
        }
    }

    return dep[t] < Inf;
}

int dfs_mcmf(int u, int flow)
{
    if (u == t || !flow) return flow;

    int sum = 0;

    for (int &i = cur[u]; ~i; i = e[i].nx) {
        int v = e[i].v;

        if (e[i].f && dep[v] == dep[u] + 1 &&
            e[i].c + h[u] - h[v] == 0) {
            int f = dfs_mcmf(v, min(flow - sum, e[i].f));

            e[i].f -= f;
            e[i ^ 1].f += f;
            sum += f;

            if (sum == flow) return sum;
        }
    }

    return sum;
}

void addedge(int u, int v, int f, int c) {
    e[cnt] = {v, fir[u], f, c};
    fir[u] = cnt++;
    e[cnt] = {u, fir[v], 0, -c};
    fir[v] = cnt++;
}

pair<int, int> MCMF() {
    spfa();

    int maxf = 0, cost = 0;

    while (dijkstra()) {
        while (bfs_mcmf()) {
            int f = dfs_mcmf(s, Inf);
            maxf += f;
            cost += f * h[t];
        }
    }
    return {maxf, cost};
}

void init(int n) {
    cnt = 0;
    rep(i, 1, n) fir[i] = -1;
}
```

### 上下界网络流

#### 一、 最大流

> **依赖提示**：需要预先定义好原有的 `e`, `fir`, `cnt` (初始值为 `0`)，以及 `addedge(u, v, c)` 和 `dinic(S, T)` 函数。
>
> 使用上下界网络流时，超级源汇默认取 `SS = n + 1, TT = n + 2`，因此初始化时直接 `init(n + 2)`。

##### 1. 核心变量与建图辅助函数

```cpp
int bal[N];     // 盈亏数组：流入下界和 - 流出下界和
int sum_demand; // 需要补给的总流量

// 初始化 (每次求解前调用)
void init_bound() {
    sum_demand = 0;
    // memset(bal, 0, sizeof(bal)); // 多测时视情况清空
}

// 添加上下界边 (u -> v, 容量 [l, r])
inline void add_bound_edge(int u, int v, int l, int r) {
    bal[v] += l;
    bal[u] -= l;
    addedge(u, v, r - l); // 原图连自由流
}

// 建立超级源汇 SS 到 TT 的补给网络
inline void build_SS_TT() {
    int SS = n + 1, TT = n + 2;
    for (int i = 1; i <= n; i++) {
        if (bal[i] > 0) {
            addedge(SS, i, bal[i]);
            sum_demand += bal[i];
        } else if (bal[i] < 0) {
            addedge(i, TT, -bal[i]);
        }
    }
}
```

##### 2. 无源汇上下界可行流

```cpp
bool solve_feasible_flow() {
    int SS = n + 1, TT = n + 2;

    build_SS_TT();
    n = TT;

    return dinic(SS, TT) == sum_demand; 
}
```

##### 3. 有源汇上下界可行流

```cpp
bool solve_ST_feasible_flow(int S, int T) {
    int SS = n + 1, TT = n + 2;

    addedge(T, S, inf);
    build_SS_TT();

    n = TT;

    return dinic(SS, TT) == sum_demand;
}
```

##### 4. 有源汇上下界最大流 / 最小流

```cpp
// is_max_flow = true (最大流), false (最小流)
void solve_ST_max_min_flow(int S, int T, bool is_max_flow) {
    int SS = n + 1, TT = n + 2;

    int cycle_edge_idx = cnt;
    addedge(T, S, inf); 
    
    build_SS_TT();

    n = TT;
    
    if (dinic(SS, TT) != sum_demand) {
        // 无解处理
        return; 
    }
    
    int ans = e[cycle_edge_idx ^ 1].c; // 当前可行流大小
    e[cycle_edge_idx].c = e[cycle_edge_idx ^ 1].c = 0; // 拆除循环边
    
    if (is_max_flow) ans += dinic(S, T); // 最大流：残量网络正向增广
    else ans -= dinic(T, S);             // 最小流：残量网络反向退流
    
    // ans 即为最终答案
}
```

---

#### 二、 费用流

> **依赖提示**：需要预先定义好原有的 `e`, `fir`, `cnt` (初始值为 `0`)，全局起点终点 `s`, `t`，以及 `addedge(u, v, f, c)`，`spfa()`，`dijkstra()`，`bfs_mcmf()`，`dfs_mcmf()` 函数。
>
> 使用上下界费用流时同样默认 `SS = n + 1, TT = n + 2`，初始化时直接 `init(n + 2)`。

##### 1. 核心变量与建图辅助函数

```cpp
int bal[N];     // 盈亏数组：流入下界和 - 流出下界和
int sum_demand; // 需要补给的总流量
int total_base_cost; // 记录必须流的下界所产生的保底费用

// 初始化 (每次求解前调用)
void init_bound_cost() {
    sum_demand = total_base_cost = 0;
    // memset(bal, 0, sizeof(bal)); 
}

// 添加上下界费用边 (u -> v, 容量 [l, r], 单位费用 w)
inline void add_bound_edge_cost(int u, int v, int l, int r, int w) {
    bal[v] += l;
    bal[u] -= l;
    total_base_cost += l * w; // 累计保底费用
    addedge(u, v, r - l, w);  // 残量网络连自由流
}

// 建立补给网络
inline void build_SS_TT_cost() {
    int SS = n + 1, TT = n + 2;
    for (int i = 1; i <= n; i++) {
        if (bal[i] > 0) {
            addedge(SS, i, bal[i], 0);
            sum_demand += bal[i];
        } else if (bal[i] < 0) {
            addedge(i, TT, -bal[i], 0);
        }
    }
}
```

##### 2. 封装 MCMF 跑流逻辑

```cpp
// 内部函数：势能 + Dijkstra，并在零 reduced-cost 边上跑 blocking flow
// 返回 pair {当前增广的总流量, 当前增广的总费用}
pair<int, int> run_mcmf(int src, int dst) {
    s = src;
    t = dst;

    spfa();

    int flow = 0, cost = 0;

    while (dijkstra()) {
        while (bfs_mcmf()) {
            int f = dfs_mcmf(s, Inf);
            flow += f;
            cost += f * h[t];
        }
    }

    return {flow, cost};
}
```

##### 3. 上下界无源汇最小费用可行流

```cpp
int BDMCMF() {
    int SS = n + 1, TT = n + 2;

    buildST();

    n = TT, s = SS, t = TT;
    
    auto [maxf, cost] = MCMF();
    
    if (maxf != sum_demand) {
        return -1;
    }

    int ans_cost = total_base_cost + cost;
    return ans_cost;
}
```

##### 4. 上下界有源汇最小费用最大流

```cpp
void solve_ST_mincost_maxflow(int S, int T) {
    int SS = n + 1, TT = n + 2;

    int cycle_edge_idx = cnt;
    addedge(T, S, Inf, 0);
    
    build_SS_TT_cost();

    n = TT;
    
    // 第一阶段：填补下界亏空
    auto [maxf, cost1] = run_mcmf(SS, TT);
    
    if (maxf != sum_demand) {
        // 无解处理
        return;
    }
    
    int initial_flow = e[cycle_edge_idx ^ 1].f;
    e[cycle_edge_idx].f = e[cycle_edge_idx ^ 1].f = 0;
    
    // 第二阶段：在残量网络上追求最大流
    auto [add_flow, cost2] = run_mcmf(S, T);
    
    int final_maxflow = initial_flow + add_flow;
    int final_mincost = total_base_cost + cost1 + cost2;
    
    // 输出或返回 final_maxflow 和 final_mincost
}
```

### 计算几何 (dls version)

```
#include <bits/stdc++.h>

#define rep(i, a, b) for(int i = a; i <= b; i++)
#define db double
//#define tp double
#define tp long long
//#define int long long

using namespace std;

const db EPS = 1e-9;

inline int sign(db a) { return a < -EPS ? -1 : (a > EPS); }
inline int cmp(db a, db b) { return sign(a - b); }

inline int sign(int a) { return a < 0 ? -1 : (a > 0); }
inline int cmp(int a, int b) { return (a > b) - (a < b); }

struct P {
    tp x, y;

    P() {}
    P(tp _x, tp _y) : x(_x), y(_y) {}
    
    P operator + (P p) { return {x + p.x, y + p.y}; }
    P operator - (P p) { return {x - p.x, y - p.y}; }
    P operator * (tp d) { return {x * d, y * d}; }
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
    tp abs2() { return x * x + y * y; }
    P rot90() { return P(-y, x); }
    P unit() { return *this / abs(); }
    int quad() const { return sign(y) == 1 || (sign(y) == 0 && sign(x) >= 0); }

    tp dot(P p) { return x * p.x + y * p.y; }
    tp det(P p) { return x * p.y - y * p.x; }

    P rot(db an) { return {x * cos(an) - y * sin(an), x * sin(an) + y * cos(an)}; }

};
tp cross(P p1, P p2) { return p1.x * p2.y - p1.y * p2.x; }
tp cross(P p1, P p2, P p3) { return ((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y)); }
int crossOp(P p1, P p2) { return sign(cross(p1, p2)); } 
int crossOp(P p1, P p2, P p3) { return sign(cross(p1, p2, p3)); }

bool chkLL(P p1, P p2, P q1, P q2) {//直线平行
    db a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return sign(a1 + a2) != 0;
}

P isLL(P p1, P p2, P q1, P q2) { // 直线交点
    tp a1 = cross(q1, q2, p1), a2 = -cross(q1, q2, p2);
    return (p1 * a2 + p2 * a1) / (a1 + a2);
}

bool intersect(tp l1, tp r1, tp l2, tp r2) { // 跨立实验
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

bool isMiddle(tp a, tp m, tp b) {
    return sign(a - m) == 0 || sign(b - m) == 0 || (a < m != b < m);
}

bool isMiddle(P a, P m, P b) {
    return isMiddle(a.x, m.x, b.x) && isMiddle(a.y, m.y, b.y);
}

bool onSeg(P p1, P p2, P q) { // 点在线段上
    return crossOp(p1, p2, q) == 0 && isMiddle(p1, q, p2);
}

bool onSeg_strict(P p1, P p2, P q) {
    return crossOp(p1, p2, q) == 0 && sign((q - p1).dot(q - p2)) < 0;
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

tp area(vector<P> ps) { // 多边形面积 * 2
    tp ret = 0; rep(i, 0, ps.size() - 1) ret += ps[i].det(ps[(i + 1) % ps.size()]);

    return ret;
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

tp dist2(P a, P b) { return (a - b).abs2(); }

db convexDiameter(vector<P> ps) { // 旋转卡尺
    int n = ps.size();
    if (n <= 1) return 0;
    if (n == 2) return ps[0].distTo(ps[1]);
    
    tp maxd2 = 0;
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
struct L {
    P a, b;
    db ang;
    
    L() {}
    L(P _a, P _b) : a(_a), b(_b) { ang = (b - a).alpha(); } 
    
    bool operator < (const L &y) const {
        int c = cmp(ang, y.ang);
        if(c != 0) return c == -1;
        return crossOp(a, y.a, y.b) > 0;
    }
};

/*
struct L {
    P p1, p2; 
    
    P v() const { return p2 - p1; }
    
    bool operator < (const L& b) const {
        P v1 = p2 - p1, v2 = b.p2 - b.p1;
        if (v1.quad() != v2.quad()) 
            return v1.quad() < v2.quad();
        if (sign(v1.det(v2)) != 0) 
            return sign(v1.det(v2)) > 0;
        
        return crossOp(p1, p2, b.p1) < 0; 
    }
};
*/

tp bpmj(vector<L> s) {
    sort(s.begin(), s.end());
    int n = s.size();
    
    vector<L> q(n + 2);
    vector<P> t(n + 2); 
    
    // 现代写法的核心：从 1 开始，彻底抛弃错位的下标
    int l = 1, r = 0;   
    
    for (int i = 0; i < n; ++i) {
        // 极角相同时，由于排序已经把最严的排在最前面了，后续平行的直接忽略
        if (i > 0 && cmp(s[i].ang, s[i - 1].ang) == 0) continue;
        
        // 弹出队尾：如果队尾交点 t[r-1] 被当前半平面 s[i] 切掉了（在右侧或共线）
        while (l < r && crossOp(s[i].a, s[i].b, t[r - 1]) <= 0) r--;
        
        // 弹出队头：如果队头交点 t[l] 被当前半平面 s[i] 切掉了
        while (l < r && crossOp(s[i].a, s[i].b, t[l]) <= 0) l++;
        
        q[++r] = s[i]; // 当前直线入队
        
        // 队列里至少有两条线时，求最新两条线的交点，存入 t[r-1]
        if (l < r) {
            t[r - 1] = isLL(q[r].a, q[r].b, q[r - 1].a, q[r - 1].b); 
        }
    }
    
    // 闭环处理：最后用队头的一刀，再检查并淘汰一次不合格的队尾
    while (l < r && crossOp(q[l].a, q[l].b, t[r - 1]) <= 0) r--;
    
    // 如果有效半平面数量不足以围成多边形（少于3个），面积为 0
    if (r - l <= 1) return 0; 
    
    // 补齐多边形最后首尾相连的那个交点
    t[r] = isLL(q[l].a, q[l].b, q[r].a, q[r].b);
    
    vector<P> poly;
    for (int i = l; i <= r; ++i) {
        poly.push_back(t[i]);
    }
    
    return area(poly); 
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

---

### 集合操作（位运算）

#### 1. 子集枚举

##### 枚举集合 S 的所有子集

时间复杂度：$O(2^{|S|})$，其中 $|S|$ 是集合中元素个数（即二进制中 1 的个数）

```cpp
// 枚举 S 的所有子集（包括空集和 S 本身）
for (int sub = S; ; sub = (sub - 1) & S) {
    // 处理子集 sub
    if (sub == 0) break;
}

// 或者从大到小枚举（包括 S，不包括空集）
for (int sub = S; sub; sub = (sub - 1) & S) {
    // 处理子集 sub
}

// 枚举所有非空子集（不包括空集）
for (int sub = S; sub; sub = (sub - 1) & S) {
    // 处理子集 sub
}
```

**原理**：`(sub - 1) & S` 会得到 S 的下一个子集。减 1 会使最低位的 1 变成 0，并将其右边的所有 0 变成 1，然后与 S 做 AND 运算保证结果仍是 S 的子集。

##### 枚举所有集合的所有子集

```cpp
// 枚举 [0, 2^n) 中每个集合的所有子集
// 总时间复杂度：O(3^n)
for (int S = 0; S < (1 << n); S++) {
    for (int sub = S; sub; sub = (sub - 1) & S) {
        // 处理 S 的子集 sub
    }
}
```

##### 枚举包含某个子集的所有超集

```cpp
// 枚举 [0, 2^n) 中所有包含 S 的集合（S 的所有超集）
for (int sup = S; sup < (1 << n); sup = (sup + 1) | S) {
    // 处理超集 sup
}
```

---

#### 2. 高维前缀和（SOS DP）

**Sum over Subsets Dynamic Programming**

##### 问题

给定函数 $g: 2^n \to \mathbb{R}$，对所有 $S \subseteq \{0, 1, \ldots, n-1\}$，计算：

$$
f(S) = \sum_{T \subseteq S} g(T)
$$

即：对每个集合 $S$，求其所有子集对应的 $g$ 值之和。


##### SOS DP 优化

时间复杂度：$O(n \cdot 2^n)$

**思路**：逐位考虑。定义 $dp[mask][i]$ 表示只考虑 mask 的前 $i$ 位的所有子集的 $g$ 值之和。

**状态转移**：
- 如果 mask 的第 $i$ 位是 0：$dp[mask][i] = dp[mask][i-1]$
- 如果 mask 的第 $i$ 位是 1：$dp[mask][i] = dp[mask][i-1] + dp[mask \oplus 2^i][i-1]$

```cpp
int f[1 << n];

// 初始化
for (int mask = 0; mask < (1 << n); mask++) {
    f[mask] = g[mask];
}

// SOS DP
for (int i = 0; i < n; i++) {
    for (int mask = 0; mask < (1 << n); mask++) {
        if (mask & (1 << i)) {
            f[mask] += f[mask ^ (1 << i)];
        }
    }
}

// 现在 f[mask] = sum of g[submask] for all submask ⊆ mask
```

**注意**：内层循环必须按 mask 从小到大枚举，保证 `mask ^ (1 << i)` 已经被更新。

##### 应用示例

- 给定每个子集的权值，快速查询任意集合的所有子集权值和
- 统计有多少个集合是给定集合的子集
- 计算集合的子集中满足某种性质的个数

---

#### 3. 高维后缀和

##### 问题

给定函数 $g: 2^n \to \mathbb{R}$，对所有 $S \subseteq \{0, 1, \ldots, n-1\}$，计算：

$$
f(S) = \sum_{T \supseteq S} g(T)
$$

即：对每个集合 $S$，求其所有超集对应的 $g$ 值之和。

##### SOS DP 实现

时间复杂度：$O(n \cdot 2^n)$

```cpp
int f[1 << n];

// 初始化
for (int mask = 0; mask < (1 << n); mask++) {
    f[mask] = g[mask];
}

// 高维后缀和
for (int i = 0; i < n; i++) {
    for (int mask = 0; mask < (1 << n); mask++) {
        if (!(mask & (1 << i))) {  // 第 i 位是 0
            f[mask] += f[mask ^ (1 << i)];
        }
    }
}

// 现在 f[mask] = sum of g[supermask] for all supermask ⊇ mask
```

### 三元环计数

```
for(int i=1;i<=m;i++) deg[u[i]]++,deg[v[i]]++;
for(int i=1;i<=m;i++){
	if(make_pair(deg[u[i]],u[i]) > make_pair(deg[v[i]],v[i])) swap(u[i],v[i]);
	g[u[i]].pb(v[i]);
}
int ans=0;
for(int i=1;i<=n;i++){
	for(int j:g[i]) tim[j]=i;
	for(int j:g[i]) for(int k:g[j]) if(tim[k]==i)
        ans++;
}
```

### 四元环计数

```
for(int i=1;i<=m;i++) deg[u[i]]++,deg[v[i]]++;
for(int i=1;i<=m;i++){
	if(make_pair(deg[u[i]],u[i]) > make_pair(deg[v[i]],v[i])) swap(u[i],v[i]);
	g[u[i]].pb(v[i]);rg[v[i]].pb(u[i]);
}
int ans=0;
for(int i=1;i<=n;i++){
	int B1=0,B2=0,B3=0;
	for(int j:g[i]) for(int k:g[j]) if(k^i) (B1+=cnt1[k])%=MOD,cnt1[k]++;
	for(int j:rg[i]) for(int k:g[j]) if(k^i) (B2+=cnt1[k])%=MOD,(B3+=cnt2[k])%=MOD,cnt2[k]++;
	for(int j:g[i]) for(int k:g[j]) cnt1[k]=0;
	for(int j:rg[i]) for(int k:g[j]) cnt2[k]=0;
	ans=(0ll+ans+B1+B2+1ll*INV2*B3)%MOD;
}
```

### 弱周期引理
```
int fail[N], diff[N], top[N];

void build(const string& s) {
    int n = s.length();
    fail[0] = fail[1] = 0;

    // 1. 求基础 fail 数组
    for (int i = 1, j = 0; i < n; i++) {
        while (j > 0 && s[i] != s[j]) {
            j = fail[j];
        }
        if (s[i] == s[j]) {
            j++;
        }
        fail[i + 1] = j;
    }

    // 2. 划分 O(log N) 段等差数列
    for (int i = 1; i <= n; i++) {
        diff[i] = i - fail[i];
        if (fail[i] > 0 && diff[i] == diff[fail[i]]) {
            top[i] = top[fail[i]];
        } else {
            top[i] = fail[i];
        }
    }
}

// 演示：按等差数列“块”快速跳跃
void jump(int len) {
    for (int x = len; x > 0; x = top[x]) {
        // 当前块（等差数列）的公差为 diff[x]
        // 包含的 Border 长度集为: x, x - diff[x], x - 2*diff[x] ... 直到大于 top[x]
        
        // 批量 O(1) 转移 DP 等逻辑写在这里...
        
        /* 展开当前块的代码
        for (int i = x; i > top[x]; i -= diff[x]) {
            // ...
        }
        */
    }
}
```

### fail 链路径压缩
fail[x]：普通 KMP，跳到下一个 border。
fa[x]：在 fail 链上压缩掉一批“下一个字符相同”的 border，直接跳到第一个需要处理的位置。
```
int n;
int fail[N], fa[N];
char s[N];

void solve()
{
    fail[1] = fa[1] = 0;

    for(int i = 2; i <= n; i++) {
        int p = fail[i - 1];

        if(p && s[p + 1] == s[i])
            fa[i - 1] = fa[p];
        else
            fa[i - 1] = p;

        p = fail[i - 1];

        while(p) {
            int q = fa[p];

            if(s[p + 1] != s[i]) {
                while(p != q) {
                    // 处理 border p
                    p = fail[p];
                }
            }

            p = q;
        }

        p = fail[i - 1];

        while(p && s[p + 1] != s[i])
            p = fail[p];

        if(s[p + 1] == s[i]) p++;

        fail[i] = p;
    }
}
```
### ntt
```
#include <bits/stdc++.h>
using namespace std;

const int mod = 998244353, G = 3;

int qpow(int a, int k) {
    int r = 1;
    for (; k; k >>= 1, a = 1LL * a * a % mod)
        if (k & 1) r = 1LL * r * a % mod;
    return r;
}

void ntt(vector<int> &a, bool inv) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int b = n >> 1;
        for (; j & b; b >>= 1) j ^= b;
        j ^= b;
        if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        int w0 = qpow(G, (mod - 1) / len);
        if (inv) w0 = qpow(w0, mod - 2);
        for (int i = 0; i < n; i += len) {
            int w = 1;
            for (int j = 0; j < len / 2; j++) {
                int u = a[i + j], v = 1LL * a[i + j + len / 2] * w % mod;
                a[i + j] = (u + v) % mod;
                a[i + j + len / 2] = (u - v + mod) % mod;
                w = 1LL * w * w0 % mod;
            }
        }
    }
    if (inv) {
        int iv = qpow(n, mod - 2);
        for (int &x : a) x = 1LL * x * iv % mod;
    }
}

vector<int> mul(vector<int> a, vector<int> b) {
    if (a.empty() || b.empty()) return {};
    int need = a.size() + b.size() - 1, n = 1;
    while (n < need) n <<= 1;
    a.resize(n); b.resize(n);
    ntt(a, false); ntt(b, false);
    for (int i = 0; i < n; i++) a[i] = 1LL * a[i] * b[i] % mod;
    ntt(a, true);
    a.resize(need);
    return a;
}

// Returns a(x)^k, retaining only degrees 0..lim-1.
vector<int> poly_pow(vector<int> a, long long k, int lim) {
    if (lim <= 0) return {};
    if ((int)a.size() > lim) a.resize(lim);
    vector<int> res = {1};
    while (k) {
        if (k & 1) {
            res = mul(res, a);
            if ((int)res.size() > lim) res.resize(lim);
        }
        k >>= 1;
        if (k) {
            a = mul(a, a);
            if ((int)a.size() > lim) a.resize(lim);
        }
    }
    res.resize(lim);
    return res;
}


```
