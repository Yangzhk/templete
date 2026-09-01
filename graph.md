# 图论算法模板

## 目录

- [图的存储](#图的存储)
- [最短路](#最短路)
  - [Dijkstra（堆优化）](#dijkstra堆优化)
  - [SPFA / Bellman-Ford](#spfa--bellman-ford)
  - [Floyd-Warshall](#floyd-warshall)
  - [0-1 BFS](#0-1-bfs)
  - [Johnson 全源最短路](#johnson-全源最短路)
  - [K 短路（A\*）](#k-短路a)
- [最小生成树](#最小生成树)
  - [Kruskal](#kruskal)
  - [Prim（堆优化）](#prim堆优化)
  - [Borůvka](#borůvka)
  - [次小生成树](#次小生成树)
  - [朱刘算法（有向最小树形图）](#朱刘算法有向最小树形图)
- [连通性](#连通性)
  - [Tarjan 强连通分量](#tarjan-强连通分量)
  - [割点与桥](#割点与桥)
  - [边双连通分量](#边双连通分量)
  - [点双连通分量与圆方树](#点双连通分量与圆方树)
  - [2-SAT](#2-sat)
- [拓扑排序与 DAG](#拓扑排序与-dag)
- [树上算法](#树上算法)
  - [树哈希](#树哈希)
  - [倍增 LCA](#倍增-lca)
  - [Tarjan 离线 LCA](#tarjan-离线-lca)
  - [树链剖分（HLD）](#树链剖分hld)
  - [DSU on Tree](#dsu-on-tree)
  - [长链剖分](#长链剖分)
  - [点分治与点分树](#点分治与点分树)
  - [虚树](#虚树)
- [网络流](#网络流)
  - [Dinic 最大流](#dinic-最大流)
  - [HLPP 最高标号预流推进](#hlpp-最高标号预流推进)
  - [MCMF 最小费用最大流](#mcmf-最小费用最大流)
  - [上下界网络流](#上下界网络流)
  - [最大权闭合子图](#最大权闭合子图)
- [图匹配](#图匹配)
  - [匈牙利算法](#匈牙利算法)
  - [Hopcroft-Karp](#hopcroft-karp)
  - [KM 算法（二分图最大权完美匹配）](#km-算法二分图最大权完美匹配)
  - [带花树（一般图最大匹配）](#带花树一般图最大匹配)
- [欧拉与哈密顿](#欧拉与哈密顿)
  - [Hierholzer 欧拉路径](#hierholzer-欧拉路径)
  - [状压 Hamilton DP](#状压-hamilton-dp)
- [其他高阶算法](#其他高阶算法)
  - [Stoer-Wagner 全局最小割](#stoer-wagner-全局最小割)
  - [Gomory-Hu 树](#gomory-hu-树)
  - [支配树（Lengauer-Tarjan）](#支配树lengauer-tarjan)
  - [Prüfer 序列](#prüfer-序列)
  - [Kirchhoff 矩阵树定理](#kirchhoff-矩阵树定理)
  - [BEST 定理](#best-定理)

---

## 图的存储

**链式前向星**：节省内存、常用于网络流（成对存边方便 `i ^ 1` 取反向边）。

```cpp
const int N = 1e5 + 10, M = 5e5 + 10;
int head[N], cnt = 1; // cnt 从 1 开始, 这样反向边 i ^ 1 成对
struct Edge { int to, nxt, w; } e[M << 1];

void addEdge(int u, int v, int w = 1) {
    e[++cnt] = {v, head[u], w}; head[u] = cnt;
}
void addUndirected(int u, int v, int w = 1) {
    addEdge(u, v, w); addEdge(v, u, w);
}

// 遍历 u 的出边
for (int i = head[u]; i; i = e[i].nxt) {
    int v = e[i].to, w = e[i].w;
}
```

**vector 邻接表**：写起来直观，常数略大。

```cpp
vector<vector<pair<int,int>>> g; // g[u] = {(v, w)}
g.assign(n + 1, {});
g[u].push_back({v, w});
```

---

## 最短路

### Dijkstra（堆优化）

复杂度 $O((n+m)\log m)$，**只适用于非负边权**。

```cpp
const ll INF = 1e18;
vector<ll> dijkstra(int s, int n, vector<vector<pair<int,int>>>& g) {
    vector<ll> dis(n + 1, INF);
    priority_queue<pair<ll,int>, vector<pair<ll,int>>, greater<>> pq;
    dis[s] = 0; pq.push({0, s});
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dis[u]) continue; // 关键剪枝
        for (auto [v, w] : g[u]) {
            if (dis[u] + w < dis[v]) {
                dis[v] = dis[u] + w;
                pq.push({dis[v], v});
            }
        }
    }
    return dis;
}
```

**应用**：单源最短路；BFS 推广；分层图最短路（拆点：`dis[u][k]` 表示到 $u$ 用了 $k$ 次特殊操作）。

### SPFA / Bellman-Ford

可处理**负权边**；存在负环时无解。SPFA 最坏 $O(nm)$，被卡时改用 Bellman-Ford 或 Johnson。

```cpp
// SPFA + SLF 优化
bool spfa(int s, int n, vector<vector<pair<int,int>>>& g, vector<ll>& dis) {
    dis.assign(n + 1, INF); dis[s] = 0;
    vector<int> cnt(n + 1, 0), inq(n + 1, 0);
    deque<int> q; q.push_back(s); inq[s] = 1;
    while (!q.empty()) {
        int u = q.front(); q.pop_front(); inq[u] = 0;
        for (auto [v, w] : g[u]) {
            if (dis[u] + w < dis[v]) {
                dis[v] = dis[u] + w;
                if (++cnt[v] >= n) return false; // 负环
                if (!inq[v]) {
                    if (!q.empty() && dis[v] < dis[q.front()]) q.push_front(v);
                    else q.push_back(v);
                    inq[v] = 1;
                }
            }
        }
    }
    return true;
}
```

**应用**：差分约束（$x_i - x_j \le c$ 转化为 $j \to i$ 边权 $c$，跑最长/最短路）。

### Floyd-Warshall

$O(n^3)$，全源最短路；处理小图、传递闭包、最小环。

```cpp
for (int k = 1; k <= n; k++)
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
```

**最小环**：枚举 $k$ 时，先用 $1..k-1$ 的最短路加上 $i \to k, k \to j$，再更新 $d[i][j]$。

```cpp
ll ans = INF;
for (int k = 1; k <= n; k++) {
    for (int i = 1; i < k; i++)
        for (int j = i + 1; j < k; j++)
            ans = min(ans, d[i][j] + g[i][k] + g[k][j]);
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n; j++)
            d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
}
```

### 0-1 BFS

边权只有 0/1 时，用 deque 替代堆，复杂度 $O(n + m)$。

```cpp
deque<int> q; q.push_front(s); dis[s] = 0;
while (!q.empty()) {
    int u = q.front(); q.pop_front();
    for (auto [v, w] : g[u]) {
        if (dis[u] + w < dis[v]) {
            dis[v] = dis[u] + w;
            if (w == 0) q.push_front(v);
            else q.push_back(v);
        }
    }
}
```

### Johnson 全源最短路

带负权边时求全源最短路。复杂度 $O(nm \log m)$。

**思路**：新增超级源 $0$ 连向所有点（边权 $0$），跑 Bellman-Ford 得到势函数 $h[u]$；将每条边 $(u,v,w)$ 改为 $w + h[u] - h[v] \ge 0$，再以每点为源跑 Dijkstra；最后还原 $dis[u][v] - h[u] + h[v]$。

### K 短路（A\*）

单起点单终点 K 短路：先反向图跑 Dijkstra 得 $h[u]$（启发函数 = 终点到 $u$ 的最短路），再正向 A\* 用 `f = g + h` 出堆，第 $k$ 次出终点即第 $k$ 短。

```cpp
struct Node { ll g, f; int u; bool operator>(const Node& o) const { return f > o.f; } };
ll kth(int s, int t, int K) {
    auto h = dijkstra_rev(t);
    if (h[s] == INF) return -1;
    priority_queue<Node, vector<Node>, greater<>> pq;
    pq.push({0, h[s], s});
    vector<int> cnt(n + 1, 0);
    while (!pq.empty()) {
        auto [g, f, u] = pq.top(); pq.pop();
        if (++cnt[u] > K) continue;
        if (u == t && cnt[u] == K) return g;
        for (auto [v, w] : G[u])
            if (cnt[v] < K && h[v] != INF)
                pq.push({g + w, g + w + h[v], v});
    }
    return -1;
}
```

更高效的金牌做法：**可持久化左偏树**优化，复杂度 $O((n + m) \log n + K \log K)$。

---

## 最小生成树

### Kruskal

按边权排序 + 并查集，复杂度 $O(m \log m)$。

```cpp
struct DSU {
    vector<int> fa, sz;
    DSU(int n) : fa(n + 1), sz(n + 1, 1) { iota(fa.begin(), fa.end(), 0); }
    int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }
    bool merge(int x, int y) {
        x = find(x); y = find(y);
        if (x == y) return false;
        if (sz[x] < sz[y]) swap(x, y);
        fa[y] = x; sz[x] += sz[y]; return true;
    }
};

ll kruskal(int n, vector<tuple<int,int,int>>& edges) {
    sort(edges.begin(), edges.end(), [](auto& a, auto& b) {
        return get<2>(a) < get<2>(b);
    });
    DSU d(n);
    ll res = 0; int cnt = 0;
    for (auto [u, v, w] : edges) {
        if (d.merge(u, v)) { res += w; if (++cnt == n - 1) break; }
    }
    return cnt == n - 1 ? res : -1;
}
```

### Prim（堆优化）

复杂度 $O((n+m)\log n)$，稠密图比 Kruskal 慢，常用 Kruskal。

### Borůvka

每轮为每个连通块找最短出边并合并，$O(m \log n)$。**适用场景**：边由某种规则隐式生成（如完全图按位运算定义边权），无法显式建图时仍可逐位/逐组找最优出边。

### 次小生成树

先求 MST；对每条非树边 $(u,v,w)$，求树上 $u \to v$ 路径上的最大边权 $\max1$ 和严格次大边权 $\max2$；用 $w - \max1$（或 $w - \max2$ 当 $w = \max1$）更新答案。倍增维护即可。

```cpp
// 倍增维护路径上最大、次大边权
int up[N][LOG], mx1[N][LOG], mx2[N][LOG];
void dfsLCA(int u, int p) {
    for (int k = 1; k < LOG; k++) {
        int m = up[u][k-1];
        up[u][k] = up[m][k-1];
        mx1[u][k] = max(mx1[u][k-1], mx1[m][k-1]);
        mx2[u][k] = max(mx2[u][k-1], mx2[m][k-1]);
        if (mx1[u][k-1] != mx1[m][k-1])
            mx2[u][k] = max(mx2[u][k], min(mx1[u][k-1], mx1[m][k-1]));
    }
    // ...
}
```

### 朱刘算法（有向最小树形图）

固定根 $r$，求所有点到 $r$ 可达的最小有向生成树，$O(VE)$。

```cpp
struct Edge { int u, v; ll w; };
ll zhuLiu(int r, int n, vector<Edge>& es) {
    ll res = 0;
    while (true) {
        vector<ll> in(n + 1, INF);
        vector<int> pre(n + 1, -1);
        for (auto& e : es) if (e.u != e.v && e.w < in[e.v])
            in[e.v] = e.w, pre[e.v] = e.u;
        in[r] = 0;
        for (int i = 1; i <= n; i++) if (in[i] == INF) return -1;
        // 找环
        vector<int> id(n + 1, -1), vis(n + 1, -1);
        int cnt = 0;
        for (int i = 1; i <= n; i++) {
            res += in[i];
            int v = i;
            while (vis[v] != i && id[v] == -1 && v != r)
                vis[v] = i, v = pre[v];
            if (v != r && id[v] == -1) {
                ++cnt;
                for (int u = pre[v]; u != v; u = pre[u]) id[u] = cnt;
                id[v] = cnt;
            }
        }
        if (cnt == 0) return res; // 无环则结束
        for (int i = 1; i <= n; i++) if (id[i] == -1) id[i] = ++cnt;
        // 缩点
        for (auto& e : es) {
            ll w = e.w - in[e.v];
            e.u = id[e.u]; e.v = id[e.v]; e.w = w;
        }
        n = cnt; r = id[r];
    }
}
```

**应用**：每个节点选一条入边的最小代价（如某些"分配问题"）。

---

## 连通性

### Tarjan 强连通分量

$O(n + m)$ 求 SCC。缩点后得到 DAG。

```cpp
int dfn[N], low[N], dfsClock, scc[N], sccCnt;
int stk[N], top, inStk[N];
vector<int> g[N];

void tarjan(int u) {
    dfn[u] = low[u] = ++dfsClock;
    stk[++top] = u; inStk[u] = 1;
    for (int v : g[u]) {
        if (!dfn[v]) { tarjan(v); low[u] = min(low[u], low[v]); }
        else if (inStk[v]) low[u] = min(low[u], dfn[v]);
    }
    if (low[u] == dfn[u]) {
        ++sccCnt;
        int v;
        do { v = stk[top--]; inStk[v] = 0; scc[v] = sccCnt; } while (v != u);
    }
}
```

**应用**：缩点后 DAG 上 DP；判可达性；2-SAT 建图后判可行性。

### 割点与桥

无向图。`low[v] >= dfn[u]` ⇒ $u$ 是割点（根需出度 ≥ 2）；`low[v] > dfn[u]` ⇒ $(u,v)$ 是桥。

```cpp
// 割点
void tarjanCut(int u, int root) {
    dfn[u] = low[u] = ++dfsClock;
    int child = 0;
    for (int v : g[u]) {
        if (!dfn[v]) {
            tarjanCut(v, root); ++child;
            low[u] = min(low[u], low[v]);
            if (u != root && low[v] >= dfn[u]) cut[u] = 1;
        } else low[u] = min(low[u], dfn[v]);
    }
    if (u == root && child >= 2) cut[u] = 1;
}

// 桥（注意重边：传入边的编号 ein 而非父节点）
void tarjanBridge(int u, int ein) {
    dfn[u] = low[u] = ++dfsClock;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (!dfn[v]) {
            tarjanBridge(v, i);
            low[u] = min(low[u], low[v]);
            if (low[v] > dfn[u]) bridge[i] = bridge[i ^ 1] = 1;
        } else if (i != (ein ^ 1)) low[u] = min(low[u], dfn[v]);
    }
}
```

### 边双连通分量

去掉所有桥后的连通分量。可在 Tarjan 求桥时同时维护栈：

```cpp
void tarjanEBCC(int u, int ein) {
    dfn[u] = low[u] = ++dfsClock;
    stk[++top] = u;
    for (int i = head[u]; i; i = e[i].nxt) {
        int v = e[i].to;
        if (!dfn[v]) {
            tarjanEBCC(v, i);
            low[u] = min(low[u], low[v]);
        } else if (i != (ein ^ 1)) low[u] = min(low[u], dfn[v]);
    }
    if (low[u] == dfn[u]) {
        ++ebccCnt; int v;
        do { v = stk[top--]; ebcc[v] = ebccCnt; } while (v != u);
    }
}
```

### 点双连通分量与圆方树

点双：极大不含割点的子图（割点可同属多个点双）。

**圆方树**：原图每个点为"圆点"，每个点双新建一个"方点"，方点连接该点双的所有原点。圆方树是树，且任意两点 $u \to v$ 路径上的方点恰好是 $u \to v$ 所有简单路径并集中的点双。

```cpp
int dfn[N], low[N], dfsClock, top, stk[N], idx; // idx = n + 已建方点数
vector<int> T[N << 1]; // 圆方树

void tarjanBCC(int u) {
    dfn[u] = low[u] = ++dfsClock;
    stk[++top] = u;
    for (int v : g[u]) {
        if (!dfn[v]) {
            tarjanBCC(v);
            low[u] = min(low[u], low[v]);
            if (low[v] >= dfn[u]) {
                ++idx; int w;
                do {
                    w = stk[top--];
                    T[idx].push_back(w); T[w].push_back(idx);
                } while (w != v);
                T[idx].push_back(u); T[u].push_back(idx);
            }
        } else low[u] = min(low[u], dfn[v]);
    }
}
```

**应用**：仙人掌图问题；"必经点"统计；广义圆方树上跑树形 DP。

### 2-SAT

每个变量 $x$ 拆成 $x$ 和 $\neg x$ 两点；子句 $a \lor b$ 等价于 $\neg a \to b$ 与 $\neg b \to a$；跑 SCC，若 $x$ 与 $\neg x$ 同 SCC 则无解；否则按 SCC 拓扑序逆序赋值（编号大的为真）。

```cpp
// 变量 i 的两个点: 2*i, 2*i+1 (真/假)
void addClause(int x, bool xv, int y, bool yv) {
    // x = xv 或 y = yv
    addEdge(2 * x + !xv, 2 * y + yv); // !x -> y
    addEdge(2 * y + !yv, 2 * x + xv);
}
// 跑 Tarjan SCC; 答案: x 取 (scc[2*x] > scc[2*x+1])
```

---

## 拓扑排序与 DAG

```cpp
vector<int> topo(int n, vector<vector<int>>& g) {
    vector<int> in(n + 1), res;
    for (int u = 1; u <= n; u++) for (int v : g[u]) in[v]++;
    queue<int> q;
    for (int i = 1; i <= n; i++) if (!in[i]) q.push(i);
    while (!q.empty()) {
        int u = q.front(); q.pop(); res.push_back(u);
        for (int v : g[u]) if (--in[v] == 0) q.push(v);
    }
    return (int)res.size() == n ? res : vector<int>{};
}
```

**应用**：DAG 最长路（DP）；字典序最小拓扑序换 priority_queue；判环。

---

## 树上算法

### 树哈希

```
#include <cctype>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using ull = unsigned long long;

const ull mask = std::mt19937_64(time(nullptr))();

ull shift(ull x) {
  x ^= mask;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  x ^= mask;
  return x;
}

constexpr int N = 1e6 + 10;

int n;
ull hash[N];
std::vector<int> edge[N];
std::set<ull> trees;

void getHash(int x, int p) {
  hash[x] = 1;
  for (int i : edge[x]) {
    if (i == p) {
      continue;
    }
    getHash(i, x);
    hash[x] += shift(hash[i]);
  }
  trees.insert(hash[x]);
}

using std::cin;
using std::cout;

int main() {
  cin.tie(nullptr)->sync_with_stdio(false);
  cin >> n;
  for (int i = 1; i < n; i++) {
    int u, v;
    cin >> u >> v;
    edge[u].push_back(v);
    edge[v].push_back(u);
  }
  getHash(1, 0);
  cout << trees.size();
}
```

### 倍增 LCA

$O((n+q)\log n)$ 预处理 + 查询。

```cpp
const int LOG = 20;
int dep[N], up[N][LOG];

void dfs(int u, int p) {
    up[u][0] = p; dep[u] = dep[p] + 1;
    for (int k = 1; k < LOG; k++) up[u][k] = up[up[u][k-1]][k-1];
    for (int v : g[u]) if (v != p) dfs(v, u);
}

int lca(int u, int v) {
    if (dep[u] < dep[v]) swap(u, v);
    int d = dep[u] - dep[v];
    for (int k = 0; k < LOG; k++) if (d >> k & 1) u = up[u][k];
    if (u == v) return u;
    for (int k = LOG - 1; k >= 0; k--)
        if (up[u][k] != up[v][k]) u = up[u][k], v = up[v][k];
    return up[u][0];
}
```

### Tarjan 离线 LCA

$O((n+q)\alpha)$，把所有询问按端点挂载，DFS 回溯时合并并查集。

```cpp
int fa[N], vis[N], ans[Q];
vector<pair<int,int>> qry[N]; // {另一端, 询问编号}

int find(int x) { return fa[x] == x ? x : fa[x] = find(fa[x]); }

void tarjanLCA(int u, int p) {
    fa[u] = u;
    for (int v : g[u]) if (v != p) { tarjanLCA(v, u); fa[v] = u; }
    vis[u] = 1;
    for (auto [v, id] : qry[u]) if (vis[v]) ans[id] = find(v);
}
```

### 树链剖分（HLD）

按重儿子剖分，每条重链 DFS 序连续，可用线段树维护路径/子树信息。$O(\log^2 n)$ 单次查询/修改。

```cpp
int sz[N], hson[N], top[N], dfn[N], rk[N], dep[N], fa[N], dfsClock;

void dfs1(int u, int p) {
    fa[u] = p; dep[u] = dep[p] + 1; sz[u] = 1; hson[u] = 0;
    for (int v : g[u]) if (v != p) {
        dfs1(v, u);
        sz[u] += sz[v];
        if (sz[v] > sz[hson[u]]) hson[u] = v;
    }
}

void dfs2(int u, int t) {
    top[u] = t; dfn[u] = ++dfsClock; rk[dfsClock] = u;
    if (hson[u]) dfs2(hson[u], t);
    for (int v : g[u]) if (v != fa[u] && v != hson[u]) dfs2(v, v);
}

// 路径 u-v 操作 (例如线段树区间加)
void updatePath(int u, int v, ll w) {
    while (top[u] != top[v]) {
        if (dep[top[u]] < dep[top[v]]) swap(u, v);
        seg.update(dfn[top[u]], dfn[u], w);
        u = fa[top[u]];
    }
    if (dep[u] > dep[v]) swap(u, v);
    seg.update(dfn[u], dfn[v], w);
}
// 子树 u 操作: seg.update(dfn[u], dfn[u] + sz[u] - 1, w);
```

**应用**：路径求和/最值/染色；子树询问；LCA（求路径时顺便得到）。

### DSU on Tree

子树询问、不带修：$O(n \log n)$。轻儿子算完清空，重儿子保留信息。

```cpp
int sz[N], hson[N], cnt[VAL], cur;
ll ans[N];

void add(int u, int p, int x) {
    // 把 u 子树（除被 keep 的部分外）的贡献加入/删除
}

void dfs(int u, int p, bool keep) {
    for (int v : g[u]) if (v != p && v != hson[u]) dfs(v, u, false);
    if (hson[u]) dfs(hson[u], u, true);
    // 加入 u 自己 + 所有轻儿子子树
    ans[u] = cur;
    if (!keep) /* 清空 u 子树贡献 */;
}
```

**应用**：子树众数、子树颜色种数等离线问题。

### 长链剖分

按"最深链"剖分。$O(n)$ 解决"与子树最深距离相关"的 DP（如：每个点的 $f[u][k]$ = 子树内距 $u$ 为 $k$ 的点数），用指针在数组上偏移传给重儿子。

**应用**：$O(n \log n)$ 求 K 级祖先；树上以深度为下标的 DP。

### 点分治与点分树

**点分治**：每次找重心作根，统计经过该重心的路径，再递归子树。$O(n \log n)$ 路径统计。

```cpp
int sz[N], mx[N], rt, totSz; bool vis[N];

void getSize(int u, int p) {
    sz[u] = 1; mx[u] = 0;
    for (auto [v, w] : g[u]) if (v != p && !vis[v]) {
        getSize(v, u); sz[u] += sz[v];
        mx[u] = max(mx[u], sz[v]);
    }
}
void getRoot(int u, int p) {
    int cur = max(mx[u], totSz - sz[u]);
    if (cur < mx[rt]) rt = u;
    for (auto [v, w] : g[u]) if (v != p && !vis[v]) getRoot(v, u);
}
void solve(int u) {
    vis[u] = 1;
    // ... 统计经过 u 的路径
    for (auto [v, w] : g[u]) if (!vis[v]) {
        totSz = sz[v]; mx[rt = 0] = INT_MAX;
        getSize(v, u); getRoot(v, u);
        solve(rt);
    }
}
```

**点分树**：把每层重心连成一棵树（深度 $O(\log n)$），支持带修路径/距离询问。

### 虚树

给定关键点集合 $S$，构造仅包含 $S$ 与其两两 LCA 的"虚树"，规模 $O(|S|)$。

```cpp
bool cmpDfn(int x, int y) { return dfn[x] < dfn[y]; }

vector<vector<int>> vt;
int stk[N], top;

void buildVT(vector<int>& key) {
    sort(key.begin(), key.end(), cmpDfn);
    auto add = [&](int u, int v) { vt[u].push_back(v); };
    stk[top = 1] = 1; vt[1].clear();
    for (int u : key) {
        if (u == 1) continue;
        int l = lca(u, stk[top]);
        if (l != stk[top]) {
            while (top > 1 && dfn[l] < dfn[stk[top - 1]]) {
                add(stk[top - 1], stk[top]); --top;
            }
            if (l != stk[top - 1]) { vt[l].clear(); add(l, stk[top]); stk[top] = l; }
            else { add(l, stk[top--]); }
        }
        vt[u].clear(); stk[++top] = u;
    }
    while (top > 1) { add(stk[top - 1], stk[top]); --top; }
}
```

**应用**：多组询问每次仅涉及若干关键点的树形 DP（如"覆盖关键点的最小代价"）。

---

## 网络流

### Dinic 最大流

复杂度 $O(V^2 E)$，单位容量图 $O(E\sqrt V)$，二分图 $O(E\sqrt V)$。

```cpp
struct Dinic {
    struct E { int to, nxt; ll c; };
    vector<E> e; vector<int> head, cur, d;
    int n, S, T;
    void init(int n_) { n = n_; head.assign(n + 1, -1); e.clear(); }
    void addEdge(int u, int v, ll c) {
        e.push_back({v, head[u], c}); head[u] = e.size() - 1;
        e.push_back({u, head[v], 0}); head[v] = e.size() - 1;
    }
    bool bfs() {
        d.assign(n + 1, -1); d[S] = 0;
        queue<int> q; q.push(S);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int i = head[u]; ~i; i = e[i].nxt)
                if (e[i].c && d[e[i].to] < 0) {
                    d[e[i].to] = d[u] + 1; q.push(e[i].to);
                }
        }
        return d[T] >= 0;
    }
    ll dfs(int u, ll f) {
        if (u == T || !f) return f;
        ll res = 0;
        for (int& i = cur[u]; ~i; i = e[i].nxt) {
            int v = e[i].to;
            if (e[i].c && d[v] == d[u] + 1) {
                ll x = dfs(v, min(f, e[i].c));
                e[i].c -= x; e[i ^ 1].c += x;
                res += x; f -= x;
                if (!f) break;
            }
        }
        if (!res) d[u] = -1;
        return res;
    }
    ll maxflow(int s, int t) {
        S = s; T = t; ll flow = 0;
        while (bfs()) { cur = head; flow += dfs(S, LLONG_MAX); }
        return flow;
    }
};
```

**应用**：二分图匹配、最小割、项目选择、棋盘覆盖、区间调度。

### HLPP 最高标号预流推进

复杂度 $O(V^2 \sqrt E)$，理论与实际都比 Dinic 更快，常用于大规模流图（金牌常用）。

```cpp
// 核心：维护高度 h[u]，按"高度桶 + gap 优化"循环 push & relabel
// 实现较长，比赛中通常用 Dinic，遇到 1e5 节点的稠密流时再考虑 HLPP
```

### MCMF 最小费用最大流

SPFA 求最短增广路 + Dinic 多路增广。

```cpp
struct MCMF {
    struct E { int to, nxt; ll c, w; };
    vector<E> e; vector<int> head;
    vector<ll> dis; vector<int> vis, cur;
    int n, S, T; ll totFlow, totCost;
    void init(int n_) { n = n_; head.assign(n + 1, -1); e.clear(); totFlow = totCost = 0; }
    void addEdge(int u, int v, ll c, ll w) {
        e.push_back({v, head[u], c,  w}); head[u] = e.size() - 1;
        e.push_back({u, head[v], 0, -w}); head[v] = e.size() - 1;
    }
    bool spfa() {
        dis.assign(n + 1, LLONG_MAX); vis.assign(n + 1, 0);
        queue<int> q; dis[S] = 0; q.push(S); vis[S] = 1;
        while (!q.empty()) {
            int u = q.front(); q.pop(); vis[u] = 0;
            for (int i = head[u]; ~i; i = e[i].nxt) {
                int v = e[i].to;
                if (e[i].c && dis[u] + e[i].w < dis[v]) {
                    dis[v] = dis[u] + e[i].w;
                    if (!vis[v]) { q.push(v); vis[v] = 1; }
                }
            }
        }
        return dis[T] != LLONG_MAX;
    }
    ll dfs(int u, ll f) {
        if (u == T) return f;
        vis[u] = 1;
        ll res = 0;
        for (int& i = cur[u]; ~i; i = e[i].nxt) {
            int v = e[i].to;
            if (!vis[v] && e[i].c && dis[u] + e[i].w == dis[v]) {
                ll x = dfs(v, min(f, e[i].c));
                if (x) {
                    e[i].c -= x; e[i ^ 1].c += x;
                    totCost += x * e[i].w;
                    res += x; f -= x;
                    if (!f) break;
                }
            }
        }
        vis[u] = 0;
        return res;
    }
    pair<ll,ll> solve(int s, int t) {
        S = s; T = t;
        while (spfa()) { cur = head; totFlow += dfs(S, LLONG_MAX); }
        return {totFlow, totCost};
    }
};
```

**优化**：Primal-Dual（势函数 + Dijkstra 替代 SPFA）能避免负权 SPFA 退化。

### 上下界网络流

边 $(u,v)$ 容量在 $[L,R]$。**无源汇可行流**：建辅助源汇 $S',T'$；将 $L$ 视为"必须流"，建 $u \to v$ 容量 $R-L$ 的边，并对 $u$ 减 $L$、$v$ 加 $L$ 流量需求；最后从 $S'$ 向"需求流入"的点连边、"流出"的点连向 $T'$；满流则可行。

**有源汇最大流**：$T \to S$ 加无穷边转化为无源汇；可行后再在原 $S \to T$ 上跑 Dinic。

### 最大权闭合子图

带点权图，选一些点要求若选 $u$ 则必须选 $u$ 的所有"后继"。

**建图**：源 $S$ 向正权点连容量 $w_u$ 的边，负权点向汇 $T$ 连容量 $-w_u$ 的边，原图 $u \to v$ 连无穷边。**答案 = 正权和 - 最小割**。

---

## 图匹配

### 匈牙利算法

二分图最大匹配，$O(VE)$。

```cpp
int matchL[N], matchR[N], vis[N];
vector<int> g[N];

bool dfs(int u) {
    for (int v : g[u]) if (!vis[v]) {
        vis[v] = 1;
        if (!matchR[v] || dfs(matchR[v])) {
            matchL[u] = v; matchR[v] = u;
            return true;
        }
    }
    return false;
}

int hungary(int n) {
    int res = 0;
    for (int u = 1; u <= n; u++) {
        fill(vis + 1, vis + n + 1, 0);
        if (dfs(u)) ++res;
    }
    return res;
}
```

### Hopcroft-Karp

二分图匹配 $O(E\sqrt V)$，等价于在二分图上跑 Dinic。

### KM 算法（二分图最大权完美匹配）

$O(n^3)$。维护可行顶标 $l_x, l_y$（$l_x[i] + l_y[j] \ge w[i][j]$），用 slack 与 BFS 加速。

```cpp
const ll INF = 1e18;
int n; ll w[N][N], lx[N], ly[N], slack[N];
int matchY[N], pre[N]; bool visx[N], visy[N];

void bfsKM(int u) {
    fill(slack + 1, slack + n + 1, INF);
    fill(visx + 1, visx + n + 1, 0);
    fill(visy + 1, visy + n + 1, 0);
    int y = 0; matchY[0] = u;
    do {
        visy[y] = 1; int x = matchY[y], ny = 0; ll d = INF;
        for (int j = 1; j <= n; j++) if (!visy[j]) {
            ll t = lx[x] + ly[j] - w[x][j];
            if (t < slack[j]) { slack[j] = t; pre[j] = y; }
            if (slack[j] < d) { d = slack[j]; ny = j; }
        }
        for (int j = 0; j <= n; j++) {
            if (visy[j]) lx[matchY[j]] -= d, ly[j] += d;
            else slack[j] -= d;
        }
        y = ny;
    } while (matchY[y]);
    while (y) { matchY[y] = matchY[pre[y]]; y = pre[y]; }
}

ll km() {
    fill(matchY, matchY + n + 1, 0);
    fill(lx, lx + n + 1, 0); fill(ly, ly + n + 1, 0);
    for (int i = 1; i <= n; i++) bfsKM(i);
    ll res = 0;
    for (int j = 1; j <= n; j++) res += w[matchY[j]][j];
    return res;
}
```

**应用**：分配问题（n 个工人 n 个任务，最大化收益）。**注意**：求最小权时取负或用 INF - w。

### 带花树（一般图最大匹配）

非二分图最大匹配，$O(VE)$ 或 $O(V^3)$。**核心**：BFS 时遇到奇环（"花"），将整条花路径上的点缩成一个点继续 BFS。

实现长（200+ 行），仅在题目明确"非二分图最大匹配"时使用，否则尽量转二分图。

---

## 欧拉与哈密顿

### Hierholzer 欧拉路径

无向图存在欧拉路径：连通且度数为奇数的点 $\le 2$；欧拉回路：度数全偶。
有向图欧拉路径：弱连通且最多一个点 $|in - out| = 1$ 起点、一个 $|out - in| = 1$ 终点。

```cpp
// 无向图，链式前向星 cnt 从 2 开始
vector<int> path; int cur[N];
void dfs(int u) {
    for (int& i = cur[u]; i; ) {
        int v = e[i].to;
        if (vis[i / 2]) { i = e[i].nxt; continue; }
        vis[i / 2] = 1; i = e[i].nxt;
        dfs(v);
    }
    path.push_back(u);
}
// 最后 reverse(path)，即欧拉路径
```

**应用**：一笔画问题、Eulerian 序列拼接。

### 状压 Hamilton DP

$f[S][i]$ = 已访问集合 $S$、当前在 $i$ 的最小代价；$O(2^n n^2)$。

```cpp
vector<vector<ll>> f(1 << n, vector<ll>(n, INF));
f[1][0] = 0;
for (int S = 1; S < (1 << n); S++)
    for (int i = 0; i < n; i++) if ((S >> i & 1) && f[S][i] < INF)
        for (int j = 0; j < n; j++) if (!(S >> j & 1))
            f[S | (1 << j)][j] = min(f[S | (1 << j)][j], f[S][i] + w[i][j]);
ll ans = INF;
for (int i = 1; i < n; i++) ans = min(ans, f[(1 << n) - 1][i] + w[i][0]); // 哈回路
```

---

## 其他高阶算法

### Stoer-Wagner 全局最小割

无向带权图，求两点不指定的最小割，$O(V^3)$ 或 $O(VE + V^2 \log V)$。

**思路**：每轮做"最大邻接序"找出最后两点 $s,t$，记录 $cut(t)$ 即 $s$-$t$ 的割值，候选答案；将 $s,t$ 合并继续。

```cpp
int n; ll w[N][N], dis[N]; bool vis[N], merged[N];
ll stoerWagner() {
    ll res = LLONG_MAX;
    for (int i = 0; i < n - 1; i++) {
        memset(vis, 0, sizeof vis);
        memset(dis, 0, sizeof dis);
        int prev = -1, last = 0;
        for (int j = 0; j < n - i; j++) {
            int sel = -1;
            for (int k = 0; k < n; k++)
                if (!merged[k] && !vis[k] && (sel == -1 || dis[k] > dis[sel])) sel = k;
            if (j == n - i - 1) {
                res = min(res, dis[sel]);
                merged[sel] = 1;
                for (int k = 0; k < n; k++) if (!merged[k])
                    w[prev][k] += w[sel][k], w[k][prev] += w[k][sel];
                break;
            }
            vis[sel] = 1; prev = sel;
            for (int k = 0; k < n; k++) if (!merged[k] && !vis[k]) dis[k] += w[sel][k];
        }
    }
    return res;
}
```

### Gomory-Hu 树

无向图，构造一棵 $n-1$ 条边的等价"最小割树"：树上 $u,v$ 路径最小边权 = 原图 $u,v$ 最小割。$n-1$ 次最大流即可构造。

```cpp
int parent[N]; ll cutVal[N];
void buildGH(int n) {
    for (int i = 1; i <= n; i++) parent[i] = 1;
    for (int s = 2; s <= n; s++) {
        // 复原网络容量
        ll f = maxflow(s, parent[s]);
        cutVal[s] = f;
        // BFS 标记 s 一侧
        for (int i = s + 1; i <= n; i++)
            if (sideOfSrc(i) && parent[i] == parent[s]) parent[i] = s;
    }
}
```

**应用**：多组 $u,v$ 询问最小割，建树后路径最小值即可。

### 支配树（Lengauer-Tarjan）

有向图固定起点 $r$，$u$ 的"支配点" = 所有从 $r$ 到 $u$ 的路径必经的点；支配树是这些支配关系构成的树。

复杂度 $O((n+m)\alpha(n))$ 或 $O((n+m)\log n)$。常用于编译器分析，竞赛中较少出现，遇到即"金牌题"。

### Prüfer 序列

$n$ 个标号节点的有标号无根树 ↔ 长度 $n-2$ 序列（每位 $\in [1,n]$）。

**编码**：每次找编号最小的叶子，写出其父节点编号，删除该叶子；重复 $n-2$ 次。
**结论**：节点 $i$ 在序列出现次数 = $\deg(i) - 1$；标号树数 = $n^{n-2}$（Cayley 公式）。

### Kirchhoff 矩阵树定理

无向图生成树个数 = 度数矩阵 - 邻接矩阵（Laplacian）的任意一个余子式的行列式。

```cpp
// L[i][i] = deg(i); L[i][j] = -(i,j 的边数)
// 删去任意一行一列后取行列式 (mod 意义下用 Gauss-Jordan)
```

**应用**：生成树计数；带权矩阵树计数所有生成树边权乘积之和（把 $-1$ 改为 $-w_{ij}$）。

### BEST 定理

有向欧拉图的欧拉回路数：

$$
ec(G) = t_w(G) \cdot \prod_{v \in V}(\deg^+(v) - 1)!
$$

其中 $t_w(G)$ 为以任意点 $w$ 为根的内向（或外向）生成树数（用矩阵树求）。

**应用**：De Bruijn 序列计数、特殊有向图欧拉回路计数。

整个算法流程分三步：

1. **特判**：检查所有节点的入度是否等于出度。如果有一个不相等，欧拉回路数量直接为 0。
2. **构建基尔霍夫矩阵（Laplacian Matrix）**：
    * 主对角线 $L[i][i]$ 填节点 $i$ 的出度。
    * 非主对角线 $L[i][j]$ 填节点 $i$ 到 $j$ 的有向边数的相反数（即 $-cnt(i \to j)$）。
3. **求行列式**：随便挑一个度数不为 0 的节点 $w$（比如节点 1），把矩阵的第 $w$ 行和第 $w$ 列删掉。对剩下大小为 $(n-1) \times (n-1)$ 的矩阵做高斯消元求行列式，得到的就是 $T^{root}_w(G)$。
4. **注意**: 求欧拉回路的时候, 一定要手动删掉 0 度点, 否则算出的行列式不对.
```
#define int long long // 习惯全开 long long 防止乘法溢出
const int MOD = 998244353;
const int MAXN = 505;

int in_deg[MAXN], out_deg[MAXN];
int L[MAXN][MAXN]; // 基尔霍夫矩阵
int fact[MAXN * MAXN]; // 预处理阶乘

// 快速幂求逆元（如果需要）
int qpow(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1) res = res * a % MOD;
        a = a * a % MOD;
        b >>= 1;
    }
    return res;
}

// 高斯消元求行列式 (大小为 n 的方阵)
int det(int n) {
    int res = 1;
    for (int i = 1; i <= n; i++) {
        // 找主元
        int pivot = i;
        for (int j = i + 1; j <= n; j++) {
            if (L[j][i]) {
                pivot = j;
                break;
            }
        }
        if (pivot != i) {
            swap(L[i], L[pivot]);
            res = (MOD - res) % MOD; // 交换两行，行列式变号
        }
        if (!L[i][i]) return 0; // 满秩判断

        // 消元
        int inv = qpow(L[i][i], MOD - 2);
        for (int j = i + 1; j <= n; j++) {
            if (L[j][i]) {
                int rate = L[j][i] * inv % MOD;
                for (int k = i; k <= n; k++) {
                    L[j][k] = (L[j][k] - rate * L[i][k]) % MOD;
                    if (L[j][k] < 0) L[j][k] += MOD;
                }
            }
        }
        res = res * L[i][i] % MOD;
    }
    return (res % MOD + MOD) % MOD;
}

signed main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // 1. 预处理阶乘
    fact[0] = 1;
    for (int i = 1; i < MAXN * MAXN; i++) {
        fact[i] = fact[i - 1] * i % MOD;
    }

    int n, m; // n 个点, m 条边
    if (!(cin >> n >> m)) return 0;

    for (int i = 1; i <= m; i++) {
        int u, v;
        cin >> u >> v;
        out_deg[u]++;
        in_deg[v]++;
        
        // 构建基尔霍夫矩阵
        L[u][u]++;             // 对角线加上出度
        L[u][v]--;             // 有向边 u->v
    }

    // 2. 检查欧拉图基本条件
    for (int i = 1; i <= n; i++) {
        if (in_deg[i] != out_deg[i]) {
            cout << 0 << "\n";
            return 0;
        }
    }

    // 3. 处理矩阵 (由于可能有负数，统一步长取模)
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            L[i][j] = (L[i][j] % MOD + MOD) % MOD;
        }
    }

    // 4. 去掉第一行第一列 (假设节点1在连通块内)，求 n-1 阶行列式
    // 注意：实现时直接将 det 函数传入 n-1，且矩阵从下标 2 开始传入或者直接覆盖平移。
    // 为了方便，这里把第2行到第n行，平移到1到n-1
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < n; j++) {
            L[i][j] = L[i + 1][j + 1];
        }
    }
    
    int tree_cnt = det(n - 1); // T_w 的数量

    // 5. 乘上阶乘部分
    int ans = tree_cnt;
    for (int i = 1; i <= n; i++) {
        if (out_deg[i] > 0) {
            ans = ans * fact[out_deg[i] - 1] % MOD;
        }
    }

    cout << ans << "\n";
    return 0;
}
```
