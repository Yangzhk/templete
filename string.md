<img width="2440" height="1244" alt="image" src="https://github.com/user-attachments/assets/dc37b5eb-5513-45c9-9bb2-1272afd2e0ee" /># 字符串算法模板

## 目录

- [字符串算法模板](#字符串算法模板)
  - [目录](#目录)
  - [1. 字符串哈希](#1-字符串哈希)
  - [2. KMP — 前缀函数 / 模式匹配](#2-kmp--前缀函数--模式匹配)
  - [3. Z 函数](#3-z-函数)
  - [4. 字典树](#4-字典树)
  - [5. AC 自动机](#5-ac-自动机)
  - [6. Manacher 算法 — 回文](#6-manacher-算法--回文)
  - [7. 后缀数组 - O(n log n)](#7-后缀数组---on-log-n)
  - [压缩后缀树](#压缩后缀树)
  - [8. 后缀自动机 (SAM) — O(nA)](#8-后缀自动机-sam--ona)
  - [9. 最长公共子序列 (LCS) — DP](#9-最长公共子序列-lcs--dp)
  - [10. 最小循环表示 (Booth 算法)](#10-最小循环表示-booth-算法)
  - [11. Lyndon 分解 (Duval 算法)](#11-lyndon-分解-duval-算法)
  - [12. 编辑距离 (Levenshtein)](#12-编辑距离-levenshtein)
  - [13. 回文树 / Eertree — O(nA)](#13-回文树--eertree--ona)
  - [14. 常用字符串工具](#14-常用字符串工具)

---

## 1. 字符串哈希

```cpp
ull Hash(string s)
{
    ull res = 0;
    for(char c : s)
        res = res * 131 + c;
    return res;
}
```

---

## 2. KMP — 前缀函数 / 模式匹配

```cpp
// pi[i] = s[0..i] 最长真前缀 = 后缀的长度
vector<int> prefix_function(const string &s) {
    int n = sz(s);
    vector<int> pi(n);
    for (int i = 1; i < n; i++) {
        int j = pi[i-1];
        while (j > 0 && s[i] != s[j]) j = pi[j-1];
        if (s[i] == s[j]) j++;
        pi[i] = j;
    }
    return pi;
}
```

**前缀函数应用：**

- 全部 border：`pi[n-1], pi[pi[n-1]-1], ...`
- 最小循环节：`n - pi[n-1]`，整除 n 时即为周期

**每个前缀的最短 border：**

```cpp
vector<int> min_borders(const string &s) {
    auto pi = prefix_function(s);
    int n = sz(s);
    vector<int> mb(n);
    for (int i = 0; i < n; i++) {
        if (pi[i] == 0) mb[i] = 0;
        else mb[i] = mb[pi[i] - 1] ? mb[pi[i] - 1] : pi[i];
    }
    return mb;
}
```

**弱周期引理:**
```
#include <iostream>
#include <string>

using namespace std;

const int N = 1000005; // 根据实际题目修改大小

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

---

## 3. Z 函数

```cpp
// z[i] = s 与 s[i..] 的最长公共前缀长度
vector<int> z_function(const string &s) {
    int n = sz(s);
    vector<int> z(n);
    for (int i = 1, l = 0, r = 0; i < n; i++) {
        if (i < r) z[i] = min(r - i, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
        if (i + z[i] > r) l = i, r = i + z[i];
    }
    return z;
}
```

---

## 4. 字典树

```cpp
struct Trie {
    static constexpr int A = 26;
    struct Node {
        array<int, A> ch{};
        int end_cnt = 0;   // 以此节点结尾的串数
        int pass_cnt = 0;  // 经过此节点的串数
    };
    vector<Node> t = {Node{}};

    static int idx(char c) { return c - 'a'; }

    void insert(const string &s) {
        int v = 0;
        for (char c : s) {
            int x = idx(c);
            if (!t[v].ch[x]) { t[v].ch[x] = sz(t); t.emplace_back(); }
            v = t[v].ch[x];
            t[v].pass_cnt++;
        }
        t[v].end_cnt++;
    }

    int find(const string &s) const {
        int v = 0;
        for (char c : s) {
            int x = idx(c);
            if (!t[v].ch[x]) return -1;
            v = t[v].ch[x];
        }
        return v;
    }
};
```

---

## 5. AC 自动机

```cpp
struct AhoCorasick {
    static constexpr int A = 26;
    struct Node {
        array<int, A> ch{}, go{};
        int fail = 0;       // 失配指针
        int dict_link = 0;  // fail 链上最近的有 ids 的祖先
        vector<int> ids;    // 此处结尾的模式串编号
    };
    vector<Node> t = {Node{}};

    static int idx(char c) { return c - 'a'; }

    int insert(const string &s, int id) {
        int v = 0;
        for (char c : s) {
            int x = idx(c);
            if (!t[v].ch[x]) { t[v].ch[x] = sz(t); t.emplace_back(); }
            v = t[v].ch[x];
        }
        t[v].ids.push_back(id);
        return v;
    }

    void build() {
        queue<int> q;
        for (int x = 0; x < A; x++)
            if (t[0].ch[x]) { t[0].go[x] = t[0].ch[x]; q.push(t[0].ch[x]); }
        while (!q.empty()) {
            int v = q.front(); q.pop();
            int f = t[v].fail;
            t[v].dict_link = !t[f].ids.empty() ? f : t[f].dict_link;
            for (int x = 0; x < A; x++) {
                if (t[v].ch[x]) {
                    int u = t[v].ch[x];
                    t[u].fail = t[f].go[x];
                    t[v].go[x] = u;
                    q.push(u);
                } else {
                    t[v].go[x] = t[f].go[x];
                }
            }
        }
    }

    // 在 text 上匹配，命中 (pattern_id, end_pos) 调用 on_match
    template<class F>
    void match(const string &text, F&& on_match) const {
        int v = 0;
        for (int i = 0; i < sz(text); i++) {
            v = t[v].go[idx(text[i])];
            for (int u = v; u; u = t[u].dict_link)
                for (int id : t[u].ids) on_match(id, i);
        }
    }
};
```

---

## 6. Manacher 算法 — 回文

把 s 变换为 `^#a#b#c#$`，所有奇/偶回文都变成变换串中的奇回文。
**`d[i]` 直接为以 t[i] 为中心的原回文长度**。

```cpp
vector<int> manacher(const string &s) {
    string t = "^";
    for (char c : s) { t += '#'; t += c; }
    t += "#$";
    int n = sz(t);
    vector<int> d(n);
    for (int i = 1, l = 0, r = 0; i + 1 < n; i++) {
        if (i < r) d[i] = min(r - i, d[l + r - i]);
        while (t[i + d[i] + 1] == t[i - d[i] - 1]) d[i]++;
        if (i + d[i] > r) { l = i - d[i]; r = i + d[i]; }
    }
    return d;
}
```

**位置映射（原串长 n，t 长 2n+3）：**

- t[2k+2] = s[k]，奇回文中心
- t[2k+1] 为字符间空隙，偶回文中心
- 原回文左端点 = `(i - d[i]) / 2`，长度 = `d[i]`

---

## 7. 后缀数组 - O(n log n)

```
const int N = 1000005; 

int sa[N];
int rk[N];
int oldrk[N << 1]; 
int tmp[N];
int cnt[N];
int height[N];

void build_sa(const string& s) {
    int n = s.length();
    int m = n;

    fill(cnt, cnt + m + 5, 0);
    fill(oldrk, oldrk + 2 * n + 5, 0);
    fill(height, height + n + 5, 0);
 

    for (int i = 1; i <= n; i++) cnt[rk[i] = s[i - 1]]++;
    for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
    for (int i = n; i >= 1; i--) sa[cnt[rk[i]]--] = i;

    for (int w = 1, p = 0; w < n; w <<= 1, m = p) {
        p = 0;
        for (int i = n - w + 1; i <= n; i++) tmp[++p] = i; 
        for (int i = 1; i <= n; i++) {
            if (sa[i] > w) tmp[++p] = sa[i] - w;
        }

        for (int i = 1; i <= m; i++) cnt[i] = 0;
        for (int i = 1; i <= n; i++) cnt[rk[i]]++;
        for (int i = 1; i <= m; i++) cnt[i] += cnt[i - 1];
        for (int i = n; i >= 1; i--) sa[cnt[rk[tmp[i]]]--] = tmp[i];

        for (int i = 1; i <= n; i++) oldrk[i] = rk[i];
        p = 0;
        for (int i = 1; i <= n; i++) {
            if (oldrk[sa[i]] == oldrk[sa[i - 1]] && oldrk[sa[i] + w] == oldrk[sa[i - 1] + w])
                rk[sa[i]] = p;
            else
                rk[sa[i]] = ++p;
        }
        if (p == n) break; 
    }
}

void build_height(const string& s) {
    int n = s.length();
    int k = 0;
    for (int i = 1; i <= n; i++) {
        if (rk[i] == 1) continue;
        if (k) k--;
        int j = sa[rk[i] - 1];
        while (i + k <= n && j + k <= n && s[i + k - 1] == s[j + k - 1]) k++;
        height[rk[i]] = k;
    }
}
```

**应用：**

- 本质不同子串数：`n*(n+1)/2 - sum(lcp)`
- 双串 LCS：拼接 + 分隔符，扫相邻不同源后缀的最大 LCP
- 最长重复子串：LCP 数组最大值

---

## 压缩后缀树

### 1. 基本概念

将字符串 `s` 的所有后缀插入 Trie，得到后缀 Trie，节点数最坏为 $O(n^2)$。

**压缩后缀树**：把没有分叉的连续路径压缩成一条边；边用原字符串中的区间表示，不实际复制字符。节点数、边数均为 $O(n)$。

构建时在字符串末尾加入唯一终止符 `$`（不出现在原串中），保证每个后缀有独立叶子。

- 内部节点：至少两个后缀的公共前缀，非根内部节点至少有两个儿子。
- 叶子：一个完整后缀。
- **隐式节点**：压缩边中间的某个位置，也代表一个不同子串。

记：

- `dep[u]`：从根到 `u` 的字符串长度。
- `len(u,v)`：边 `u -> v` 的长度。

$$
len(u,v)=dep[v]-dep[u].
$$

> 注意：一个不同子串可能落在某条边的中间，不能只统计显式节点。

### 2. SA + Height 建树

采用 **SA + Height（LCP）+ 栈**，已知 SA、Height 后可在线性时间建树。

- `sa[i]`：字典序第 `i` 小的后缀起点（1-index）。
- `height[i]`：`sa[i-1]` 与 `sa[i]` 的 LCP 长度，`height[1]=0`。

按 SA 顺序插入叶子，维护最右侧路径上的内部节点栈。处理第 `i` 个后缀，令 `h=height[i]`：

1. 弹出深度大于 `h` 的内部节点。
2. 若栈顶深度小于 `h`，就在栈顶最后一个儿子的边上插入深度为 `h` 的新内部节点，将原儿子移到新节点下面。
3. 新建深度为当前后缀长度的叶子，挂到栈顶内部节点下。

**示例**：`s = "abab$"`

| i | sa[i] | 后缀 | height[i] |
|---|---:|---|---:|
| 1 | 5 | `$` | 0 |
| 2 | 3 | `ab$` | 0 |
| 3 | 1 | `abab$` | 2 |
| 4 | 4 | `b$` | 0 |
| 5 | 2 | `bab$` | 1 |

其中 `height[3]=2` 意味着两个相邻后缀有公共前缀 `ab`，需要创建深度为 2 的分叉节点。

#### C++ 模板

以下 `n` **包含 `$`**，`sa[1..n]`、`height[1..n]` 已经求好。

```cpp
struct SuffixTree {
    struct Node {
        int dep = 0, pos = -1, cnt = 0;
        vector<int> son;
    };

    vector<Node> tr;
    int rt = 0;

    int newnode(int dep, int pos = -1) {
        int id = tr.size();
        tr.push_back({});
        tr[id].dep = dep;
        tr[id].pos = pos;
        return id;
    }

    void build(int n, int sa[], int height[]) {
        tr.clear();
        rt = newnode(0);
        vector<int> st = {rt};

        for (int i = 1; i <= n; i++) {
            int h = height[i];
            while (tr[st.back()].dep > h) st.pop_back();

            if (tr[st.back()].dep < h) {
                int p = st.back();
                int v = newnode(h);
                int u = tr[p].son.back();
                tr[p].son.back() = v;
                tr[v].son.push_back(u);
                st.push_back(v);
            }

            int u = newnode(n - sa[i] + 1, sa[i]);
            tr[st.back()].son.push_back(u);
        }
    }

    void dfs(int u) {
        if (tr[u].son.empty()) {
            tr[u].cnt = 1;
            return;
        }
        tr[u].cnt = 0;
        for (int v : tr[u].son) {
            dfs(v);
            tr[u].cnt += tr[v].cnt;
            if (tr[u].pos == -1) tr[u].pos = tr[v].pos;
        }
    }
};
```

- `pos[u]`：`u` 子树某个后缀的起点，可用于定位边字符串。
- `cnt[u]`：子树叶子数，也就是该节点字符串的出现次数（以非终止符结尾的子串为准）。
- 构建树：$O(n)$ 时间、$O(n)$ 空间（不含 SA 的预处理）。

> 注：上述 `dfs` 使用递归，极端长链可能爆栈；大数据可改迭代后序遍历。

### 3. 子串信息怎么维护？

#### 3.1 出现次数

对于边 `u -> v`，在这条边中间任意位置结束的子串，都有相同的出现次数：

$$
\boxed{cnt[v]}.
$$

原因：它们对应的后缀集合完全相同，恰好是 `v` 子树的所有叶子。

#### 3.2 去掉终止符

若只统计原串的非空子串，记有效深度：

$$
D[u]=\begin{cases}
dep[u]-1,&u\text{ 是叶子},\\
dep[u],&u\text{ 是内部节点}.
\end{cases}
$$

对于边 `u -> v`，有效长度为：

$$
\boxed{k=D[v]-D[u]}.
$$

此处原串后缀树的非根内部节点深度均不含终止符；`$` 叶子的有效深度是 0，贡献自然为 0。

#### 3.3 不同子串数量

一条边代表 `k` 个不同子串（对应这条边的 `k` 个结束位置），所以：

$$
\boxed{\text{distinct}=\sum_{u\to v}(D[v]-D[u])}.
$$

#### 3.4 出现至少 / 恰好 t 次的不同子串数

$$
\boxed{\text{atLeast}(t)=\sum_{u\to v,\ cnt[v]\ge t}(D[v]-D[u])}.
$$

$$
\boxed{\text{exact}(t)=\sum_{u\to v,\ cnt[v]=t}(D[v]-D[u])}.
$$

#### 3.5 所有不同子串的长度之和

边 `u -> v` 对应的子串长度为 `D[u]+1, ..., D[v]`，贡献为等差数列：

$$
\boxed{\text{sumLen}=\sum_{u\to v}\frac{(D[u]+1+D[v])(D[v]-D[u])}{2}}.
$$

#### 3.6 定位边对应的字符串

若 `pos[v]` 是 `v` 子树某个后缀的起点，则边 `u -> v` 对应原串区间（1-index）：

$$
\boxed{[pos[v]+dep[u],\ pos[v]+dep[v]-1]}.
$$

因此只需存 `pos` 和 `dep`，不需要复制字符串。叶子末端可能包含 `$`，如要输出原串子串，请用有效深度裁掉末尾终止符。

### 4. 后缀 Trie 的祖先互斥选点 DP 如何优化？

问题：选一些**非空且不同的子串**，要求任意两个选中的子串不存在真前缀关系。它们在后缀 Trie 上对应一个反链（任意两个节点互不为祖先）。

记：`dp[u]` 为 `u` 子树中的合法选择方案数，包含空集；先假设 `u` 对应一个可以选的原串子串。

原始 Trie 的转移：

$$
\boxed{dp[u]=1+\prod_{v\in son(u)}dp[v]}.
$$

- 选 `u`：只有 1 种，其子树后代都不能选。
- 不选 `u`：各个儿子子树独立，方案数相乘。

压缩后，边 `u -> v` 的有效长度为 `k=D[v]-D[u]`。

- **若 `v` 是内部节点**：从 `u` 到 `v` 的边上有 `k-1` 个被压缩的中间节点；选其中任意一个只有一种对应方案，不选这些中间节点则有 `dp[v]` 种，因此贡献为 `dp[v]+k-1`。
- **若 `v` 是叶子**：边上共有 `k` 个有效位置，可以任选其中一个，或者都不选，贡献为 `k+1`。当 `k=0` 时贡献自然为 1。

统一记每个儿子分支的贡献：

$$
\boxed{g(u,v)=\begin{cases}
k+1,&v\text{ 为叶子},\\
dp[v]+k-1,&v\text{ 为内部节点}.
\end{cases}}
$$

非根内部节点代表非空子串，其转移是：

$$
\boxed{dp[u]=1+\prod_{v\in son(u)}g(u,v)}.
$$

根代表空串，**不能选**，因此最终答案是：

$$
\boxed{ans=\prod_{v\in son(root)}g(root,v)}.
$$

这里 `ans` 包含什么都不选的空集方案；若不允许空集，最后减 1。需要取模时，按题目模数计算。


## 8. 后缀自动机 (SAM) — O(nA)

```cpp
struct SAM {
    static constexpr int A = 26;
    struct State {
        int len = 0, link = -1;
        array<int, A> trans{};
        ll cnt = 0; // endpos 大小（出现次数）
    };
    vector<State> st = {{}};
    int last = 0;

    static int idx(char c) { return c - 'a'; }

    void extend(char c) {
        int x = idx(c);
        int cur = sz(st);
        st.emplace_back();
        st[cur].len = st[last].len + 1;
        st[cur].cnt = 1;

        int p = last;
        while (p != -1 && !st[p].trans[x]) {
            st[p].trans[x] = cur;
            p = st[p].link;
        }
        if (p == -1) {
            st[cur].link = 0;
        } else {
            int q = st[p].trans[x];
            if (st[q].len == st[p].len + 1) {
                st[cur].link = q;
            } else {
                int cl = sz(st);
                st.push_back(st[q]);
                st[cl].len = st[p].len + 1;
                st[cl].cnt = 0;
                while (p != -1 && st[p].trans[x] == q) {
                    st[p].trans[x] = cl;
                    p = st[p].link;
                }
                st[q].link = st[cur].link = cl;
            }
        }
        last = cur;
    }

    // 拓扑序累计 endpos 大小
    void compute_cnt() {
        int n = sz(st);
        vector<int> ord(n);
        iota(all(ord), 0);
        sort(all(ord), [&](int a, int b) { return st[a].len > st[b].len; });
        for (int v : ord)
            if (st[v].link != -1) st[st[v].link].cnt += st[v].cnt;
    }

    ll distinct_substrings() const {
        ll ans = 0;
        for (int i = 1; i < sz(st); i++) ans += st[i].len - st[st[i].link].len;
        return ans;
    }
};
```

---

## 9. 最长公共子序列 (LCS) — DP

```cpp
int lcs(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<int> pre(m + 1), cur(m + 1);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++)
            cur[j] = max({pre[j], cur[j-1], pre[j-1] + (a[i-1] == b[j-1])});
        pre = cur;
    }
    return pre[m];
}
```

---

## 10. 最小循环表示 (Booth 算法)

```cpp
// 最小循环表示的起始下标
int min_rotation(const string &s) {
    int n = sz(s);
    string t = s + s;
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        char a = t[i+k], b = t[j+k];
        if (a == b) { k++; continue; }
        if (a < b) j += k + 1;
        else       i += k + 1;
        if (i == j) j++;
        k = 0;
    }
    return min(i, j);
}
```

---

## 11. Lyndon 分解 (Duval 算法)

```cpp
// 把 s 划分为若干 Lyndon 串
vector<string> lyndon_factorize(const string &s) {
    vector<string> res;
    int n = sz(s), i = 0;
    while (i < n) {
        int j = i, k = i + 1;
        while (k < n && s[j] <= s[k]) {
            j = (s[j] < s[k]) ? i : j + 1;
            k++;
        }
        while (i <= j) {
            res.push_back(s.substr(i, k - j));
            i += k - j;
        }
    }
    return res;
}
```

---

## 12. 编辑距离 (Levenshtein)

```cpp
int edit_distance(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<int> pre(m + 1), cur(m + 1);
    iota(all(pre), 0);
    for (int i = 1; i <= n; i++) {
        cur[0] = i;
        for (int j = 1; j <= m; j++)
            cur[j] = min({pre[j] + 1, cur[j-1] + 1,
                          pre[j-1] + (a[i-1] != b[j-1])});
        pre = cur;
    }
    return pre[m];
}
```

---

## 13. 回文树 / Eertree — O(nA)

```cpp
struct PalindromeTree {
    static constexpr int A = 26;
    struct Node {
        int len, link, occ = 0;
        array<int, A> ch{};
    };
    string s;
    vector<Node> t;
    int last;

    static int idx(char c) { return c - 'a'; }

    PalindromeTree() {
        t.push_back({-1, 0, 0, {}}); // 奇根（虚拟，len=-1）
        t.push_back({0, 0, 0, {}});  // 偶根
        last = 1;
    }

    // 沿 fail 链找最长 v，使 s[pos - len(v) - 1] == s[pos]
    int find_link(int v, int pos) {
        while (pos - t[v].len - 1 < 0 || s[pos - t[v].len - 1] != s[pos])
            v = t[v].link;
        return v;
    }

    void extend(int pos) {
        int x = idx(s[pos]);
        int cur = find_link(last, pos);
        if (t[cur].ch[x]) {
            last = t[cur].ch[x];
            t[last].occ++;
            return;
        }
        int now = sz(t);
        t.push_back({t[cur].len + 2, 0, 1, {}});
        t[cur].ch[x] = now;
        t[now].link = (t[now].len == 1) ? 1 : t[find_link(t[cur].link, pos)].ch[x];
        last = now;
    }

    void build(const string &str) {
        s = str;
        for (int i = 0; i < sz(s); i++) extend(i);
        // 累计 occ：拓扑序按创建逆序
        for (int i = sz(t) - 1; i >= 2; i--) t[t[i].link].occ += t[i].occ;
    }

    int distinct_palindromes() const { return sz(t) - 2; }
};
```

---

## 14. 常用字符串工具

```cpp
// t 是否为 s 的子序列
bool is_subsequence(const string &t, const string &s) {
    int j = 0;
    for (int i = 0; i < sz(s) && j < sz(t); i++)
        if (s[i] == t[j]) j++;
    return j == sz(t);
}

// 二分 + 哈希求两后缀的 LCP，模板支持 SingleHash / DoubleHash
template<class H>
int suffix_lcp(const H &h, int i, int j) {
    int lo = 0, hi = min(h.n - i, h.n - j);
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (h.get(i, i + mid) == h.get(j, j + mid)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

// 字典序比较 s[l1..r1) 与 s[l2..r2)
template<class H>
bool substring_less(const H &h, const string &s,
                    int l1, int r1, int l2, int r2) {
    int len1 = r1 - l1, len2 = r2 - l2, m = min(len1, len2);
    int lo = 0, hi = m;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (h.get(l1, l1 + mid) == h.get(l2, l2 + mid)) lo = mid;
        else hi = mid - 1;
    }
    if (lo == m) return len1 < len2;
    return s[l1 + lo] < s[l2 + lo];
}
```
