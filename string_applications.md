# 字符串算法应用详解

本文针对 [string.md](string.md) 中的算法模板，给出每个算法的**典型应用场景**与**完整实现代码**，便于在比赛中快速套用。

---

## 目录

1. [字符串哈希](#1-字符串哈希)
2. [KMP 算法](#2-kmp-算法)
3. [Z-Function](#3-z-function)
4. [Trie 字典树](#4-trie-字典树)
5. [AC 自动机](#5-ac-自动机)
6. [Manacher 算法](#6-manacher-算法)
7. [后缀数组 SA](#7-后缀数组-sa)
8. [后缀自动机 SAM](#8-后缀自动机-sam)
9. [LCS 最长公共子序列](#9-lcs-最长公共子序列)
10. [最小表示法 / Lyndon 分解](#10-最小表示法--lyndon-分解)
11. [编辑距离](#11-编辑距离)
12. [回文树 PAM](#12-回文树-pam)

---

## 1. 字符串哈希

**应用场景：**

- **任意子串相等性判断**：O(1) 比较两段是否相等
- **本质不同子串数量统计**：所有子串哈希塞入 set
- **回文判断**：正反两个哈希比对
- **二分 + 哈希求 LCP**：求两后缀最长公共前缀
- **字典序比较子串**：二分 LCP + 比较下一字符

### 应用 1.1：判断两子串是否相等

```cpp
SingleHash h(s);
bool equal = (h.get(l1, r1 + 1) == h.get(l2, r2 + 1));
```

### 应用 1.2：求两后缀的最长公共前缀

模板已提供 `suffix_lcp(h, i, j)`，二分 + 哈希。

### 应用 1.3：判断子串是否为回文

```cpp
struct PalCheck {
    int n;
    SingleHash fwd, bwd;
    PalCheck(const string &s) : n(sz(s)), fwd(s), bwd(string(s.rbegin(), s.rend())) {}
    bool is_pal(int l, int r) const { // [l, r] 闭区间
        return fwd.get(l, r + 1) == bwd.get(n - 1 - r, n - l);
    }
};
```

### 应用 1.4：统计本质不同子串数量

```cpp
ll count_distinct_brute(const string &s) {
    int n = sz(s);
    DoubleHash dh(s);
    set<pair<ll,ll>> st;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j <= n; j++)
            st.insert(dh.get(i, j));
    return sz(st);
}
// 大数据应改用 SAM 或 SA
```

---

## 2. KMP 算法

**应用场景：**

- **模式串匹配**：在文本 t 中找模式串 p 的所有出现位置
- **求字符串的最小循环节**：`n - pi[n-1]` 整除 n 即为周期
- **求所有 border**：`pi[n-1], pi[pi[n-1]-1], ...`
- **统计前缀出现次数**：每个前缀作为子串出现次数
- **字符串压缩**：判断能否表示为 k 次重复

### 应用 2.1：模式串多次出现位置

```cpp
auto pos = find_all(p, t);
```

### 应用 2.2：求最小循环节

```cpp
int min_period(const string &s) {
    auto pi = prefix_function(s);
    int n = sz(s), p = n - pi[n - 1];
    return (n % p == 0) ? p : n;
}
```

### 应用 2.3：每个前缀作为子串的出现次数

```cpp
// ans[i]: 长度为 i 的前缀在 s 中出现次数
vector<int> count_prefix_occurrences(const string &s) {
    int n = sz(s);
    auto pi = prefix_function(s);
    vector<int> ans(n + 1);
    for (int i = 0; i < n; i++) ans[pi[i]]++;
    for (int i = n - 1; i > 0; i--) ans[pi[i - 1]] += ans[i];
    for (int i = 1; i <= n; i++) ans[i]++; // 前缀本身
    return ans;
}
```

### 应用 2.4：两串的公共 border 长度

```cpp
// 求 t 与 s 的 (后缀, 前缀) 最长匹配长度，限长 ≤ min(|s|,|t|)
int common_border(const string &s, const string &t) {
    string c = t + '\1' + s;
    auto pi = prefix_function(c);
    int v = pi.back();
    while (v > min(sz(s), sz(t))) v = pi[v - 1];
    return v;
}
```

---

## 3. Z-Function

**应用场景：**

- **模式匹配**：与 KMP 等价，实现更简
- **最长公共扩展 LCE**：`z[i]` 直接给出 s 与 s[i..] 的 LCP
- **字符串压缩**

### 应用 3.1：基于 Z 的模式匹配

```cpp
vector<int> find_all_z(const string &p, const string &t) {
    string s = p + '\1' + t;
    auto z = z_function(s);
    vector<int> res;
    int m = sz(p);
    for (int i = m + 1; i < sz(s); i++)
        if (z[i] >= m) res.push_back(i - m - 1);
    return res;
}
```

### 应用 3.2：最大重复次数（s = t^k 的最大 k）

```cpp
int max_repeat(const string &s) {
    int n = sz(s);
    auto z = z_function(s);
    for (int len = 1; len <= n; len++) {
        if (n % len) continue;
        bool ok = true;
        for (int i = len; i < n; i += len)
            if (z[i] < len) { ok = false; break; }
        if (ok) return n / len;
    }
    return 1;
}
```

---

## 4. Trie 字典树

**应用场景：**

- **多模式串前缀查询**：给定字典，查询某串是否为字典中某词前缀
- **01-Trie 异或最值**：`a^b` 最大/最小值
- **统计前缀出现次数**：`pass_cnt` 字段

### 应用 4.1：01-Trie 求异或最大值

```cpp
struct XorTrie {
    static constexpr int B = 30; // 位数
    struct Node {
        array<int, 2> ch{};
        int cnt = 0;
    };
    vector<Node> t = {Node{}};

    void insert(int x, int delta = 1) {
        int v = 0;
        for (int b = B; b >= 0; b--) {
            int d = (x >> b) & 1;
            if (!t[v].ch[d]) { t[v].ch[d] = sz(t); t.emplace_back(); }
            v = t[v].ch[d];
            t[v].cnt += delta;
        }
    }
    void erase(int x) { insert(x, -1); }

    // 与 x 异或最大的值
    int query_max(int x) const {
        int v = 0, res = 0;
        for (int b = B; b >= 0; b--) {
            int d = (x >> b) & 1;
            int u = t[v].ch[d ^ 1];
            if (u && t[u].cnt) { res |= 1 << b; v = u; }
            else v = t[v].ch[d];
        }
        return res;
    }
};
```

### 应用 4.2：判断字典中是否存在某串的前缀

```cpp
bool has_prefix_of(const Trie &tr, const string &q) {
    int v = 0;
    for (char c : q) {
        if (tr.t[v].end_cnt > 0) return true;
        int x = Trie::idx(c);
        if (!tr.t[v].ch[x]) return false;
        v = tr.t[v].ch[x];
    }
    return tr.t[v].end_cnt > 0;
}
```

---

## 5. AC 自动机

**应用场景：**

- **多模式串同时匹配**：在文本中找出多个模式串的全部出现位置
- **病毒检测、敏感词过滤**
- **AC 自动机上 DP**：例如统计不含任何模式串的字符串数量

### 应用 5.1：统计每个模式串出现次数

```cpp
AhoCorasick ac;
for (int i = 0; i < sz(patterns); i++) ac.insert(patterns[i], i);
ac.build();

vector<int> cnt(sz(patterns), 0);
ac.match(text, [&](int id, int /*pos*/) { cnt[id]++; });
```

### 应用 5.2：长度 L 且不含任何禁词的字符串数量

```cpp
const int MOD = 1e9 + 7;

ll count_safe_strings(const vector<string> &patterns, int L) {
    AhoCorasick ac;
    for (int i = 0; i < sz(patterns); i++) ac.insert(patterns[i], i);
    ac.build();

    int N = sz(ac.t);
    // bad[v]: 节点 v 经 fail 链可达终点（即包含模式串）
    vector<bool> bad(N, false);
    for (int v = 1; v < N; v++)
        bad[v] = !ac.t[v].ids.empty() || bad[ac.t[v].fail];

    vector<ll> dp(N, 0); dp[0] = 1;
    for (int step = 0; step < L; step++) {
        vector<ll> nd(N, 0);
        for (int v = 0; v < N; v++) {
            if (!dp[v] || bad[v]) continue;
            for (int x = 0; x < AhoCorasick::A; x++) {
                int u = ac.t[v].go[x];
                if (!bad[u]) (nd[u] += dp[v]) %= MOD;
            }
        }
        dp = move(nd);
    }
    ll ans = 0;
    for (ll x : dp) (ans += x) %= MOD;
    return ans;
}
```

---

## 6. Manacher 算法

**应用场景：**

- **求最长回文子串**
- **统计回文子串总数**
- **以每个位置结尾的最长回文长度**
- **任意子串是否为回文（配合预处理）**

> 模板返回的 `d[i]` 即为以变换串中 t[i] 为中心的**原回文长度**。

### 应用 6.1：最长回文子串

```cpp
pair<int,int> longest_palindrome(const string &s) {
    auto d = manacher(s);
    int best = 0, center = 1;
    for (int i = 1; i + 1 < sz(d); i++)
        if (d[i] > best) { best = d[i]; center = i; }
    int l = (center - best) / 2;
    return {l, l + best - 1};
}
```

### 应用 6.2：统计回文子串总数

```cpp
ll count_palindromes(const string &s) {
    auto d = manacher(s);
    ll ans = 0;
    for (int i = 1; i + 1 < sz(d); i++) ans += (d[i] + 1) / 2;
    return ans;
}
```

> 推导：长度 d 的最长回文中嵌套 (d+1)/2 个同心回文。
> 奇中心 d 必为奇数，(d+1)/2 = 半径数；偶中心 d 必为偶数，(d+1)/2 = d/2。

### 应用 6.3：以每个位置结尾的最长回文长度

```cpp
vector<int> longest_pal_ending_at(const string &s) {
    int n = sz(s);
    auto d = manacher(s);
    vector<int> res(n, 1);
    // 对变换串中每个回文中心 [l, r]，原右端点 = (i + d[i] - 2) / 2
    for (int i = 1; i + 1 < sz(d); i++) {
        int len = d[i];
        int r = (i + len) / 2 - 1;
        if (r >= 0 && r < n) res[r] = max(res[r], len);
    }
    // 长度递推：以 i 结尾的最长回文不会比以 i-1 结尾的少超过 2
    for (int i = 1; i < n; i++) res[i] = max(res[i], res[i-1] - 2);
    return res;
}
```

---

## 7. 后缀数组 SA

**应用场景：**

- **本质不同子串数**：`n*(n+1)/2 - sum(lcp)`
- **两字符串最长公共子串**
- **最长重复子串**：LCP 数组最大值
- **字典序第 k 小子串**

### 应用 7.1：本质不同子串数量

```cpp
ll distinct_substrings_sa(const string &s) {
    SuffixArray sa(s);
    ll total = (ll)sa.n * (sa.n + 1) / 2;
    for (int x : sa.lcp) total -= x;
    return total;
}
```

### 应用 7.2：两串最长公共子串

```cpp
int longest_common_substring(const string &a, const string &b) {
    string s = a + char(1) + b;
    SuffixArray sa(s);
    int na = sz(a), ans = 0;
    for (int i = 1; i < sa.n; i++) {
        int p = sa.sa[i], q = sa.sa[i - 1];
        // 一个在 a 中、一个在 b 中
        if ((p < na) != (q < na)) ans = max(ans, sa.lcp[i - 1]);
    }
    return ans;
}
```

### 应用 7.3：最长重复子串

```cpp
int longest_repeated_substring(const string &s) {
    SuffixArray sa(s);
    return sa.lcp.empty() ? 0 : *max_element(all(sa.lcp));
}
```

### 应用 7.4：字典序第 k 小本质不同子串

```cpp
string kth_distinct_substring(const string &s, ll k) {
    SuffixArray sa(s);
    ll cur = 0;
    for (int i = 0; i < sa.n; i++) {
        int prev_lcp = (i == 0) ? 0 : sa.lcp[i - 1];
        ll add = (sa.n - sa.sa[i]) - prev_lcp;
        if (cur + add >= k) {
            int len = prev_lcp + (int)(k - cur);
            return s.substr(sa.sa[i], len);
        }
        cur += add;
    }
    return ""; // k 越界
}
```

---

## 8. 后缀自动机 SAM

**应用场景：**

- **本质不同子串数**：`sum(len[v] - len[link[v]])`
- **子串出现次数**：`cnt` 字段
- **两串最长公共子串（在线）**
- **字典序第 k 大子串**

### 应用 8.1：本质不同子串数

```cpp
SAM sam;
for (char c : s) sam.extend(c);
ll ans = sam.distinct_substrings();
```

### 应用 8.2：每个子串的出现次数

```cpp
ll occurrences(const SAM &sam, const string &t) {
    int v = 0;
    for (char c : t) {
        int x = SAM::idx(c);
        if (!sam.st[v].trans[x]) return 0;
        v = sam.st[v].trans[x];
    }
    return sam.st[v].cnt;
}
// 使用前需调用 sam.compute_cnt()
```

### 应用 8.3：两串最长公共子串

```cpp
int lcs_sam(const string &a, const string &b) {
    SAM sam;
    for (char c : a) sam.extend(c);
    int v = 0, l = 0, ans = 0;
    for (char c : b) {
        int x = SAM::idx(c);
        while (v && !sam.st[v].trans[x]) {
            v = sam.st[v].link;
            l = sam.st[v].len;
        }
        if (sam.st[v].trans[x]) { v = sam.st[v].trans[x]; l++; }
        ans = max(ans, l);
    }
    return ans;
}
```

### 应用 8.4：字典序第 k 小子串

```cpp
// distinct=true: 本质不同；false: 计算重复（需先 compute_cnt）
string kth_substring_sam(const SAM &sam, ll k, bool distinct = true) {
    int n = sz(sam.st);
    vector<int> ord(n);
    iota(all(ord), 0);
    sort(all(ord), [&](int a, int b) { return sam.st[a].len > sam.st[b].len; });

    vector<ll> dp(n);
    for (int v : ord) {
        dp[v] = (v == 0) ? 0 : (distinct ? 1 : sam.st[v].cnt);
        for (int x = 0; x < SAM::A; x++) {
            int u = sam.st[v].trans[x];
            if (u) dp[v] += dp[u];
        }
    }
    if (k > dp[0]) return ""; // 越界

    string res;
    int v = 0;
    while (k > 0) {
        for (int x = 0; x < SAM::A; x++) {
            int u = sam.st[v].trans[x];
            if (!u) continue;
            if (dp[u] >= k) { res += char('a' + x); v = u; k -= distinct ? 1 : sam.st[u].cnt; break; }
            k -= dp[u];
        }
    }
    return res;
}
```

---

## 9. LCS 最长公共子序列

**应用场景：**

- **DNA 比对、文件 diff**
- **字符串相似度度量**
- **特殊场景**：当其中一个串元素互不相同时，LCS → LIS（O(n log n)）

### 应用 9.1：还原 LCS 串

```cpp
string lcs_string(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<vector<int>> dp(n + 1, vector<int>(m + 1));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dp[i][j] = (a[i-1] == b[j-1]) ? dp[i-1][j-1] + 1
                                          : max(dp[i-1][j], dp[i][j-1]);
    string res;
    for (int i = n, j = m; i && j; ) {
        if (a[i-1] == b[j-1]) { res += a[i-1]; i--; j--; }
        else if (dp[i-1][j] >= dp[i][j-1]) i--;
        else j--;
    }
    reverse(all(res));
    return res;
}
```

### 应用 9.2：排列 LCS → LIS（O(n log n)）

```cpp
// 当 a, b 都是 1..n 的排列时，LCS(a, b) = LIS(map(b))
int lcs_permutation(const vector<int> &a, const vector<int> &b) {
    int n = sz(a);
    vector<int> pos(n + 1);
    for (int i = 0; i < n; i++) pos[a[i]] = i;

    vector<int> tails;
    for (int x : b) {
        int p = pos[x];
        auto it = lower_bound(all(tails), p);
        if (it == tails.end()) tails.push_back(p);
        else *it = p;
    }
    return sz(tails);
}
```

---

## 10. 最小表示法 / Lyndon 分解

**应用场景：**

- **循环字符串归一化**：判断两个串是否循环同构
- **最小循环表示**

### 应用 10.1：判断两串循环同构

```cpp
bool cyclic_equal(const string &a, const string &b) {
    if (sz(a) != sz(b)) return false;
    int n = sz(a), i = min_rotation(a), j = min_rotation(b);
    for (int k = 0; k < n; k++)
        if (a[(i + k) % n] != b[(j + k) % n]) return false;
    return true;
}
```

### 应用 10.2：求字符串的最小循环表示

```cpp
string min_cyclic(const string &s) {
    int idx = min_rotation(s);
    return s.substr(idx) + s.substr(0, idx);
}
```

### 应用 10.3：Duval 分解的另一种应用 —— 求最小表示起点

```cpp
// 在 s+s 上做 Lyndon 分解，第一个长度 ≥ n 的段起点即为答案
int min_rotation_duval(const string &s) {
    int n = sz(s);
    string t = s + s;
    int ans = 0, i = 0;
    while (i < n) {
        ans = i;
        int j = i + 1, k = i;
        while (j < sz(t) && t[k] <= t[j]) {
            k = (t[k] < t[j]) ? i : k + 1;
            j++;
        }
        while (i <= k) i += j - k;
    }
    return ans;
}
```

---

## 11. 编辑距离

**应用场景：**

- **拼写纠错、DNA 比对**
- **字符串相似度**
- **带权变种**：增/删/替换不同代价

### 应用 11.1：还原编辑操作序列

```cpp
vector<string> edit_ops(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<vector<int>> dp(n + 1, vector<int>(m + 1));
    for (int i = 0; i <= n; i++) dp[i][0] = i;
    for (int j = 0; j <= m; j++) dp[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dp[i][j] = min({dp[i-1][j] + 1, dp[i][j-1] + 1,
                            dp[i-1][j-1] + (a[i-1] != b[j-1])});
    vector<string> ops;
    for (int i = n, j = m; i || j; ) {
        if (i && j && a[i-1] == b[j-1] && dp[i][j] == dp[i-1][j-1]) {
            ops.push_back("KEEP " + string(1, a[i-1])); i--; j--;
        } else if (i && j && dp[i][j] == dp[i-1][j-1] + 1) {
            ops.push_back("REPLACE " + string(1, a[i-1]) + "->" + string(1, b[j-1])); i--; j--;
        } else if (i && dp[i][j] == dp[i-1][j] + 1) {
            ops.push_back("DELETE " + string(1, a[i-1])); i--;
        } else {
            ops.push_back("INSERT " + string(1, b[j-1])); j--;
        }
    }
    reverse(all(ops));
    return ops;
}
```

### 应用 11.2：带权编辑距离

```cpp
int weighted_edit(const string &a, const string &b,
                  int c_ins, int c_del, int c_rep) {
    int n = sz(a), m = sz(b);
    vector<int> pre(m + 1), cur(m + 1);
    for (int j = 0; j <= m; j++) pre[j] = j * c_ins;
    for (int i = 1; i <= n; i++) {
        cur[0] = i * c_del;
        for (int j = 1; j <= m; j++)
            cur[j] = min({pre[j] + c_del, cur[j-1] + c_ins,
                          pre[j-1] + (a[i-1] == b[j-1] ? 0 : c_rep)});
        pre = cur;
    }
    return pre[m];
}
```

---

## 12. 回文树 PAM

**应用场景：**

- **本质不同回文子串数**
- **每种回文出现次数**
- **以每个位置结尾的回文数量**

### 应用 12.1：本质不同回文 + 每种出现次数

```cpp
PalindromeTree pt;
pt.build(s);
cout << pt.distinct_palindromes() << "\n";
ll max_pal_value = 0;
for (int i = 2; i < sz(pt.t); i++)
    max_pal_value = max(max_pal_value, (ll)pt.t[i].occ * pt.t[i].len);
```

### 应用 12.2：以每个位置结尾的回文数量

```cpp
vector<int> palindromes_ending(const string &s) {
    PalindromeTree pt;
    pt.s = s;
    int n = sz(s);
    vector<int> ans(n);
    for (int i = 0; i < n; i++) {
        pt.extend(i);
        // 沿 last 的 fail 链遍历到偶根：经过的节点数即为回文数
        int v = pt.last, cnt = 0;
        while (v > 1) { cnt++; v = pt.t[v].link; }
        ans[i] = cnt;
    }
    return ans;
}
// 注：暴力上行 O(n*depth)。需 O(n log n) 时维护 diff/series 链。
```

---

## 综合应用：竞赛题型对应

| 题型 | 推荐算法 |
| ---- | -------- |
| 多次询问子串相等 | Hash |
| 子串字典序排序 | SA |
| 单串本质不同子串 | SAM 或 SA |
| 多模式串匹配 | AC 自动机 |
| 最长回文子串 | Manacher |
| 回文相关复杂统计 | PAM |
| 循环串归一化 | 最小表示 |
| DNA 比对 | LCS / Edit Distance |
| 异或最值（数组中） | 01-Trie |
| AC 自动机上 DP | AC + DP |

---

## 复杂度与选择参考

- $n \le 10^3$：DP（LCS、Edit Distance）足够
- $n \le 10^5$：Hash / KMP / Z / Manacher / SA / SAM 任选
- $n \le 10^6$：避免 SA（log 因子大），优先 SAM / KMP / Hash
- 多模式串：AC 自动机
- 在线扩展字符串：SAM
- 回文专题：Manacher（统计）+ PAM（高级查询）
