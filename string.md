# String Algorithms Template for Competitive Programming

## 0. Fast I/O & Utility

```cpp
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ull = unsigned long long;
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
```

---

## 1. String Hashing (Rolling Hash)

```cpp
template<int MOD, int BASE>
struct Hash {
    int n;
    vector<ll> h, p;
    Hash(const string &s) : n(sz(s)), h(n + 1), p(n + 1) {
        p[0] = 1;
        for (int i = 0; i < n; i++) {
            p[i+1] = p[i] * BASE % MOD;
            h[i+1] = (h[i] * BASE + s[i]) % MOD;
        }
    }
    // 子串 s[l..r) 的哈希值
    ll get(int l, int r) const {
        return (h[r] - h[l] * p[r-l] % MOD + MOD) % MOD;
    }
};

using SingleHash = Hash<1'000'000'007, 911382323>;

// 双哈希：两套 (MOD, BASE) 组合，几乎不可能冲突
struct DoubleHash {
    int n;
    Hash<1'000'000'007, 911382323> h1;
    Hash<1'000'000'009, 972663749> h2;
    DoubleHash(const string &s) : n(sz(s)), h1(s), h2(s) {}
    pair<ll,ll> get(int l, int r) const { return {h1.get(l, r), h2.get(l, r)}; }
};
```

---

## 2. KMP — Prefix Function / Pattern Matching

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

// 在 text 中查找 pattern 的所有出现起点
vector<int> find_all(const string &pattern, const string &text) {
    string s = pattern + '\1' + text;
    auto pi = prefix_function(s);
    int m = sz(pattern);
    vector<int> occ;
    for (int i = m + 1; i < sz(s); i++)
        if (pi[i] == m) occ.push_back(i - 2 * m);
    return occ;
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

---

## 3. Z-Function

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

## 4. Trie

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

## 5. Aho-Corasick Automaton

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

## 6. Manacher's Algorithm — Palindromes

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

## 7. Suffix Array — O(n log² n)

```cpp
struct SuffixArray {
    int n;
    string s;
    vector<int> sa, rk, lcp; // lcp[i] = LCP(sa[i], sa[i+1])

    SuffixArray(const string &_s) : n(sz(_s)), s(_s) {
        string t = s + char(0);
        int N = n + 1;
        sa.resize(N); iota(all(sa), 0);
        sort(all(sa), [&](int i, int j) { return t[i] < t[j]; });
        rk.assign(N, 0);
        for (int i = 1; i < N; i++)
            rk[sa[i]] = rk[sa[i-1]] + (t[sa[i]] != t[sa[i-1]]);

        for (int k = 1; k < N && rk[sa[N-1]] < N - 1; k <<= 1) {
            auto cmp = [&](int i, int j) {
                if (rk[i] != rk[j]) return rk[i] < rk[j];
                int ri = i + k < N ? rk[i+k] : -1;
                int rj = j + k < N ? rk[j+k] : -1;
                return ri < rj;
            };
            sort(all(sa), cmp);
            vector<int> nrk(N);
            for (int i = 1; i < N; i++)
                nrk[sa[i]] = nrk[sa[i-1]] + cmp(sa[i-1], sa[i]);
            rk = move(nrk);
        }
        sa.erase(sa.begin()); // 去掉哨兵

        // Kasai 求 LCP
        if (n > 1) lcp.assign(n - 1, 0);
        rk.assign(n, 0);
        for (int i = 0; i < n; i++) rk[sa[i]] = i;
        for (int i = 0, k = 0; i < n; i++) {
            if (rk[i] == 0) { k = 0; continue; }
            int j = sa[rk[i] - 1];
            while (i + k < n && j + k < n && s[i+k] == s[j+k]) k++;
            lcp[rk[i] - 1] = k;
            if (k) k--;
        }
    }

    // 后缀 i 与 j 的 LCP（暴力 min，配合 RMQ 可优化为 O(1)）
    int get_lcp(int i, int j) const {
        if (i == j) return n - i;
        int a = rk[i], b = rk[j];
        if (a > b) swap(a, b);
        return *min_element(lcp.begin() + a, lcp.begin() + b);
    }
};
```

**应用：**

- 本质不同子串数：`n*(n+1)/2 - sum(lcp)`
- 双串 LCS：拼接 + 分隔符，扫相邻不同源后缀的最大 LCP
- 最长重复子串：LCP 数组最大值

---

## 8. Suffix Automaton (SAM) — O(n A)

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

## 9. Longest Common Subsequence (LCS) — DP

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

需要还原 LCS 字符串时要保留二维 dp，参见 [string_applications.md](string_applications.md)。

---

## 10. Lexicographically Minimal Rotation (Booth's Algorithm)

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

## 11. Lyndon Factorization (Duval's Algorithm)

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

## 12. Edit Distance (Levenshtein)

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

## 13. Palindrome Tree / Eertree — O(n A)

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

## 14. Common String Utilities

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

---

## Complexity Summary

| Algorithm | Build Time | Query/Match | Memory |
| --------- | ---------- | ----------- | ------ |
| Hashing | O(n) | O(1) | O(n) |
| KMP | O(n) | O(n+m) | O(n) |
| Z-Function | O(n) | O(n+m) | O(n) |
| Trie | O(sum\|s\|) | O(\|s\|) | O(sum\|s\| · A) |
| Aho-Corasick | O(sum\|p\|) | O(\|t\| + #matches) | O(sum\|p\| · A) |
| Manacher | O(n) | — | O(n) |
| Suffix Array | O(n log² n) | O(区间) (RMQ→O(1)) | O(n) |
| SAM | O(n A) | O(\|s\|) | O(n A) |
| Palindrome Tree | O(n A) | — | O(n A) |

---

## When to Use What

| Problem type | Recommended |
| ------------ | ----------- |
| 子串相等 / 字典序比较 | Rolling hash |
| 单模式匹配 | KMP / Z |
| 多模式匹配 | Aho-Corasick |
| 回文查询 | Manacher / Palindrome Tree |
| 双串最长公共子串 | Suffix Array / SAM |
| 本质不同子串数 | SAM / Suffix Array |
| 字典序 K 小子串 | SAM (DAG 上 DP) |
| 子串出现次数 | SAM (cnt) |
| 编辑距离（小串） | DP O(n·m) |
