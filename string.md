# String Algorithms Template for Competitive Programming

## 0. Fast I/O & Utility

```cpp
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using ull = unsigned long long;
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()

const int MOD = 1e9 + 7;
const int BASE = 911382323;  // large random prime
const int BASE2 = 972663749;  // second base for double hash
```

---

## 1. String Hashing (Rolling Hash) — Single & Double

Single modulus (faster, collision possible on adversarial tests):

```cpp
struct StringHash {
    int n;
    vector<ll> h, p;
    
    StringHash(const string &s, int base = BASE, int mod = MOD) : n(sz(s)), h(n+1), p(n+1) {
        p[0] = 1;
        for (int i = 0; i < n; i++) {
            p[i+1] = p[i] * base % mod;
            h[i+1] = (h[i] * base + s[i]) % mod;
        }
    }
    
    ll get(int l, int r) { // [l, r)  0-indexed
        return (h[r] - h[l] * p[r-l] % MOD + MOD) % MOD;
    }
};
```

Double hash (virtually collision-proof):

```cpp
struct DoubleHash {
    int n;
    vector<pair<ll,ll>> h, p;
    static const int M1 = 1e9+7, M2 = 1e9+9;
    static const int B1 = 911382323, B2 = 972663749;
    
    DoubleHash(const string &s) : n(sz(s)), h(n+1), p(n+1) {
        p[0] = {1, 1};
        for (int i = 0; i < n; i++) {
            p[i+1] = {p[i].first * B1 % M1, p[i].second * B2 % M2};
            h[i+1] = {
                (h[i].first * B1 + s[i]) % M1,
                (h[i].second * B2 + s[i]) % M2
            };
        }
    }
    
    pair<ll,ll> get(int l, int r) { // [l, r)
        auto sub = [&](ll hr, ll hl, ll pw, int M) {
            return (hr - hl * pw % M + M) % M;
        };
        return {
            sub(h[r].first, h[l].first, p[r-l].first, M1),
            sub(h[r].second, h[l].second, p[r-l].second, M2)
        };
    }
};
```

---

## 2. KMP — Prefix Function / Pattern Matching

```cpp
// pi[i] = longest proper prefix of s[0..i] that is also a suffix
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

// Find all occurrences of pattern p in text t
vector<int> kmp(const string &p, const string &t) {
    string s = p + '#' + t;
    auto pi = prefix_function(s);
    vector<int> occ;
    for (int i = sz(p)+1; i < sz(s); i++)
        if (pi[i] == sz(p))
            occ.push_back(i - 2*sz(p));
    return occ;
}
```

**Applications of prefix function:**
- Border of a string: `pi[n-1]`, next border: `pi[pi[n-1]-1]`, ...
- Minimum period: `n - pi[n-1]` (divides n iff string is periodic)
- Compress string: repeat of `n - pi[n-1]` length prefix

---

## 3. Z-Function

```cpp
// z[i] = LCP of s and s[i..]
vector<int> z_function(const string &s) {
    int n = sz(s);
    vector<int> z(n);
    int l = 0, r = 0;
    for (int i = 1; i < n; i++) {
        if (i <= r) z[i] = min(r - i + 1, z[i - l]);
        while (i + z[i] < n && s[z[i]] == s[i + z[i]]) z[i]++;
        if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
    }
    return z;
}

// Pattern matching with Z: build Z on p + '#' + t, then find z[i] == |p|
```

---

## 4. Trie

```cpp
struct Trie {
    struct Node {
        int nxt[26] = {};
        int cnt = 0;   // how many words end here
        int pass = 0;  // how many words pass through
    };
    vector<Node> t = {Node{}};
    
    void add(const string &s) {
        int v = 0;
        for (char c : s) {
            int x = c - 'a';
            if (!t[v].nxt[x]) {
                t[v].nxt[x] = sz(t);
                t.emplace_back();
            }
            v = t[v].nxt[x];
            t[v].pass++;
        }
        t[v].cnt++;
    }
    
    int find(const string &s) {
        int v = 0;
        for (char c : s) {
            int x = c - 'a';
            if (!t[v].nxt[x]) return -1;
            v = t[v].nxt[x];
        }
        return v; // return node index
    }
};
```

---

## 5. Aho-Corasick Automaton

```cpp
struct AhoCorasick {
    struct Node {
        int nxt[26] = {}, go[26] = {}; // nxt: trie edges, go: transition
        int link = 0, out = 0;         // out: first output link
        vector<int> end_ids;           // pattern ids ending here
    };
    vector<Node> t = {Node{}};
    
    int add(const string &s, int id) {
        int v = 0;
        for (char c : s) {
            int x = c - 'a';
            if (!t[v].nxt[x]) {
                t[v].nxt[x] = sz(t);
                t.emplace_back();
            }
            v = t[v].nxt[x];
        }
        t[v].end_ids.push_back(id);
        return v;
    }
    
    void build() {
        queue<int> q;
        for (int x = 0; x < 26; x++) {
            if (t[0].nxt[x]) {
                t[0].go[x] = t[0].nxt[x];
                q.push(t[0].nxt[x]);
            } else {
                t[0].go[x] = 0;
            }
        }
        while (!q.empty()) {
            int v = q.front(); q.pop();
            // merge output of suffix link
            int sl = t[v].link;
            t[v].out = t[sl].end_ids.empty() ? t[sl].out : sl;
            for (int x = 0; x < 26; x++) {
                if (t[v].nxt[x]) {
                    int u = t[v].nxt[x];
                    t[u].link = t[sl].go[x];
                    t[v].go[x] = u;
                    q.push(u);
                } else {
                    t[v].go[x] = t[sl].go[x];
                }
            }
        }
    }
    
    // Traverse text and call lambda for each match (pattern_id, end_pos)
    template<typename F>
    void match(const string &text, F&& on_match) {
        int v = 0;
        for (int i = 0; i < sz(text); i++) {
            v = t[v].go[text[i] - 'a'];
            for (int u = v; u; u = t[u].out) {
                for (int id : t[u].end_ids)
                    on_match(id, i);
            }
        }
    }
};
```

---

## 6. Manacher's Algorithm — Palindromes

```cpp
// Returns array d where d[i] = palindrome radius at center i
// d1 for odd-length, d2 for even-length
// Number of palindromes starting/ending at i can be derived from d.

// Odd-length palindromes: center i, radius d1[i], length = 2*d1[i]-1
vector<int> manacher_odd(const string &s) {
    int n = sz(s);
    vector<int> d1(n);
    for (int i = 0, l = 0, r = -1; i < n; i++) {
        int k = (i > r) ? 1 : min(d1[l + r - i], r - i + 1);
        while (i - k >= 0 && i + k < n && s[i-k] == s[i+k]) k++;
        d1[i] = k--;
        if (i + k > r) l = i - k, r = i + k;
    }
    return d1;
}

// Even-length palindromes: center between i-1 and i, radius d2[i], length = 2*d2[i]
vector<int> manacher_even(const string &s) {
    int n = sz(s);
    vector<int> d2(n);
    for (int i = 0, l = 0, r = -1; i < n; i++) {
        int k = (i > r) ? 0 : min(d2[l + r - i + 1], r - i + 1);
        while (i - k - 1 >= 0 && i + k < n && s[i-k-1] == s[i+k]) k++;
        d2[i] = k--;
        if (i + k > r) l = i - k - 1, r = i + k;
    }
    return d2;
}
```

---

## 7. Suffix Array — O(n log n)

```cpp
struct SuffixArray {
    int n;
    string s;
    vector<int> sa, rk, lcp; // sa: suffix array, rk: rank, lcp: Kasai
    
    SuffixArray(const string &_s) : s(_s), n(sz(s)) {
        s += '$'; n++; // sentinel, smaller than all chars
        sa.resize(n); iota(all(sa), 0);
        
        // Sort by first character
        sort(all(sa), [&](int i, int j) { return s[i] < s[j]; });
        rk.resize(n); rk[sa[0]] = 0;
        for (int i = 1; i < n; i++)
            rk[sa[i]] = rk[sa[i-1]] + (s[sa[i]] != s[sa[i-1]]);
        
        for (int k = 1; k < n; k <<= 1) {
            auto cmp = [&](int i, int j) {
                if (rk[i] != rk[j]) return rk[i] < rk[j];
                int ri = i + k < n ? rk[i+k] : -1;
                int rj = j + k < n ? rk[j+k] : -1;
                return ri < rj;
            };
            sort(all(sa), cmp);
            vector<int> nrk(n);
            nrk[sa[0]] = 0;
            for (int i = 1; i < n; i++)
                nrk[sa[i]] = nrk[sa[i-1]] + cmp(sa[i-1], sa[i]);
            rk.swap(nrk);
            if (rk[sa[n-1]] == n-1) break;
        }
        
        // Kasai's algorithm for LCP
        lcp.resize(n-1);
        for (int i = 0, k = 0; i < n-1; i++) { // n-1: ignore sentinel
            int j = sa[rk[i] - 1];
            while (s[i+k] == s[j+k]) k++;
            lcp[rk[i]-1] = k;
            if (k) k--;
        }
        // Remove sentinel
        sa.erase(sa.begin());
        s.pop_back(); n--;
    }
    
    // LCP of suffixes starting at i and j (require rk access)
    int get_lcp(int i, int j) {
        if (i == j) return n - i;
        int a = rk[i], b = rk[j];
        if (a > b) swap(a, b);
        return *min_element(lcp.begin() + a, lcp.begin() + b);
    }
};
```

**Applications:**
- Count distinct substrings: `n*(n+1)/2 - sum(lcp)`
- Longest common substring of 2 strings: concatenate with separator `#`, build SA, max LCP where suffixes from different strings
- Longest repeated substring: max value in LCP array

---

## 8. Suffix Automaton (SAM) — O(n)

Compact representation of all substrings. Each state represents a set of endpos-equivalent substrings.

```cpp
struct SAM {
    struct State {
        int len, link;            // max length, suffix link
        map<char,int> next;       // transitions (use array for speed)
        ll cnt = 0;               // occurrences (need to compute)
    };
    vector<State> st = {State{0, -1, {}}};
    int last = 0;
    
    void extend(char c) {
        int cur = sz(st);
        st.emplace_back();
        st[cur].len = st[last].len + 1;
        st[cur].cnt = 1;
        
        int p = last;
        while (p != -1 && !st[p].next.count(c)) {
            st[p].next[c] = cur;
            p = st[p].link;
        }
        if (p == -1) {
            st[cur].link = 0;
        } else {
            int q = st[p].next[c];
            if (st[q].len == st[p].len + 1) {
                st[cur].link = q;
            } else {
                int clone = sz(st);
                st.push_back(st[q]); // copy
                st[clone].len = st[p].len + 1;
                st[clone].cnt = 0;
                while (p != -1 && st[p].next[c] == q) {
                    st[p].next[c] = clone;
                    p = st[p].link;
                }
                st[q].link = st[cur].link = clone;
            }
        }
        last = cur;
    }
    
    // Compute occurrence count for each state (topological order by len)
    void build_cnt() {
        int n = sz(st);
        vector<int> ord(n);
        iota(all(ord), 0);
        sort(all(ord), [&](int a, int b) { return st[a].len > st[b].len; });
        for (int v : ord)
            if (st[v].link != -1)
                st[st[v].link].cnt += st[v].cnt;
    }
    
    // Count distinct substrings: sum(len[v] - len[link[v]])
    ll distinct_substrings() {
        ll ans = 0;
        for (int i = 1; i < sz(st); i++)
            ans += st[i].len - st[st[i].link].len;
        return ans;
    }
};
```

---

## 9. Longest Common Subsequence (LCS) — DP

```cpp
int lcs(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<vector<int>> dp(n+1, vector<int>(m+1));
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dp[i][j] = max({dp[i-1][j], dp[i][j-1], 
                            dp[i-1][j-1] + (a[i-1] == b[j-1])});
    return dp[n][m];
}

// Space-optimized: keep only 2 rows
int lcs_opt(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<int> pre(m+1), cur(m+1);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            cur[j] = max({pre[j], cur[j-1], pre[j-1] + (a[i-1] == b[j-1])});
        }
        pre.swap(cur);
    }
    return pre[m];
}
```

---

## 10. Lexicographically Minimal Rotation (Booth's Algorithm)

```cpp
// Returns starting index of minimal rotation
int minimal_rotation(const string &s) {
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
// Splits string into Lyndon words (lexicographically minimal among rotations)
vector<string> lyndon_factorize(const string &s) {
    vector<string> res;
    int n = sz(s), i = 0;
    while (i < n) {
        int j = i, k = i + 1;
        while (k < n && s[j] <= s[k]) {
            if (s[j] < s[k]) j = i;
            else             j++;
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
    vector<vector<int>> dp(n+1, vector<int>(m+1));
    for (int i = 0; i <= n; i++) dp[i][0] = i;
    for (int j = 0; j <= m; j++) dp[0][j] = j;
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            dp[i][j] = min({
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
                dp[i-1][j-1] + (a[i-1] != b[j-1])
            });
    return dp[n][m];
}

// Space optimized
int edit_distance_opt(const string &a, const string &b) {
    int n = sz(a), m = sz(b);
    vector<int> pre(m+1), cur(m+1);
    iota(all(pre), 0);
    for (int i = 1; i <= n; i++) {
        cur[0] = i;
        for (int j = 1; j <= m; j++)
            cur[j] = min({pre[j]+1, cur[j-1]+1, pre[j-1] + (a[i-1]!=b[j-1])});
        pre.swap(cur);
    }
    return pre[m];
}
```

---

## 13. Palindrome Tree / Eertree — O(n)

Compact representation of all palindromic substrings:

```cpp
struct PalindromeTree {
    struct Node {
        int len, link, occ;
        map<char,int> nxt; // use int nxt[26] for speed
    };
    string s;
    vector<Node> t;
    int last;
    
    PalindromeTree() {
        t.push_back({-1, 0, 0, {}}); // root for odd
        t.push_back({0, 0, 0, {}});  // root for even
        last = 1;
    }
    
    void extend(int pos) {
        char c = s[pos];
        int cur = last;
        while (pos - t[cur].len - 1 < 0 || s[pos - t[cur].len - 1] != c)
            cur = t[cur].link;
        if (t[cur].nxt.count(c)) {
            last = t[cur].nxt[c];
            t[last].occ++;
            return;
        }
        int now = sz(t);
        t.push_back({t[cur].len + 2, 0, 1, {}});
        t[cur].nxt[c] = now;
        
        if (t[now].len == 1) {
            t[now].link = 1; // empty string
        } else {
            int link = t[cur].link;
            while (pos - t[link].len - 1 < 0 || s[pos - t[link].len - 1] != c)
                link = t[link].link;
            t[now].link = t[link].nxt[c];
        }
        last = now;
    }
    
    void build(const string &_s) {
        s = _s;
        for (int i = 0; i < sz(s); i++) extend(i);
        // Accumulate occurrences (process in reverse order of creation)
        for (int i = sz(t)-1; i >= 0; i--)
            t[t[i].link].occ += t[i].occ;
    }
    
    // Number of distinct palindromes = sz(t) - 2
    int distinct_pal() { return sz(t) - 2; }
};
```

---

## 14. Common String Utilities

```cpp
// Check if t is a subsequence of s
bool is_subsequence(const string &t, const string &s) {
    int j = 0;
    for (int i = 0; i < sz(s) && j < sz(t); i++)
        if (s[i] == t[j]) j++;
    return j == sz(t);
}

// Longest Common Prefix of two suffixes (using hashing, binary search)
int longest_common_prefix(const StringHash &sh, int i, int j) {
    int lo = 0, hi = min(sh.n - i, sh.n - j);
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (sh.get(i, i+mid) == sh.get(j, j+mid)) lo = mid;
        else hi = mid - 1;
    }
    return lo;
}

// Compare two substrings in O(log n) using binary search + hashing
bool substring_less(const StringHash &sh, int l1, int r1, int l2, int r2) {
    int len1 = r1 - l1, len2 = r2 - l2;
    int lcp = longest_common_prefix(sh, l1, l2);
    if (lcp >= min(len1, len2)) return len1 < len2;
    char c1 = ...; // original string access needed
    char c2 = ...;
    return c1 < c2;
}
```

---

## Complexity Summary

| Algorithm | Build Time | Query/Match | Memory |
|-----------|-----------|-------------|--------|
| Hashing | O(n) | O(1) | O(n) |
| KMP | O(n) | O(n+m) | O(n) |
| Z-Function | O(n) | O(n+m) | O(n) |
| Trie | O(sum\|s\|) | O(\|s\|) | O(sum\|s\| * alphabet) |
| Aho-Corasick | O(sum\|p\|) | O(\|t\| + #matches) | O(sum\|p\| * alphabet) |
| Manacher | O(n) | — | O(n) |
| Suffix Array | O(n log n) | O(log n) with RMQ | O(n) |
| SAM | O(n) | O(\|s\|) | O(n * alphabet) |
| Palindrome Tree | O(n) | — | O(n * alphabet) |

---

## When to Use What

| Problem type | Recommended |
|-------------|-------------|
| Substring equality / lexicographic compare | Rolling hash |
| Pattern matching in text | KMP or Z-function |
| Multiple pattern matching | Aho-Corasick |
| Palindrome queries | Manacher or Palindrome Tree |
| Longest common substring | Suffix Array + LCP |
| Distinct substrings count | SAM or Suffix Array |
| Lexicographically K-th substring | SAM (DP on DAG) |
| Shortest non-appearing substring | SAM (BFS on transitions) |
| Substring occurrences count | SAM (cnt on states) |
| Edit distance (small strings) | DP O(n*m) |
