### xor 线性基

统一约定：

- `LOG = 60`，表示 bit `0 ~ 59`
- `p[i]` 表示最高位为 `i` 的基向量
- 统一从高位向低位插入
- 本质上是 \(GF(2)\) 上的高斯消元

---

#### 1. 基础线性基

- 判断线性相关
- 求线性基的秩
- 求最大 xor

```cpp
// 线性基 (60 位, 适用于 1e18)
struct XORBasis {
    static const int LOG = 60;
    long long p[LOG] = {};

    bool insert(long long x) {
        for(int i = LOG - 1; ~i; i--) {
            if(!(x >> i & 1)) continue;

            if(!p[i]) {
                p[i] = x;
                return true;
            }

            x ^= p[i];
        }

        return false;
    }

    long long maxXor(long long init = 0) {
        long long res = init;

        for(int i = LOG - 1; ~i; i--)
            if((res ^ p[i]) > res)
                res ^= p[i];

        return res;
    }

    int rank() {
        int res = 0;

        for(int i = 0; i < LOG; i++)
            res += !!p[i];

        return res;
    }
};
```

---

#### 2. 第 \(k\) 小异或值

- 普通线性基插入后是行阶梯形
- 先通过 `rebuild()` 转为行最简形
- 此时每个主元位只在对应基向量中出现
- 按主元从低到高枚举即可得到第 \(k\) 小

行最简形类似：

```text
        5 4 3 2 1 0

p[5] = 1 x 0 x 0 x
p[3] = 0 0 1 x 0 x
p[1] = 0 0 0 0 1 x
```

主元列构成单位矩阵。

```cpp
struct XORBasis {
    static const int LOG = 60;
    long long p[LOG] = {};

    bool insert(long long x) {
        for(int i = LOG - 1; ~i; i--) {
            if(!(x >> i & 1)) continue;

            if(!p[i]) {
                p[i] = x;
                return true;
            }

            x ^= p[i];
        }

        return false;
    }

    void rebuild() {
        for(int i = LOG - 1; ~i; i--) {
            if(!p[i]) continue;

            for(int j = i - 1; ~j; j--)
                if((p[i] >> j & 1) && p[j])
                    p[i] ^= p[j];
        }
    }

    // 第 k 小不同异或值
    // 需先 rebuild()
    // k 从 1 开始, 0 为第 1 小
    long long kth(long long k) {
        k--;

        long long res = 0;

        for(int i = 0; i < LOG; i++) {
            if(!p[i]) continue;

            if(k & 1)
                res ^= p[i];

            k >>= 1;
        }

        return res;
    }
};
```

若线性基秩为 \(r\)，不同异或值一共有：

\[
2^r
\]

个。

---

#### 3. 在线维护行最简型

边插入边维护行最简形：

- `p[i]` 的最高位为 `i`
- 每个主元位只在对应 `p[i]` 中出现

```cpp
struct ReducedBasis {
    static const int LOG = 60;
    long long p[LOG] = {};

    bool insert(long long x) {
        for(int i = LOG - 1; ~i; i--) {
            if(!(x >> i & 1)) continue;

            if(p[i]) {
                x ^= p[i];
                continue;
            }

            // 消掉新行中的已有低位主元
            for(int j = i - 1; ~j; j--)
                if((x >> j & 1) && p[j])
                    x ^= p[j];

            p[i] = x;

            // 消掉已有高位行中的新主元
            for(int j = i + 1; j < LOG; j++)
                if(p[j] >> i & 1)
                    p[j] ^= p[i];

            return true;
        }

        return false;
    }
};
```

---

### 图论 xor

#### 1. 带环图路径 xor

先任选一棵生成树。

定义：

\[
dis[u]
\]

为生成树上根到 \(u\) 的路径 xor。

则生成树上 \(u \to v\) 的路径 xor 为：

\[
dis[u] \oplus dis[v]
\]

对于一条非树边：

\[
(u,v,w)
\]

它对应的环 xor 为：

\[
dis[u] \oplus dis[v] \oplus w
\]

将所有环 xor 插入线性基。

于是任意 \(u \to v\) 路径 xor 都可以写成：

\[
dis[u] \oplus dis[v] \oplus \text{环空间中的任意元素}
\]

因此最大路径 xor：

```cpp
cycleBase.maxXor(dis[u] ^ dis[v]);
```

模板：

```cpp
struct Edge {
    int v, id;
    long long w;
};

vector<Edge> g[N];

XORBasis cycleBase;
long long dis[N];
bool vis[N];

void dfs(int u, int pe = -1) {
    vis[u] = true;

    for(auto [v, id, w] : g[u]) {
        if(id == pe) continue;

        if(!vis[v]) {
            dis[v] = dis[u] ^ w;
            dfs(v, id);
        }
        else {
            cycleBase.insert(dis[u] ^ dis[v] ^ w);
        }
    }
}

// 查询:
// cycleBase.maxXor(dis[u] ^ dis[v]);
```

如果保证无重边，也可以写：

```cpp
void dfs(int u, int fa) {
    vis[u] = true;

    for(auto [v, w] : g[u]) {
        if(v == fa) continue;

        if(!vis[v]) {
            dis[v] = dis[u] ^ w;
            dfs(v, u);
        }
        else {
            cycleBase.insert(dis[u] ^ dis[v] ^ w);
        }
    }
}
```

---

### 区间 / 前缀线性基

#### 1. 前缀线性基

每个右端点 \(r\) 维护一个线性基。

额外维护：

```cpp
pos[i]
```

表示当前 `b[i]` 对应基向量中，能够保证使用的最靠右位置。

插入时优先保留位置更靠右的向量。

这样 `pref[r]` 可以查询区间：

\[
[l,r]
\]

中任取若干元素的最大 xor。

```cpp
struct PrefixBasis {
    static const int LOG = 60;

    long long b[LOG] = {};
    int pos[LOG];

    PrefixBasis() {
        fill(pos, pos + LOG, -1);
    }

    void insert(long long x, int id) {
        for(int i = LOG - 1; ~i; i--) {
            if(!(x >> i & 1)) continue;

            if(!b[i]) {
                b[i] = x;
                pos[i] = id;
                return;
            }

            if(id > pos[i]) {
                swap(x, b[i]);
                swap(id, pos[i]);
            }

            x ^= b[i];
        }
    }

    // 当前版本对应右端点 r
    // 查询 [l, r] 中任取若干数的最大 xor
    long long maxXor(int l, long long init = 0) {
        long long res = init;

        for(int i = LOG - 1; ~i; i--)
            if(pos[i] >= l && (res ^ b[i]) > res)
                res ^= b[i];

        return res;
    }
};
```

使用：

```cpp
vector<PrefixBasis> pref(n);

for(int i = 0; i < n; i++) {
    if(i)
        pref[i] = pref[i - 1];

    pref[i].insert(a[i], i);
}

// 查询 a[l..r] 中任取若干元素的最大 xor
long long ans = pref[r].maxXor(l);
```

---

### 异或高斯消元

异或线性基本质上就是 \(GF(2)\) 上的高斯消元。

对应关系：

```text
线性基                  高斯消元

p[i]                 <-> 矩阵的一行
x ^= p[i]            <-> Ri ^= Rj
最高位               <-> 主元
insert               <-> 前向高斯消元
rebuild              <-> Gauss-Jordan 回消
线性基大小           <-> rank
```

---

#### 1. 增广矩阵

求解：

\[
Ax=b
\]

其中系数矩阵中的元素为 `0 / 1`。

```cpp
constexpr int N = 1010;

bitset<N> a[N];
long long b[N];

// 返回:
//  0 : 唯一解
//  1 : 多解
// -1 : 无解
int gauss(int n, int m, vector<long long>& ans) {
    int r = 0;

    for(int c = 0; c < m && r < n; c++) {
        int p = r;

        while(p < n && !a[p][c])
            p++;

        if(p == n) continue;

        swap(a[p], a[r]);
        swap(b[p], b[r]);

        // 消掉该主元列的其他所有行
        // 最终得到行最简型
        for(int i = 0; i < n; i++) {
            if(i != r && a[i][c]) {
                a[i] ^= a[r];
                b[i] ^= b[r];
            }
        }

        r++;
    }

    // 检查 0 = b[i]
    for(int i = r; i < n; i++)
        if(b[i])
            return -1;

    // 存在自由元
    if(r < m)
        return 1;

    // 唯一解
    ans.assign(b, b + m);

    return 0;
}
```

唯一解时最终类似：

```text
1 0 0 0 | b[0]
0 1 0 0 | b[1]
0 0 1 0 | b[2]
0 0 0 1 | b[3]
```

因此：

```cpp
ans[i] = b[i];
```

---

#### 2. 系数矩阵：行阶梯型

只消掉主元下面的元素。

```cpp
bitset<N> a[N];

// n 行 m 列
// 返回 rank
int gaussXor(int n, int m) {
    int r = 0;

    for(int c = 0; c < m && r < n; c++) {
        int p = r;

        while(p < n && !a[p][c])
            p++;

        if(p == n) continue;

        swap(a[p], a[r]);

        for(int i = r + 1; i < n; i++)
            if(a[i][c])
                a[i] ^= a[r];

        r++;
    }

    return r;
}
```

结果类似：

```text
1 x x x x x
0 1 x x x x
0 0 0 1 x x
0 0 0 0 0 0
```

---

#### 3. 系数矩阵：行最简型

将每个主元列中的其他 `1` 全部消掉。

```cpp
bitset<N> a[N];

// n 行 m 列
// 返回 rank
int gaussXor(int n, int m) {
    int r = 0;

    for(int c = 0; c < m && r < n; c++) {
        int p = r;

        while(p < n && !a[p][c])
            p++;

        if(p == n) continue;

        swap(a[p], a[r]);

        for(int i = 0; i < n; i++)
            if(i != r && a[i][c])
                a[i] ^= a[r];

        r++;
    }

    return r;
}
```

结果类似：

```text
1 x 0 x 0 x
0 0 1 x 0 x
0 0 0 0 1 x
```

主元列构成单位矩阵。

---

#### 4. 记录当前行由哪些原始行异或得到

维护：

```cpp
bitset<N> from[N];
```

初始化：

```cpp
for(int i = 0; i < n; i++)
    from[i][i] = 1;
```

其中：

```cpp
from[i][j] = 1
```

表示当前第 `i` 行中包含原始第 `j` 行。

每次行操作同步维护：

```cpp
swap(a[p], a[r]);
swap(from[p], from[r]);
```

以及：

```cpp
a[i] ^= a[r];
from[i] ^= from[r];
```

完整模板：

```cpp
bitset<N> a[N];
bitset<N> from[N];

int gaussXor(int n, int m) {
    for(int i = 0; i < n; i++)
        from[i][i] = 1;

    int r = 0;

    for(int c = 0; c < m && r < n; c++) {
        int p = r;

        while(p < n && !a[p][c])
            p++;

        if(p == n) continue;

        swap(a[p], a[r]);
        swap(from[p], from[r]);

        for(int i = r + 1; i < n; i++) {
            if(a[i][c]) {
                a[i] ^= a[r];
                from[i] ^= from[r];
            }
        }

        r++;
    }

    return r;
}
```

---

### 总结

统一记忆：

```text
线性基:

p[i] = 最高位为 i 的基向量

insert:
从高位向低位消元
维护行阶梯型

rebuild:
消掉低位主元
得到行最简型

ReducedBasis:
边插入边维护行最简型

PrefixBasis:
额外维护 pos[i]
优先保留位置更靠右的向量

图上 xor:
dis[u] ^ dis[v]
+
环空间线性基

高斯消元:
一行 <-> 一个二进制向量
Ri ^= Rj <-> 异或消元
```
