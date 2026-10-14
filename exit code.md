# C++ 程序异常与退出码速查（XCPC / OJ）

> 适用于 Linux（常见 OJ）、Windows 以及本地 GCC / Clang 调试。

## 1. 基本概念

| 术语 | 英文 | 含义 |
|---|---|---|
| 退出码 | Exit Code / Exit Status | 进程结束后向父进程报告的状态值 |
| 异常信号 | Signal | 操作系统发给进程的信号，例如 `SIGSEGV`、`SIGABRT` |
| 运行时错误 | Runtime Error（RE） | OJ 的判题结果，可能由多种异常引起 |
| 编译错误 | Compilation Error（CE） | 编译或链接阶段失败，程序尚未开始正常运行 |

- **退出码就是 Exit Code**，但**退出码不等于信号编号**。
- 习惯上，`Exit Code 0` 表示正常成功退出；`Exit Code 1` 等非零值的具体含义由程序决定。
- 在 Linux 的 Bash 等 shell 中，程序被信号终止时，常见退出状态是 `128 + 信号编号`。例如 `SIGSEGV = 11`，常显示 `Exit Code 139`。
- `Exit Code 139` **不一定**表示段错误：程序也可以主动 `exit(139)`。要结合 stderr、OJ 信息或调试器判断。
- 不同操作系统、shell、IDE、OJ 对状态码的展示方式可能不同。

## 2. Linux 常见异常信号与退出码

> 以下退出码按常见 Linux / Bash 的 `128 + 信号编号` 计算；部分信号编号依赖架构。

| 信号 | 常见退出码 | 含义 | 常见原因 |
|---|---:|---|---|
| `SIGSEGV` | **139** | 段错误 | 数组越界、空指针或野指针解引用、递归栈溢出 |
| `SIGABRT` | **134** | 异常中止 | `assert` 失败、未捕获 C++ 异常导致 `std::terminate`、堆错误 |
| `SIGFPE` | **136** | 算术异常 | 整数除以零、整数模零 |
| `SIGILL` | **132** | 非法指令 | CPU 不支持所使用的指令、执行到无效指令 |
| `SIGBUS` | **135** | 总线错误 | 某些未对齐访问、非法内存映射访问 |
| `SIGKILL` | **137** | 强制杀死 | OJ/系统强制终止、某些内存超限情况 |
| `SIGTERM` | **143** | 终止请求 | IDE、脚本或外部进程请求结束程序 |
| `SIGPIPE` | **141** | 管道写入错误 | 往读端已关闭的管道写数据 |
| `SIGXCPU` | **152** | CPU 时间限制 | 超出系统设置的 CPU 时间限额 |
| `SIGXFSZ` | **153** | 文件大小限制 | 文件写入量超过系统限制 |

**注意：** `SIGFPE` 虽然名字包含 *Floating Point*，在算法竞赛中常由**整数除零**触发。浮点运算 `1.0 / 0.0` 通常产生 `inf`，不一定抛出信号。

## 3. Windows 常见异常代码

> Windows 的异常状态码常以十六进制表示；部分 IDE 将其按 32 位**有符号整数**显示为负数。

| 十进制（有符号） | 十六进制 | 典型含义 |
|---:|---|---|
| **-1073741819** | `0xC0000005` | 访问冲突（非法内存读写或执行） |
| **-1073741571** | `0xC00000FD` | 栈溢出 |
| **-1073741676** | `0xC0000094` | 整数除零 |
| -1073741795 | `0xC000001D` | 非法指令 |
| -1073740791 | `0xC0000409` | 安全检查快速失败（可包括栈缓冲区越界检测） |
| -1073741675 | `0xC0000095` | 整数溢出异常 |
| -1073741801 | `0xC0000017` | 虚拟内存资源不足 |

- **`0xC0000005`：** 先查指针、下标、越界，但不能断定一定是数组越界。
- **`0xC00000FD`：** 先查递归深度、大局部数组、过深调用链。
- **`0xC0000409`：** 不只对应“栈缓冲区越界”，还可能由其他快速失败机制触发。
- **`0xC0000095`：** 普通 C++ 有符号整数溢出是未定义行为，**不保证**自动触发该异常。

## 4. OJ 常见判题结果

| 缩写 | 全称 | 含义 |
|---|---|---|
| `AC` | Accepted | 答案正确 |
| `WA` | Wrong Answer | 答案错误 |
| `TLE` | Time Limit Exceeded | 时间超限 |
| `MLE` | Memory Limit Exceeded | 内存超限 |
| `RE` | Runtime Error | 运行时错误 |
| `CE` | Compilation Error | 编译或链接错误 |
| `OLE` | Output Limit Exceeded | 输出超限 |
| `PE` | Presentation Error | 输出格式错误（部分 OJ 单独判定） |
| `NZEC` | Non-Zero Exit Code | 非零退出码（某些平台的判题说明） |

> `RE`、`TLE`、`MLE`、`OLE` 是 **OJ 判题状态**，不是固定的 OS 退出码。同一种异常在不同 OJ 上也可能以不同判题状态显示。

## 5. C++ 常见错误对应什么现象？

| 错误原因 | 常见表现 | 优先检查 |
|---|---|---|
| 访问数组下标越界 | `SIGSEGV`、访问冲突，或无明显报错 | 下标范围、数组开小、负下标 |
| 解引用空指针 / 野指针 | `SIGSEGV`、访问冲突 | 指针初始化和生命周期 |
| 无限递归 / 递归太深 | 栈溢出、`SIGSEGV` | 终止条件、DFS 链式图深度 |
| 局部数组太大 | 栈溢出、`SIGSEGV` | 将大数组改为全局 / 静态 / 堆内存 |
| 整数除零 / `% 0` | `SIGFPE` 或对应 Windows 异常 | 分母为零的边界情况 |
| `assert(false)` | `SIGABRT`（常见） | 断言为何不成立 |
| `vector::at()` 越界 | `std::out_of_range`，未捕获时常中止 | 容器有效下标 |
| `new` / `vector` 申请内存失败 | `std::bad_alloc`、MLE、被杀死 | 总内存、容器大小 |
| `vector` 尝试创建非法超长长度 | `std::length_error` | `resize()` / `reserve()` 的参数 |
| 执行 CPU 不支持的指令 | `SIGILL` | 编译选项，例如 `-march=native` |
| 死循环或复杂度过高 | TLE、外部终止 | 循环条件、复杂度 |
| 内存超限 | MLE，有时 `SIGKILL` | 数组大小、递归、动态分配 |
| 向关闭的管道写入 | `SIGPIPE` | 对拍、管道通信双方是否提前退出 |

**警告：** C++ 的未定义行为（例如 `a[i]` 越界、带符号整数溢出）**不保证立刻报错**；也可能表现为 WA、结果随机或只在某些优化等级下出问题。

## 6. 常见终端报错提示

| 信息 | 一般含义 |
|---|---|
| `Segmentation fault` | 段错误，非法内存访问 |
| `Floating point exception` | 算术异常，常见为整数除零 |
| `Aborted (core dumped)` | 进程被异常中止 |
| `stack overflow` | 栈空间耗尽 |
| `std::bad_alloc` | 动态内存分配失败 |
| `std::out_of_range` | 使用带边界检查的接口越界 |
| `std::length_error` | 申请的容器长度不合法 / 过大 |
| `undefined reference` | 链接阶段找不到符号定义 |
| `multiple definition` | 链接阶段出现重复定义 |
| `expected ';'` | 编译阶段的语法错误 |
| `exit code 1` | 非零退出，但无法仅据此判断原因 |
| `exit code 127`（Bash） | 通常是命令不存在 / 找不到可执行文件 |
| `Permission denied` | 没有相应文件的访问或执行权限 |

## 7. 本地调试方法（GCC / Clang）

### 7.1 AddressSanitizer + UndefinedBehaviorSanitizer

```bash
g++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer main.cpp -o main
./main
```

有助于查出：堆/栈的部分越界、释放后访问、整数未定义行为等。检测范围并非覆盖所有错误。

### 7.2 libstdc++ 调试模式（GCC 标准库）

```bash
g++ -std=c++17 -O0 -g -D_GLIBCXX_DEBUG main.cpp -o main
./main
```

有助于检查部分 STL 容器 / 迭代器的非法操作（例如越界、迭代器失效）。仅在本地调试使用：可能显著变慢，且会改变部分 STL 类型的 ABI。

### 7.3 查看程序退出码

Linux / macOS 终端：

```bash
./main
echo $?
```

说明：`$?` 是 shell 上一条命令的退出状态。若发生信号终止，应结合终端提示一起判断。

## 8. XCPC 赛时速查（建议优先记住）

| 看到的退出码 | 优先怀疑 |
|---|---|
| **139** | Linux `SIGSEGV`：越界、非法指针、栈溢出 |
| **134** | Linux `SIGABRT`：断言、未捕获异常 |
| **136** | Linux `SIGFPE`：整数除以零 / `% 0` |
| **132** | Linux `SIGILL`：非法 CPU 指令 |
| **137** | Linux `SIGKILL`：被强制终止，可能是 MLE |
| **143** | Linux `SIGTERM`：被请求终止（例如 IDE 停止运行） |
| **-1073741819** | Windows `0xC0000005`：访问冲突 |
| **-1073741571** | Windows `0xC00000FD`：栈溢出 |

**排查顺序：** 看 OJ 状态 → 看终端 / stderr 报错 → 对照信号或异常码 → 本地开 Sanitizer / 调试模式 → 检查边界、递归、分母、内存。
