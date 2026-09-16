### C++ 大数运算类

```
#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using boost::multiprecision::cpp_int;
```

### 对拍

** Linux **
```
#include <bits/stdc++.h>
using namespace std;

signed main()
{
    for(int T = 1; ; T++) {
        system("./gen > in.txt");
        system("./std < in.txt > std.txt");
        system("./bf < in.txt > bf.txt");

        if(system("diff -w std.txt bf.txt > /dev/null")) {
            cout << "WA " << T << '\n';
            break;
        }

        cout << "AC " << T << '\n';
    }
}
```

** windows **
```
#include <bits/stdc++.h>
using namespace std;

signed main()
{
    for(int T = 1; ; T++) {
        system("gen.exe > in.txt");
        system("std.exe < in.txt > std.txt");
        system("bf.exe < in.txt > bf.txt");

        if(system("fc std.txt bf.txt > nul")) {
            cout << "WA " << T << '\n';
            break;
        }

        cout << "AC " << T << '\n';
    }
}
```
