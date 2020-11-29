#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

#define mp make_pair
#define pb push_back
#define ll long long

#define sz (1 << 20)
#define maxN (sz + 11)
#define mod 1000000007

ll n, i, x;
ll data[maxN];

int cnt;

void moddo(ll &x) {
    while (x >= mod)
        x -= mod;
}

void to_transform(ll dim, ll *data)
{
    ll len, i, j;
    ll u, v;

    for (len = 1; 2 * len <= dim; len <<= 1) 
    {
        for (i = 0; i < dim; i += 2 * len)
         {
            for (j = 0; j < len; j++) 
            {
                u = data[i + j];
                v = data[i +len + j];

                data[i + j] = v;
                data[i + len + j] = (u + v);

                moddo(data[i + len + j]);
            }
        }
    }
}

void inv_transform(ll dim, ll *data) {
    ll len, i, j;
    ll u, v;

    for (len = 1; 2 * len <= dim; len <<= 1) {
        for (i = 0; i < dim; i += 2 * len) {
            for (j = 0; j < len; j++) {
                u = data[i + j];
                v = data[i + len + j];

                data[i + j] = mod - u + v;
                data[i + len + j] = u;

                moddo(data[i + j]);
            }
        }
    }
}

ll poww(ll a, ll b) {
    ll ans = 1;

    while (b > 0) {
        if (b & 1) ans = (ans * a) % mod;
        b >>= 1;
        a = (a * a) % mod;
    }

    return ans;
}



int main()
{   cout<<"enter size of array"<<endl;

    scanf("%lld", &n);

    cout<<"enter array elements"<<endl;
    for (i = 1; i <= n; i++) {
        scanf("%lld", &x);
        data[x]++;
    }

    data[0]++;

    to_transform(sz, data);
    for (i = 0; i < sz; i++) data[i] = poww(data[i], n);
    inv_transform(sz, data);

    for (i = 0; i < sz; i++)
        if (data[i] != 0)
            cnt++;
    cout<<"required output"<<endl;
    printf("%d\n", cnt);

    return 0;
}
