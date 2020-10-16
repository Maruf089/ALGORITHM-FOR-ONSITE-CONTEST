/// Maximum Spanning tree for 1-5 Dimension
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
const int N = 2e5 + 7;

int dsu[N];
int get(int v)
{
    if (v == dsu[v]) return v;
    else return dsu[v] = get(dsu[v]);
}

void uni(int u, int v)
{
    dsu[get(u)] = get(v);
}

int main()
{
    int n, d;
    cin >> n >> d;
    vector <vector <int> > x(n, vector <int> (d));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < d; j++)
            cin >> x[i][j];

    vector <pair <int, pair <int, int> > > e;
    auto add_edge = [&] (int u, int v)
    {
        int sum = 0;
        for (int i = 0; i < d; i++)
            sum += abs(x[u][i] - x[v][i]);
        e.push_back({sum, {u, v}});
    };


    for (int mask = 0; mask < (1 << d); mask++)
    {
        vector <pair <int, int> > e;
        for (int i = 0; i < n; i++)
        {
            int sum = 0;
            for (int j = 0; j < d; j++)
            {
                if ((mask >> j) & 1)
                {
                    sum += x[i][j];
                }
                else
                {
                    sum -= x[i][j];
                }
            }
            e.push_back({sum, i});
        }

        int s = min_element(e.begin(), e.end())->second;
        int t = max_element(e.begin(), e.end())->second;
        for (int x = 0; x < n; x++)
        {
            add_edge(s, x);
            add_edge(t, x);
        }
    }


    sort(e.rbegin(), e.rend());
    for (int i = 0; i < n; i++)
        dsu[i] = i;
    ll ans = 0;
    for (auto c : e)
    {
        int u = c.second.first, v = c.second.second;
        if (get(u) != get(v))
        {
            uni(u, v);
            ans += c.first;
        }
    }
    cout << ans << '\n';
}
