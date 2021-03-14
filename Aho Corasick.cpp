#include<bits/stdc++.h>
using namespace std;

const int N = 1e6 + 6;
const int M = 505;
const int nx = 2e5 + 5;

struct node
{
    int link[27];
    node()
    {
        memset(link, -1, sizeof link);
    }
};

node a[N];//2e5
char s[N];//1e6
char ch[M];// 505
int suffix[N];//2e5
int val[N];//2e5
int path[N]; //2e5
int en[M]; //505
int num;
int koto = 0;

void init()
{
    memset(a, -1, sizeof a);
    memset(suffix, 0, sizeof suffix);
    memset(val, 0, sizeof val);
    memset(path, 0, sizeof path);
    num = 0;
    koto = 0;
}

void makeTrie(int j, int l)
{
    int cur = 0;
    for (int i = 1; i <= l; i++)
    {
        int x = ch[i] - 96;
        if (a[cur].link[x] == (-1))
        {
            a[cur].link[x] = ++num;
        }
        cur = a[cur].link[x];
    }
    en[j] = cur;
}

void bfs()
{
    queue<int> Q;
    for (int i = 1; i <= 26; i++)
    {
        if (a[0].link[i] == -1)
        {
            a[0].link[i] = 0;
        }
        else
        {
            Q.push(a[0].link[i]);
            suffix[a[0].link[i]] = 0;
        }
    }

    while (!Q.empty())
    {
        int u = Q.front();
        Q.pop();
        for (int i = 1; i <= 26; i++)
        {
            int v = a[u].link[i];
            if (v == -1)
            {
                a[u].link[i] = a[suffix[u]].link[i];
            }
            else
            {
                Q.push(v);
                suffix[v] = a[suffix[u]].link[i];
                path[++koto] = v;
            }
        }
    }
}

void search(int l)
{
    int cur = 0;
    for (int i = 1; i <= l; i++)
    {
        int x = s[i] - 96;
        cur = a[cur].link[x];
        val[cur]++;
    }
    for (int i = koto; i >= 1; i--)
        val[suffix[path[i]]] += val[path[i]];
}



int main()
{

    int T;
    scanf("%d", &T);
    for (int cs = 1; cs <= T; cs++)
    {
        init();
        int n;
        scanf("%d", &n);
        scanf("%s", s + 1);
        for (int i = 1; i <= n; i++)
        {
            scanf("%s", ch + 1);
            int l = strlen(ch + 1);
            makeTrie(i, l);
        }

        bfs();
        int l = strlen(s + 1);
        search(l);
        printf("Case %d:\n", cs);
        for (int i = 1; i <= n; i++)
        {
            printf("%d\n", val[en[i]]);
        }
    }
    return 0;
}
