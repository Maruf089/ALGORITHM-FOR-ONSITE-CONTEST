/// https://vjudge.net/problem/HDU-3948
/// Palindromic tree
/// Ovarall Complexity  O(|S|), মানে  Linear

#include<bits/stdc++.h>
using namespace std;

const int N = 1e5+10;
int ans[N];
int tree[N][26] , idx ;
int len[N] , link[N] , t ;
int occurrence[N];
/// char s[N] ; /// 1-indexed
string s ;

void Init()
{
    memset(tree,0,sizeof tree);
    memset(link,0,sizeof link);
    memset(len,0,sizeof len);

    len[1] = -1 , link[1] = 1 ;
    len[2] = 0 , link[2] = 1 ;
    idx = t = 2 ;
}

void extend(int p)
{
    while(s[p-len[t]-1] != s[p] ) t = link[t] ;
    int x = link[t] , c = s[p] - 'a' ;
    while(s[p-len[x]-1]!=s[p]) x = link[x] ;

    if(!tree[t][c])
    {
        tree[t][c] = ++idx;
        len[idx] = len[t] + 2 ;
        link[idx] = len[idx] == 1 ? 2 : tree[x][c];
        ans[idx] = 1 + ans[link[idx]];
    }
    t = tree[t][c];
    occurrence[t]++;
}

int main()
{
        long long counT = 0 ;
        cin >> s ;
        Init();
        s = '#'+s;
        /// 1-indexed
        for(int j=1;j<s.size();j++)
        {
            extend(j);
            counT += ans[t] ;
            cout << ans[t] << endl;
        }
        cout << counT << endl;
}

