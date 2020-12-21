#include<bits/stdc++.h>
using namespace std ;
///*

/// largest prefix that matches a suffix
int kmp(string s)
{
    string t = s + "#";
    reverse(s.begin(), s.end());
    t += s;
    s = t ;
    vector<int> lps((int)s.size(), 0);
    for (int i = 1; i < (int)lps.size(); i++) {
        int prv_idx = lps[i - 1];
        while (prv_idx > 0 && s[i] != s[prv_idx]) {
            prv_idx = lps[prv_idx - 1];
        }
        lps[i] = prv_idx + (s[i] == s[prv_idx] ? 1 : 0);
    }
    return lps[(int)lps.size() - 1];
}

void createlps()
{
    for(int i = 0; i < str.size(); i++) lps.push_back(0);
    int index = 0, i = 1;

    for(; i < str.size();){
        if(str[i] == str[index]){
            lps[i] = index+1;
            i++, index++;
        }
        else{
            if(index != 0){
                index = lps[index-1];
            }
            else{
                lps[i] = 0;
                i++;
            }
        }
    }
}


void kmp_process(string text,string pattern)
{
    int n = text.size();
    int m = pattern.size();
    prefix_func(pattern);
    int now = -1;
    for(j=0;j<n;j++)
    {
      while(now!=-1 and pattern[now+1]!=text[j])
            now = phi[now];
      if(pattern[now+1]==text[j])
        now++;
      if(now==m-1)
      {
          cout << j-(m-1) << endl; /// found pattern
          now = phi[now] ;
      }
    }
}

int main()
{
    fast;
    while(cin >> n )
    {
        cin >> pattern >> text ;
        kmp_process(text,pattern);
        cout << ln ;
    }
}

*///

void computeLPSArray(char* pat, int M, int* lps)
{
    int len = 0;
    lps[0] = 0; /// lps[0] is always 0
    int i = 1;
    while (i < M)
    {
        if (pat[i] == pat[len])
        {
            len++;
            lps[i] = len;
            i++;
        }
        else /// (pat[i] != pat[len])
        {
            if (len != 0)
            {
                len = lps[len - 1];
            }
            else /// if (len == 0)
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}
void KMPSearch(char* pat, char* txt)
{
    int M = strlen(pat);
    int N = strlen(txt);
    int lps[M];
    computeLPSArray(pat, M, lps);

    int i = 0; /// index for txt[]
    int j = 0; /// index for pat[]
    while (i < N)
    {
        if (pat[j] == txt[i])
        {
            j++;
            i++;
        }
        if (j == M)   /// printf("Found pattern at index %d ", i - j);
        {
            j = lps[j - 1];
            ///      counT++; count for number of pattern match
        }
        else if (i < N && pat[j] != txt[i])   /// mismatch after j matches
        {
            if (j != 0)
                j = lps[j - 1];
            else
                i = i + 1;
        }
    }
}
void cnounT_number_of_occurrences_of_each_preffix()
{
    int ans = s.length();
    for(int i = 1; i < s.length() ; i ++ )
    {
        int j = pi[i];
        while(j > 0 )
        {
            ans ++;
            j = pi[j-1];
            ans %= MOD;
        }
    }
}

int kmp_process(string text,string pattern) ///  search for the longest prefix of pattern in text.
{
    int j=0;
    int n = text.size();
    int m = pattern.size();
    vector<int>preffix = prefix_function(pattern);
    for(int i=0;i<n;i++)
    {
        if(text[i]==pattern[j])
            j++;
        else
        {
            while(j>0)
            {
                j = preffix[j-1];
                if(text[i]==pattern[j])
                {
                    j++;
                    break;
                }
            }
        }
    }
    return j;
}

int main()
{

}

/// https://rezwanarefin01.github.io/posts/palindromic-tree-01/
/// https://vjudge.net/problem/HDU-3948
/// Palindromic tree distinct substring
#include<bits/stdc++.h>
using namespace std;

const int N = 1e5+10;
int ans[N];
int tree[N][26], idx ;
int len[N], link[N], t ;
int occurrence[N];
/// char s[N] ; /// 1-indexed
string s ;

void Init()
{
    memset(tree,0,sizeof tree);
    memset(link,0,sizeof link);
    memset(len,0,sizeof len);

    len[1] = -1, link[1] = 1 ;
    len[2] = 0, link[2] = 1 ;
    idx = t = 2 ;
}

void extend(int p)
{
    while(s[p-len[t]-1] != s[p] )
        t = link[t] ;
    int x = link[t], c = s[p] - 'a' ;
    while(s[p-len[x]-1]!=s[p])
        x = link[x] ;

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
    int tt ;
    cin >> tt;
    for(int i=1;i<=tt;i++)
    {
        cin >> s ;
        Init();
        int counT = 0 ;
        s = '#'+s;
        for(int j=1;j<s.size();j++)
        {
            extend(j);
          //  counT += ans[t] ;
        }
        printf("Case #%d: %d\n",i,idx-2); /// Distinct palindrome
  ///  cout << counT << endl ; /// Not Distinct

    for(int i=idx;i>2;i--)
    {
        occurrence[link[i]] += occurrence[i];
    }
    for(int i=3;i<=idx;i++)
        cout << occurrence[i] << " ";
    }
}



/// For Hashing
use double hash for less collision , use mod value 10^9+7 , 10^9+11 .
 0 1 2 3 4 5 .....
 if we want to get hash for (3 4) -> Hash = prefix_hash(4) - prefix_hash(2)*base^(koto length nibo) : -> ekhane koto length nibo = 2
