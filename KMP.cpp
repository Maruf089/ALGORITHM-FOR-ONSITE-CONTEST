/// Complexity O(n+m)

#include<bits/stdc++.h>
using namespace std;

vector<int>construct_lps(string pattern)
{
    vector<int>lps(pattern.size());
    int index = 0;
    for(int i=1; i<pattern.size();)
    {
        if(pattern[i]==pattern[index])
        {
            lps[i] = index + 1;
            index++, i++;
        }
        else
        {
            if(index!=0)
                index = lps[index-1];
            else
                lps[i] = index, i++;
        }
    }
    return lps;
}

void KMP(string text,string pattern)
{
    bool found = false;
    vector<int>lps = construct_lps(pattern);
    int j = 0, i = 0 ;
    /// i --> text, j --> pattern
    while(i<text.size())
    {
        if(text[i]==pattern[j])
            i++, j++;
        else
        {
            if(j!=0)
                j = lps[j-1];
            else
                i++;
        }

        if(j==pattern.size())
        {
            found = true;
            cout << "Found at index : " << (i-pattern.size()) << endl;
            j = lps[j-1];
        }
    }
    if(!found)
        cout << "Not Found \n" ;
}


int main()
{
    string text,pattern;
    cin >> text >> pattern;

    KMP(text,pattern);
    return 0;
}


/// Counting the number of occurrences of each prefix

vector<int>LPS = construct_lps(s);
vector<int>occur(n+1);

for(int i=0; i<n; i++)
    occur[LPS[i]]++;
for(i=n-1; i>0; i--)
    occur[LPS[i-1]] += occur[i];
for(i=0; i<=n; i++)
    occur[i]++;

 /// Counting the number of occurrences of each prefix also it is a suffix

 void LPS()
{
    int j=0;
    nxt[0]=nxt[1]=0;
    for(int i=2; i<=n; i++)
    {
        while(j && s[i]!=s[j+1])
            j=nxt[j];
        if(s[i]==s[j+1])
            j++;
        nxt[i]=j;
    }
}
int main()
{
    scanf("%s",s+1);
    n = strlen(s+1);
    LPS();

    for(int i=1; i<=n; i++)
        dp[i]=1;
    for(int i=n; i; i=nxt[i])
        ans[++tot]=i;
    for(int i=n; i; i--)
        dp[nxt[i]] += dp[i];

    printf("%d\n",tot);
    sort(ans+1,ans+tot+1);
    for(int i=1; i<=tot; i++)
        printf("%d %d\n",ans[i],dp[ans[i]]);
}
