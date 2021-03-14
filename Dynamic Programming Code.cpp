#include<bits/stdc++.h>
using namespace std;



/// LCS + LIS ( also with print ) for duplicate value

    f1(i,n)
    cin >> arr[i];
    f1(i,m)
    cin >> brr[i];

    for(int i=1; i<=n; i++)
    {
        int last = 0;
        for(int j=1; j<=m; j++)
        {
            if(arr[i]==brr[j])
            {
                    dp[j] = 1+dp[last];
                    par[j] = last;
            }
            if(brr[j]<arr[i] and dp[j]>dp[last])
                last = j;
        }
    }
    int res = 0,idx = -1;
    f0(j,mx)
    {
        if(dp[j]>res)
        {
            res = dp[j];
            idx = j;
        }
    }
    cout << res << ln;
    vector<ll>final_ans;
    while(res--)
    {
       final_ans.pb(brr[idx]);
       idx = par[idx];
    }
    reverse(all(final_ans));
    for(auto it : final_ans)
        cout << it << ' ';
*******///////////////////////////////


int LCSubStr(string X, string Y) /// Longest Common Substring
{
    int m = X.length(),n = Y.length(),result = 0;
    int len[2][n];
    int currRow = 0;
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            if (i == 0 || j == 0)
            {
                len[currRow][j] = 0;
            }
            else if (X[i - 1] == Y[j - 1])
            {
                len[currRow][j] = len[1 - currRow][j - 1] + 1;
                result = max(result, len[currRow][j]);
            }
            else
            {
                len[currRow][j] = 0;
            }
        }
        currRow = 1 - currRow;
    }
    return result;
}

int LCS(string X,string Y, int m, int n )   /// Longest Common Subsequence
{
    int L[m + 1][n + 1];
    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            if (i == 0 || j == 0)
                L[i][j] = 0;
            else if (X[i - 1] == Y[j - 1])
                L[i][j] = L[i - 1][j - 1] + 1;
            else
                L[i][j] = max(L[i - 1][j], L[i][j - 1]);
        }
    }
    return L[m][n];
}

/// Print all LCS

string s,ss ;
set<string>sequences;
map<string,bool>process[86][86];
ll dp[86][86];

void findseq(int i, int j,string word)
{
	if(process[i][j].find(word)!=process[i][j].end())
		return;
	if(i==0||j==0)
	{
		sequences.insert(word);
	}
	else if(s[i-1]==ss[j-1])
	{
		word = s[i-1] + word;
		findseq(i-1,j-1,word);
	}
	else if(dp[i][j-1]>dp[i-1][j])
	{
		findseq(i,j-1,word);
	}
	else if(dp[i-1][j]>dp[i][j-1])
	{
		findseq(i-1,j,word);
	}
	else
	{
		findseq(i,j-1,word);
		findseq(i-1,j,word);
	}
	process[i][j][word] = true;
}

int main()
{
    fast ;
    //  freopen("in.txt","r",stdin);
    cin >> t ;
    f1(cs,t)
    {
        cin >> s >> ss ;

        for(int i=0; i<=85; i++)
            for(int j=0; j<=85; j++)
                dp[i][j] = 0 ;

        f1(i,s.sz)  /// Find LCS
        {
            f1(j,ss.sz)
            {
                if(s[i-1]==ss[j-1])
                {
                    dp[i][j] = dp[i-1][j-1]+1;
                }
                else
                    dp[i][j] = max(dp[i-1][j],dp[i][j-1]);
            }
        }
    //    dbg(dp[s.sz][ss.sz]);

        sequences.clear();

        findseq(s.sz,ss.sz,"");

        for(auto it : sequences)
            cout << it << ln ;
        cout << ln ;

    }
}



long int UnboundedKnapsack(long int Capacity,long int n, long int weight[],long int val[],int W) /// Unbounded Knapsack
{
    long int dp[Capacity+1];
    for(int i=0; i < W+1; i++)
    {
        dp[i]=0;
    }
    for(int i=0; i < W+1; i++)
    {
        for(int j=0; j < n; j++)
        {
            if(weight[j] < i)
            {
                dp[i] = max(dp[i], dp[i-weight[j]] + val[j]);
            }
        }
    }
    return dp[Capacity];
}

/// space optimized Knapsack
ll KnapSack(ll val[],ll wt[],ll n,ll W)
{
    ll dp[W+1];
    memset(dp, 0, sizeof(dp));
    for( i=0; i < n; i++)
        for( j=W; j>=wt[i]; j--)
            dp[j] = max(dp[j] , val[i] + dp[j-wt[i]]);
    return dp[W];
}


int editDistDP(string str1, string str2, int m, int n)
{
	int dp[m+1][n+1];
	for (int i=0; i<=m; i++)
	{
		for (int j=0; j<=n; j++)
		{
			if (i==0)
				dp[i][j] = j; // Min. operations = j
			else if (j==0)
				dp[i][j] = i; // Min. operations = i
			else if (str1[i-1] == str2[j-1])
				dp[i][j] = dp[i-1][j-1];
			else
				dp[i][j] = 1 + min(dp[i][j-1], // Insert
								min(dp[i-1][j], // Remove
								dp[i-1][j-1])); // Replace
		}
	}
	return dp[m][n];
}



ll dp[100000+9][100];
ll wt[mx+9],val[mx+9];
ll W_MAX = (ll)1e9;

ll knapSack(int r,int i = 0 )  /// Knapsack for Large weight , W = 10^9 , n*max(Vi)<=10^6
{
    if(r<=0)
        return 0;
    if(i==n)
        return W_MAX ;
    if(dp[r][i]!=-1)
        return dp[r][i];
    dp[r][i] = min(knapSack(r,i+1),wt[i]+knapSack(r-val[i],i+1));
    return dp[r][i];
}
ll maxWeight()
{
    for(i=sum; i>=0; i--)
        if(knapSack(i)<=W)
            return i ;
    return 0;
}
int main()
{
    cin >> n  >> W;
    MEM(dp,-1);
    f0(i,n)
    {
        cin >> wt[i] >> val[i];
        sum += val[i];
    }
    cout << maxWeight() << ln ;
}


unsigned int nextPowerOf2(unsigned int n)
{
    n--;n |= n >> 1;n |= n >> 2;
    n |= n >> 4;n |= n >> 8;n |= n >> 16;n++;
    return n;
}
int highestPowerof2(int n)
{
   int p = (int)log2(n);
   return (int)pow(2, p);
}
int highestOneBit(int i) {
i |= (i >>  1);i |= (i >>  2);i |= (i >>  4);
    i |= (i >>  8);i |= (i >> 16);
    return i - ((unsigned)i >> 1);
}


///  largest contiguous array sum  ( kadane's algorithm )

int maxSubArraySum(int a[], int size)
{
    int max_so_far = INT_MIN, max_ending_here = 0;

    for (int i = 0; i < size; i++)
    {
        max_ending_here = max_ending_here + a[i];
        if (max_so_far < max_ending_here)
            max_so_far = max_ending_here;

        if (max_ending_here < 0)
            max_ending_here = 0;
    }
    return max_so_far;
}

/// number of Subset sum between A and B ( meet in the middle)

inline bool checkBit(ll n, int i) { return n&(1LL<<i); }
int X[2000005] , Y[2000005];
void calcsubarray(int arr[],int x[],int n,int c)
{
    /** Subset Sum taking 0/1 elements **/
    for(int i=0;i<(1<<n);i++)
    {
        ll s = 0;
        for(int j=0;j<n;j++)
            if(checkBit(i,j))
                s += arr[j+c];
        x[i] = s;
    }
}
int numberOfSubsets(int arr[],int n,int A,int B)
{
    calcsubarray(arr,X,n/2,0);
    calcsubarray(arr,Y,n-n/2,n/2);

    int size_X = 1<<(n/2);
    int size_Y = 1<<(n-n/2);
    sort(Y,Y+size_Y);

    int counT = 0;
    for(int i=0;i<size_X;i++)
    {
        if(X[i]<=B)
        {
           auto low = lower_bound(Y,Y+size_Y,A-X[i]);
           auto high = upper_bound(Y,Y+size_Y,B-X[i]);

           counT += (high-low);
        }
    }
    return counT;
}


///Sum of subsets of all the subsets of an array | O(N)

There are N – 1CK – 1 subsets of size K for every index that include it.
Contribution of an element for a subset of size K will be equal to 2^(K – 1) times its value.
Thus, total contribution of an element for all the subsets of length K will be equal to N – 1CK – 1 * 2^(K – 1)
Total contribution among all the subsets will be equal to:
N – 1CN – 1 * 2(N – 1) + N – 1CN – 2 * 2(N – 2 + N – 1CN – 3 * 2(N – 3) + … + N – 1C0 * 20.

int ncr(int n, int r)
{
    return (fact[n] / fact[r]) / fact[n - r];
}
int findSum(int* arr, int n)
{
    fact[0] = 1;
    for (int i = 1; i < n; i++)
        fact[i] = i * fact[i - 1];
    int mul = 0;
    for (int i = 0; i <= n - 1; i++)
        mul += (int)pow(2, i) * ncr(n - 1, i);
    int ans = 0;
    for (int i = 0; i < n; i++)
        ans += mul * arr[i];

    return ans;
}

/// Subset Sum


/// digit sum
ll dp[pos_max][total_possible_sum][2];
int solve(int pos,int rem,int isSmall)
{
    if(pos==s.sz)
        return (rem == 0 );
    if(dp[pos][rem][isSmall]!=-1)
        return dp[pos][rem][isSmall];
    int low = 0 , high = (isSmall==false) ? s[pos]-'0' : 9 ;
    int ret = 0;
    for(int i=low;i<=high;i++)
    {
       int new_isSmall = isSmall | (i<high);
       int val = solve(pos+1,(rem+i)%d,new_isSmall);
       ret = modAdd(ret ,val);
    }
    return dp[pos][rem][isSmall] = ret;
}
int main()
{
      cin >> s >> d;
      MEM(dp,-1);
      int ans =  solve(0,0,0) - 1; /// 000 not countable
}





/// digit dp

ll dp[2][11][99][99];
vector<int>inp;
ll func(int pos,int isSmall,int dig_rem,int num_rem)
{
    if(pos==limit)
    {
        if( (dig_rem+num_rem) == 0)
            return 1;
        else return 0;
    }

    ll &ret = dp[isSmall][pos][dig_rem][num_rem];
    if(ret!=-1)
        return ret;

    ll low = 0 , high = inp[pos], re = 0;
    if(isSmall) high = 9;
    for(int i=low;i<=high;i++)
    {
        int newIsSmall = isSmall | (i<high);
        ll val = func(pos+1,newIsSmall,(dig_rem+i)%k,(num_rem*10+i)%k);
        re = re + val;
    }
    return ret = re;
}
