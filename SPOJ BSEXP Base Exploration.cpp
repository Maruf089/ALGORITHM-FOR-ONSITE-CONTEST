#include<bits/stdc++.h>
#define REP(i,n) for (int i = 1; i <= n; i++)
#define mod 1000000007
#define pb push_back
#define ff first
#define ss second
#define ii pair<int,int>
#define vi vector<lli>
#define vii vector<ii>
#define lli long long int
#define endl '\n'

using namespace std;
vi primes;

void init(int maxN)
{
	vi ar(maxN + 1 , 0);

	for(lli i=2;i<=maxN;i++)
	if(ar[i] == 0)
	{
		primes.pb(i);
		for(lli j=2*i;j<=maxN;j+=i)
		ar[j] = 1;
	}
}

int main()
{
	lli t , n;
	init(1000000);
	ios_base::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);

	cin>>t;
	while(t--)
	{
		cin>>n;

		lli res = 2;
		lli cnt = 0;

		for(lli P : primes)
		if(P * P <= n)
		{
			if(n % P == 0)
			{
				int k = 0;
				while(n % P == 0)
				k++ , n /= P;

				if(cnt < k)
				cnt = k , res = P;
			}
		}
		else
		break;

		if(n > 1 && cnt < 1)
		res = n;

		cout<<res<<endl;
	}
}
