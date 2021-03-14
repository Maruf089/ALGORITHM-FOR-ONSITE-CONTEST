

/// Longest Increasing Subsequence NlogN
multiset < int > s;
multiset < int > :: iterator it;
FOR(i, 1, n)
{
    s.insert(a[i]);
    it = s.upper_bound(a[i]);
    if(it != s.end())
        s.erase(it);
}
cout << s.size() << endl;



/// Maximum Increasing sum O(NlogN)

vector<ll>tree(4*mx,0);
long long Query( int node, int b, int e, int i, int j)
{
    if( b >= i && e <= j)
        return tree[node];
    if( j<b || i>e )
        return 0;

    int Left = node*2;
    int Right = node*2+1;
    int mid  = (b+e)/2;
    long long p1 = Query( Left, b, mid, i,j);
    long long p2 = Query( Right, mid+1, e, i,j);
    return max(p1,p2);
}
void Update( int node, int b, int e, int i, int j,ll newvalue)
{
    if( b >= i && e <= j)
    {
        tree[node] = max(tree[node],newvalue);
        return ;
    }
    if( j<b || i>e )
        return ;

    int Left = node*2;
    int Right = node*2+1;
    int mid  = (b+e)/2;
    Update( Left, b, mid, i,j,newvalue);
    Update( Right, mid+1, e, i,j,newvalue);
    tree[node] = max(tree[Left],tree[Right]);
}
int main()
{
      cin >> n;
      f0(i,n)
          cin >> val[i];
      Max = 0;
      for(int i=0;i<n;i++)
      {
          ll value = Query(1,0,n-1,0,val[i]-1);
          Update(1,0,n-1,val[i]-1,val[i]-1,value+val[i]);
          Max = max(Max,val[i]+value);
      }
      cout << Max << ln;
}




/// maximum sum for increasing segment
pll arr[mx];
ll st[4*mx];
vector<ll>vc;
map<ll,ll>input;

ll query(int si, int ss, int se, int qs, int qe)
{
    if(qe < ss | qs> se)
        return 0;

    if(ss>=qs && se<=qe)
        return st[si];

    int mid = (ss + se)/2;
    ll l = query(2*si, ss, mid, qs, qe);
    ll r = query(2*si+1, mid+1, se, qs, qe);

    return max(l, r);
}

void update(int n, int b, int e, int idx,ll val)
{
    if (b>e || b>idx || e<idx )
        return;
    if (b==e)
    {
        st[n] = max(val,st[n]);
        return;
    }
    update( n*2, b, (b+e)/2, idx, val );
    update( n*2 + 1, (b+e)/2 + 1, e, idx, val );
    st[n] = max( st[n*2], st[n*2 + 1] );
}
int main()
{
    inp(n);
    f0(i,n)
    {
        inp2(a,b);
        arr[i] = {a,b};
        vc.pb(a);vc.pb(b);
    }
    sort(all(vc));
    SORT_UNIQUE(vc);

    f0(i,vc.sz)
    {
        input[vc[i]] = i+1;
    }

    f0(i,n)
    {
        counT = query(1,1,vc.sz,1,input[arr[i].F]-1);
        update(1,1,vc.sz,input[arr[i].S],counT+arr[i].S-arr[i].F+1);
    }
    ans = st[1];

    printf("%lld\n",ans);
}




/// Longest increasing subsequence ( with Query )

int bit[N],n;
void update(int i,int val)
{
	while(i<=n)
	{
		bit[i]=max(bit[i],val);
		i=i+(i&(-i));
	}
}
int query(int i)
{
	int ret=0;
	while(i)
	{
		ret=max(ret,bit[i]);
		i=i-(i&(-i));
	}
	return ret;
}
vector<int> v[N];
int last[N],maxx[N],a[N];
int Out[N];
vector<pair<int,int> > Q[N];
void update_DAG(int cur,int val)
{
	if(val>maxx[cur])
	{
		for(auto x:v[cur])
			update_DAG(x,val+1);
		maxx[cur]=val;
		update(cur,val);
	}
}
int main()
{
	int t;
	sd(t);
	while(t--)
	{
		int q,i,j;
		sd(n);sd(q);
		for(i=0;i<=n;i++)
		{
			maxx[i]=last[i]=bit[i]=0;
			v[i].clear();
			Q[i].clear();
		}
		for(i=1;i<=n;i++)
		{
			sd(a[i]);
			last[a[i]]=i;
			int prev=0;
			for(j=a[i]-1;j>=1;--j)
				if(last[j]>prev)
				{
					prev=last[j];
					v[last[j]].PB(i);
				}
		}
		for(i=0;i<q;i++)
		{
			int x,y;
			sd(x);sd(y);
			Q[x].PB(MP(y,i));
		}
		for(i=n;i>=1;--i)
		{
			update_DAG(i,1);
			for(j=0;j<Q[i].size();j++)
				Out[Q[i][j].S]=query(Q[i][j].F);
		}
		for(i=0;i<q;i++)
			printf("%d\n",Out[i]);
	}
}


