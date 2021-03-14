
/// Replace (i to j) -> into Val
/// print average (sum/(j-i+1) from (i to j)

ll sum[4*mx],Lazy[4*mx];
void LazYPropagrate(int idx,int left,int right)
{
    if(left!=right)
    {
        Lazy[2*idx] = Lazy[idx];
        Lazy[2*idx+1] = Lazy[idx];
    }
    sum[idx] = Lazy[idx]*((right-left)+1);
    Lazy[idx] = -1;
}

void update(int node,int left,int right,int x,int y,ll val)
{
    if(Lazy[node]!=-1)
        LazYPropagrate(node,left,right);
    if(left>y || right < x)
        return ;
    if( x<=left && right<=y )
    {

        sum[node] = val*((right-left)+1);
        if(left!=right)
        {
            Lazy[2*node] = val;
            Lazy[2*node+1] = val;
        }
        return ;
    }
    int mid =(right + left ) /2;
    update(2*node,left,mid,x,y,val);
    update(2*node+1,mid+1,right,x,y,val);
    sum[node] = sum[2*node]+ sum[2*node+1];
}
ll query(int node,int left,int right,int x,int y)
{
    if(Lazy[node]!=-1)
        LazYPropagrate(node,left,right);
    if(left>y || x>right || left>right)
        return 0;
    if(x<=left && y>=right )
        return sum[node];
    int mid = (left + right ) /2 ;
    ll res = 0;
    res = query(2*node,left,mid,x,y) + query(2*node+1,mid+1,right,x,y);
    return res;
}


int main()
{
    inp(t);
    while(t--)
    {
        MEM(Lazy,-1);
        MEM(sum,0);

        inp2(n,m);

        printf("Case %lld:\n",++cs);
        while(m--)
        {

            inp(k);
            if(k==1)
            {
                inp2(a,b);
                inp(l);
                update(1,0,n-1,a,b,l);
            }
            else
            {
                inp2(a,b);
                ans = query(1,0,n-1,a,b);
                if(ans%(b-a+1)==0)
                    printf("%lld\n",ans/(b-a+1));
                else
                {
                    ll gcd = __gcd(ans,b-a+1);
                    ans /= gcd;
                    num = b-a+1;
                    num /= gcd;
                    printf("%lld/%lld\n",ans,num);
                }
            }

        }


    }
}


