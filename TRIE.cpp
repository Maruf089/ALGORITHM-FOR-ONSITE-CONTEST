/// max,min x-or sum for an array LOJ 1269


int tree[mx*32][2];
int cnt[mx*32][2];
int nodeNum;
void Init()
{
    nodeNum = 0 ;
    tree[0][0] = tree[0][1] = 0;
}
void Insrt(ll x)
{
    ll pos=0, bit;
    for(int i=31; i>=0; i--)
    {
        bit = (x>>i)&1;
        if(bit)
        {
            if(tree[pos][1] == 0)
            {
                tree[pos][1] = ++nodeNum;
                MEM(tree[nodeNum],0);
            }
            cnt[pos][1]++;
            pos = tree[pos][1];
        }
        else
        {
            if(tree[pos][0] == 0 )
            {
                tree[pos][0] = ++nodeNum;
                MEM(tree[nodeNum],0);
            }
            cnt[pos][0]++;
            pos = tree[pos][0];
        }
    }
    return;
}

ll Qry_max(ll x)
{
    ll pos=0, bit , res = 0;
    for(int i=31; i>=0; i--)
    {
        bit = (x>>i)&1;
        if(bit)
        {
            if(tree[pos][0]!=0 and cnt[pos][0]>0)
            {
                res |= (1<<i);
                pos = tree[pos][0];
            }
            else
                pos = tree[pos][1];
        }
        else
        {
            if(tree[pos][1]!=0 and cnt[pos][1]>0)
            {
                res |= (1<<i);
                pos = tree[pos][1];
            }
            else
                pos = tree[pos][0];
        }
    }
    return res;
}

int Qry_min(int x)
{
    int pos=0, bit, res=0;
    for(int i=31; i>=0; i--)
    {
        bit = (x>>i)&1;
        if(bit==0)
        {
            if(tree[pos][0]!=0)
                pos = tree[pos][0];
            else
            {
                res |= (1<<i);
                pos=tree[pos][1];
            }
        }
        else
        {
            if(tree[pos][1]!=-1)
                pos=tree[pos][1];
            else
            {
                pos=tree[pos][0];
                res |= (1<<i);
            }
        }
    }
    return res;
}

void Delete(ll a)
{
    ll pos = 0 , bit;
    for(int i=31; i>=0; i--)
    {
         bit = (a>>i)&1;
         if(bit)
         {
             cnt[pos][1]--;
             pos = tree[pos][1];
         }
         else
         {
             cnt[pos][0]--;
             pos = tree[pos][0];
         }
    }
    return;
}


int main()
{

    int tc, kk=1, n, x, cs = 0,cumxor, mx,mn;
    cin>>tc;
    while(tc--)
    {
        cin>>n;
        nodeNum=mx=cumxor=0,mn=INT_MAX;

        Init();
        Insrt(0);
        for(int i=0; i<n; i++)
        {
            cin>>x;
            cumxor=cumxor^x;
            mx=max(mx, Qry_max(cumxor));
            mn=min(mn, Qry_min(cumxor));
            Insrt(cumxor);
        }
        printf("Case %d: %d %d\n",++cs,mx,mn);
    }
}


/// Normal 2d trie tree

const int ALPHABET = 26;
const int MAX = 111111;

int scale(char c)
{
    if(c >= 'a' && c <= 'z')
        return c - 'a';
}
int idx, root;
int tree[MAX][ALPHABET];
int data[MAX]; /// there can be more data fields

void init()
{
    root = 0, idx = 1;
    data[root] = 0;
    memset(tree[root], -1, sizeof(tree[root]));
}
void insert(string s)
{
    int curr = 0, i, k;
    for(i = 0; s[i]; i++)
    {
        k = scale(s[i]);
        if(tree[curr][k] == -1)
        {
            tree[curr][k] = idx;
            data[idx] = 0 ;
            memset(tree[idx], -1, sizeof(tree[idx]));
            idx++;
        }
        curr = tree[curr][k];
    }
    data[curr]++;
}
int search(string s)
{
    int curr = 0, i, k;
    for(i = 0; s[i]; i++)
    {
        k = scale(s[i]);
        if(tree[curr][k] == -1)
            return 0;
        curr = tree[curr][k];
    }
    return data[curr];
}





