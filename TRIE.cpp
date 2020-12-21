/// max,min x-or sum for an array LOJ 1269
struct node
{
    int lft, rgt;
    node()
    {lft=rgt=-1;}
};
node tree[3500000];
int nodeNum;

void Insrt(int x)
{
    int pos=0, bit;
    for(int i=31; i>=0; i--)
    {
        bit = (x>>i)&1;
        if(bit)
        {
            if(tree[pos].rgt==-1)
            {
                tree[pos].rgt=++nodeNum;
                tree[nodeNum]=node();
            }
            pos=tree[pos].rgt;
        }
        else
        {
            if(tree[pos].lft==-1)
            {
                tree[pos].lft=++nodeNum;
                tree[nodeNum]=node();
            }
            pos=tree[pos].lft;
        }
    }
    return;
}

int Qry_max(int x)
{
    int pos=0, bit, res=0;
    for(int i=31; i>=0; i--)
    {
        bit = (x>>i)&1;
        if(bit)
        {
            if(tree[pos].lft!=-1)
            {
                res |= (1<<i);
                pos = tree[pos].lft;
            }
            else
                pos=tree[pos].rgt;
        }
        else
        {
            if(tree[pos].rgt!=-1)
            {
                res |= (1<<i);
                pos=tree[pos].rgt;
            }
            else
                pos=tree[pos].lft;
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
            if(tree[pos].lft!=-1)
                pos = tree[pos].lft;
            else
            {
                res |= (1<<i);
                pos=tree[pos].rgt;
            }
        }
        else
        {
            if(tree[pos].rgt!=-1)
                pos=tree[pos].rgt;
            else
            {
                pos=tree[pos].lft;
                res |= (1<<i);
            }
        }
    }
    return res;
}

int main()
{

    int tc, kk=1, n, x, cs = 0,cumxor, mx,mn;
    cin>>tc;
    while(tc--)
    {
        cin>>n;
        nodeNum=mx=cumxor=0,mn=INT_MAX;

        tree[0] = node();
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



/// using pointer
struct node
{
    node* next[2];
    node()
    {
        next[0]=next[1]=NULL;
    }
}*root;
void insert(int x)
{
    node* cur=root;
    int b;
    for(int i=20; i>=0; i--)
    {
        b = (x>>i)&1;
        if(cur->next[b]==NULL)
            cur->next[b]=new node();
        cur=cur->next[b];
    }
}
int get_min(int x)
{
    node* cur=root;
    int k,ans=0;
    for(int i=20; i>=0; i--)
    {
        k = (x>>i)&1;
        if(cur->next[k])
            cur=cur->next[k],ans<<=1;
        else
            cur=cur->next[!k],ans<<=1,ans++;
    }
    return ans;
}
int main()
{
    root = new node();
}
