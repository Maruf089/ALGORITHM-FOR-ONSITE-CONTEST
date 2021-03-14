

ll pows[mx], hashs[mx];
int getHash(int l,int r)
{
    return  ((hashs[r+1]-(hashs[l]*pows[r-l+1])%MOD)+MOD)%MOD;
}
int main()
{
    int p = 31;
    pows[0] = 1;
    for(int i=1; i<=1000005; i++)
    {
        pows[i] = (pows[i-1]*p)%MOD;
    }
    while(cin >> a >> b)
    {
        for(int i=0; i<a.size(); i++)
        {
            hash_a = (hash_a*p+a[i])%MOD;
        }

        for(int i=0; i<b.size(); i++)
        {
            hashs[i+1] = (hashs[i]*p+b[i])%MOD;
        }
    }
}

--------------> another template
ll mod1=1000000007LL;
ll mod2=1000000009LL;
ll p1=43LL;
ll p2=107LL;
ll pw1[N];
ll pw2[N];
void ini()
{
    pw1[0]=p1;
    pw2[0]=p2;
    for(int i=1; i<N-1; i++)
    {
        pw1[i]=(pw1[i-1]*p1)%mod1;
        pw2[i]=(pw2[i-1]*p2)%mod2;
    }
}

ll get_hash1(string& s)
{
    ll val=0;
    for(int i=0; i<s.size(); i++)
        val=(val + ((s[i]-'a'+1)*pw1[i]))%mod1;
    return val;
}
ll get_hash2(string& s)
{
    ll val=0;
    for(int i=0; i<s.size(); i++)
        val=(val+((s[i]-'a'+1)*pw2[i]))%mod2;
    return val;
}
ll get_new_hash(ll H,char prv,int pos,char x) /// delete prv character and add x character in that position
{
    H%=mod1;
    H = (H - ((prv - 'a' + 1)*pw1[pos]))%mod1;
    while(H<0)
        H += mod1;
    H = (H + ((x - 'a' + 1) * pw1[pos]))%mod1;
    while(H<0)
        H += mod1;
    H %= mod1;
    return H;
}

ll get_new_hash2(ll H,char prv,int pos,char x)
{
    H%=mod2;
    H = (H - ((prv - 'a' + 1)*pw2[pos]))%mod2;
    while(H<0)
        H += mod2;
    H = (H + ((x - 'a' + 1) * pw2[pos]))%mod2;
    while(H<0)
        H += mod2;
    H %= mod2;
    return H;
}
int main()
{
    ini();
    cin>>n>>q;
    while(n--)
    {
        cin>>s;
        int len=s.size();
        ll now = get_hash1(s);
        ll now2 = get_hash2(s);
        for(int j = 0; j<s.size(); j++)
        {
            for(char k='a'; k<='c'; k++) /// only a,b,c change
            {
                if(k==s[j])
                    continue;
                ll e = get_new_hash(now,s[j],j,k);
                ll f = get_new_hash2(now2,s[j],j,k);
                mp[ {e,f}] = 1;
            }
        }
    }
    while(q--)
    {
        string str;
        cin>>str;
        ll x = get_hash1(str);
        ll y = get_hash2(str);
        if(mp.count({x,y}))
            cout<<"YES"<<endl;
        else
            cout<<"NO"<<endl;
    }

}





/// get hash in middle parts

-------a------b-------c
here a,b,c is length
A = 0 to a
    B = a+1 to b
        C = b+1 to C

            A = hash_val[id][a] % mod[id];
B=(mod[id]+hash_val[id][a+b]-(hash_val[id][a]*pw[id][b])%mod[id])%mod[id];
C=(mod[id]+hash_val[id][a+b+c]-(hash_val[id][a+b]*pw[id][c])%mod[id])%mod[id];




///
#define LL long long
struct Hashing
{
    LL *hash1, *hash2;
    LL *inv1, *inv2;
    int n;
    LL mod1 = (LL) 1e9 + 97, mod2 = (LL) 1e9 + 9;
    LL multiplier1 = 43, multiplier2 = 31;
    LL invMultiplier1 = 441860508, invMultiplier2 = 838709685;
    // invMultiplier = modInv(multiplier, mod) //

    Hashing()
    {

    }

    Hashing(string &s)
    {
        build_Hash(s);
    }

    void build_Hash(string &s)
    {
        n = s.size();
        hash1 = new LL[n + 1];
        hash2 = new LL[n + 1];
        inv1 = new LL[n + 1];
        inv2 = new LL[n + 1];

        hash1[0] = hash2[0] = 0;
        inv1[0] = inv2[0] = 1;

        LL p1 = 1, p2 = 1;

        for (int i = 0; i < n; i++)
        {
            hash1[i + 1] = (hash1[i] + s[i] * p1) % mod1;
            p1 = (p1 * multiplier1) % mod1;
            inv1[i + 1] = inv1[i] * invMultiplier1 % mod1;
            hash2[i + 1] = (hash2[i] + s[i] * p2) % mod2;
            p2 = (p2 * multiplier2) % mod2;
            inv2[i + 1] = inv2[i] * invMultiplier2 % mod2;
        }
    }


    LL getHash(int i, int j)   //0-based, hash of substring [i, j]
    {
        return getHash_2(i, j - i + 1);
    }

    LL getHash_2(int i, int len)   //0- based, hash of substring [i, i+len-1]
    {
        return (((hash1[i + len] - hash1[i] + mod1) * inv1[i] % mod1) << 32)
               + (hash2[i + len] - hash2[i] + mod2) * inv2[i] % mod2;
    }

    LL revHash(int i, int j)   //0-based
    {
        return getHash(n - j - 1, n - i - 1);
    }

    void clear()
    {
        delete hash1;
        delete hash2;
        delete inv1;
        delete inv2;
    }
};



/// Maximum prefix that is also a suffix

const ll P1=1e14+31,P2=1e14+67;
const int BASE1=777,BASE2=4396;
ll p1[1000001],p2[1000001];

int id(char c)
{
    if(isupper(c))
        return c-'A'+1;
    else if(islower(c))
        return c-'a'+27;
    else
        return c-'0'+53;
}
void init()
{
    p1[0]=1;
    p2[0]=1;
    for(int i=1; i<=1000000; ++i)
    {
        p1[i] = (p1[i-1]*BASE1)%P1;
        p2[i] = (p2[i-1]*BASE2)%P2;
    }
}

ll hs1=0,hs2=0;
ll ht1=0,ht2=0;
int len=0;
for(int i=0; i<s.sz and i<t.sz; ++i)
{
    hs1=(hs1+id(s[s.sz-i-1])*p1[i])%P1;
    hs2=(hs2+id(s[s.sz-i-1])*p2[i])%P2;

    ht1=(ht1*BASE1+id(t[i]))%P1;
    ht2=(ht2*BASE2+id(t[i]))%P2;

    if(hs1==ht1&&hs2==ht2)
    {
        len = max(len,i+1);
    }
}

///////////////////////////////
