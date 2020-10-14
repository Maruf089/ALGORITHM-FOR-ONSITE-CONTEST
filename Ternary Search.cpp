
ll ternary2(ll l,ll r)
{
    while(l<r-2)
    {
        ll mid1 = l+(r-l)/3;
        ll mid2 = mid1+(r-l)/3;
        if(num2(mid1)<=num2(mid2))
            r = mid2;
        else l = mid1;
    }
    double val = 1e9;
    ll idx = l;
    for(i=l;i<=r;i++)
    {
        if(val>num2(i))
        {
            idx = i;
            val = num2(i);
        }
    }
    return idx;
}



  while(abs(ac-bd)>EPS)
    {
        pair<double,double> a,b,c,d;

        a.first=(2*ar[0].first+ar[1].first)/3.0;  /// jei pase nibo oi pase 2* (2*left+right)/3 ans (2*right+left)/3
        a.second=(2*ar[0].second+ar[1].second)/3.0;
        b.first=(ar[0].first+2*ar[1].first)/3.0;
        b.second=(ar[0].second+2*ar[1].second)/3.0;

        c.first=(2*ar[2].first+ar[3].first)/3.0;
        c.second=(2*ar[2].second+ar[3].second)/3.0;
        d.first=(ar[2].first+2*ar[3].first)/3.0;
        d.second=(ar[2].second+2*ar[3].second)/3.0;

        ac=dist(a,c);
        bd=dist(b,d);

        if(ac>bd)ar[0]=a,ar[2]=c;
        else ar[1]=b,ar[3]=d;
    }
