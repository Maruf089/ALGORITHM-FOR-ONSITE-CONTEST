/// In-Out DP in Tree
void dfs(int node,int parent) /// compute in[] array
{
    in[node] = 0 ;
    for(auto next : vc[node])
    {
        if(next==parent)
            continue;
        dfs(next,node);
        in[node] = max(in[next]+1,in[node]);
    }
}
void dfs2(int node, int parent) /// compute out array
{
    int Max1 = -1;
    int Max2 = -1;
    for(auto next : vc[node])
    {
        if(next==parent)
            continue;
        if(in[next]>=Max1)
            Max2 = Max1,Max1 = in[next];
        else if(in[next]>Max2)
            Max2 = in[next];
    }
    for(auto next : vc[node])
    {
        if(next==parent)
            continue;
        int use = Max1;
        if(Max1==in[next])
            use = Max2;
        out[next] = max(1+out[node],2+use); /// 1+max(out[p],use+1)
        dfs2(next,node);
    }
}
dfs(1,-1);
dfs2(1,-1);
f1(i,vertices)
    cout << in[i] << ' ' << out[i] << ln ;
////////////////////////////////////////////////////

/// Sum of distances to all nodes from each node
ll predfs(int node,int par)  /// subtree child count
{
    Subtree_Child[node] = 1;
    Subtree_Sum[node] = 0;
    for(auto it : vc[node])
    {
        if(par==it)
            continue;
        predfs(it,node);
        Subtree_Child[node] += Subtree_Child[it];
        Subtree_Sum[node] += (Subtree_Child[it]+Subtree_Sum[it]);
    }
}
ll dfs(int node,int par) /// Sum of distances of all node
{
    for(auto it : vc[node])
    {
        if(it==par)
            continue;
        ll partial_ans = ans[node] -  (Subtree_Child[it] + Subtree_Sum[it]) ;
        ans[it] = Subtree_Sum[it] +  partial_ans + (n-Subtree_Child[it]);
        dfs(it,node);
    }
}
	predfs(1, -1);
	ans[1] = in[1];
	dfs(1,-1);
	for(int i=1;i<=n;i++)
       cout << ans[i] << ' ';
////////////////////////////////////////////////////////
