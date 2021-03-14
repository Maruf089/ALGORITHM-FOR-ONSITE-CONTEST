#include<bits/stdc++.h>
using namespace std;

struct state
{
    int len,link;
    map<char,int>next;
};
const int MAXLEN = 100000;
state st[MAXLEN * 2];
int node, last;

void sa_init()
{
    st[0].len = 0;
    st[0].link = -1;
    node++;
    last = 0;
}

void sa_extend(char c)
{
    int cur = node++; /// creat new node
    st[cur].len = st[last].len + 1;
    int parent = last;
    while (parent != -1 && !st[parent].next.count(c))
    {
        st[parent].next[c] = cur; /// connect edge(parent,new node)
        parent = st[parent].link;
    }
    if (parent == -1)
    {
        st[cur].link = 0; /// connect link to the root
    }
    else
    {
        int probable_fail_link = st[parent].next[c]; /// found probable fail node
        if (st[parent].len + 1 == st[probable_fail_link].len)
        {
            st[cur].link = probable_fail_link; /// coonect edge(new node , probable fail node)
        }
        else
        {
            int clone = node++; /// make clone node
            st[clone].len = st[parent].len + 1;
            st[clone].next = st[probable_fail_link].next;
            st[clone].link = st[probable_fail_link].link;
            while (parent != -1 && st[parent].next[c] == probable_fail_link)
            {
                st[parent].next[c] = clone; /// connect(parent,clone node)
                parent = st[parent].link;
            }
            st[probable_fail_link].link = st[cur].link = clone;
        }
    }
    last = cur;
}
int main()
{

}



/// class template of Suffix Automaton
class suffix_automaton
{
    struct state
    {
        int len,link = 0;
        bool terminal = false;
        bool is_clone = false;
        map<char,int>next;
        state(int len=0) : len(len){}
        bool have_next(char ch)
        {
            return ( next.find(ch) != next.end() );
        }
    };
    vector<state>st;
    vector<long long>dp;
    int last = 0;
public:
    suffix_automaton()
    {
        st.push_back(state());
        st[0].link = -1;
    }
    suffix_automaton(const string &s) : suffix_automaton()
    {
        for(char ch : s)
          sa_extend(ch);
        mark_terminals();
        dp = vector<long long>(st.size());
    }

    void mark_terminals()
    {
        for(int cur = last;cur>=0;cur = st[cur].link)
            st[cur].terminal = true;
    }

    void sa_extend(char c)
    {
        int cur = st.size(); /// creat new node
        st.push_back(state(st[last].len + 1));
        int parent = last;
        last = cur;

        while (parent != -1 && !st[parent].have_next(c))
        {
            st[parent].next[c] = cur; /// connect edge(parent,new node)
            parent = st[parent].link;
        }
        if(parent == -1)
        {
             st[cur].link = 0;
        }
        else
        {
            int probable_fail_link = st[parent].next[c]; /// found probable fail node
            if (st[parent].len + 1 == st[probable_fail_link].len)
            {
                st[cur].link = probable_fail_link; /// coonect edge(new node , probable fail node)
            }
            else
            {
                int clone = st.size(); /// make clone node
                st.push_back((st[probable_fail_link]));
                st[clone].is_clone = true;
                st[clone].len = st[parent].len + 1;
                st[clone].next = st[probable_fail_link].next;
                st[clone].link = st[probable_fail_link].link;
                while (parent != -1 && st[parent].next[c] == probable_fail_link)
                {
                    st[parent].next[c] = clone; /// connect(parent,clone node)
                    parent = st[parent].link;
                }
                st[probable_fail_link].link = st[cur].link = clone;
            }
        }
    }
    long long distict_substring(int vertex)
    {
        long long &ret = dp[vertex];
        if(ret)
            return ret;
        ret = 1;
        for(auto it : st[vertex].next)
            ret += distict_substring(it.second);
        return ret;
    }

};
