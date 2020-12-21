#include<bits/stdc++.h>
using namespace std;

vector<int>construct_lps(string pattern)
{
    vector<int>lps(pattern.size());
    int index = 0;
    for(int i=1;i<pattern.size();)
    {
        if(pattern[i]==pattern[index])
        {
            lps[i] = index + 1;
            index++ , i++;
        }
        else
        {
            if(index!=0) index = lps[index-1];
            else lps[i] = index , i++;
        }
    }
    return lps;
}

void KMP(string text,string pattern)
{
    bool found = false;
    vector<int>lps = construct_lps(pattern);
    int j = 0 , i = 0 ;
      /// i --> text, j --> pattern
    while(i<text.size())
    {
        if(text[i]==pattern[j])
            i++ , j++;
        else
        {
            if(j!=0) j = lps[j-1];
            else i++;
        }

        if(j==pattern.size())
        {
            found = true;
            cout << "Found at index : " << (i-pattern.size()) << endl;
            j = lps[j-1];
        }
    }
    if(!found)
        cout << "Not Found \n" ;
}


int main()
{
    string text,pattern;
    cin >> text >> pattern;

    KMP(text,pattern);
    return 0;
}
