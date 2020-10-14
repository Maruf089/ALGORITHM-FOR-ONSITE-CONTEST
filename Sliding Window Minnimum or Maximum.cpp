
vector<int> maxSlidingWindow(vector<int>& nums, int k)
{
    int n = nums.size();
    deque<int>dq;
    vector<int>ans;

    for(i=0; i<n; i++)
    {
        while(!dq.empty() and dq.front()<=i-k)
            dq.pop_front();

        while(!dq.empty() and nums[dq.back()]<nums[i])
            dq.pop_back();

        dq.push_back(i);
        if(i>=k-1)
            ans.push_back(nums[dq.front()]);
    }
    return ans;
}


vector<int> minSlidingWindow(vector<int>& nums, int k)
{
    int n = nums.size();
    deque<int>dq;
    vector<int>ans;

    for(int i=0; i<n; i++)
    {
        while(!dq.empty() and dq.front()<=i-k)
            dq.pop_front();

        while(!dq.empty() and nums[dq.back()]>nums[i])
            dq.pop_back();

        dq.push_back(i);
        if(i>=k-1)
            ans.push_back(nums[dq.front()]);
    }
    for(auto it : ans)
        cout << it << ' ';
}

/// Given a string S and a string T, find the minimum window in S
/// which will contain all the characters in T in complexity O(n).
class Solution
{
public:
    string minWindow(string s,string t)
    {
        int s_len = s.size();
        int t_len = t.size();

        vector<int>freq_t(256,0);
        vector<int>freq_s(256,0);
        for(int i=0;i<t_len;i++)
            freq_t[t[i]]++;

        int win_left = 0 , ans = INT_MAX , ans_left,counT = 0 ;
         for(int i=0;i<s_len;i++)
         {
             if(freq_t[s[i]]>0 and freq_t[s[i]]>freq_s[s[i]])
                counT++;
             freq_s[s[i]]++;
             if(counT==t_len)
             {
                 /*
				All characters required are in current window
                Now , try to minimize window by removing extra character at the start of window
                */

                 while(win_left<i and (freq_t[s[win_left]]<freq_s[s[win_left]] or freq_t[s[win_left]]==0))
                 {
                     /// remove Extra frequency or zero frequency characters
                     if(freq_t[s[win_left]]<freq_s[s[win_left]])
                        freq_s[s[win_left]]--;
                     win_left++;
                 }
                 if(ans>(i-win_left+1))
                 {
                     ans = i-win_left+1;
                     ans_left = win_left;
                 }
             }
         }
         if(ans==INT_MAX)
            return "";
         else return s.substr(ans_left,ans);
    }
};

