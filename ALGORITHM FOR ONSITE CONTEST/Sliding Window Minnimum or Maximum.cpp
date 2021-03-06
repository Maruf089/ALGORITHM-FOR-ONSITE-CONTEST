
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
