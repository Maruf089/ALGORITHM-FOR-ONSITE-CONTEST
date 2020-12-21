/// string with space convert with string stream

char s[10];
scanf("%[^\n]",s);

char s1[201], ss[100005] ;
getchar();   /// need this for input
ans = 1;
fgets(ss, sizeof(ss), stdin);
string new_str(ss,ss+strlen(ss));
stringstream get_ss;
get_ss << new_str;
while(get_ss>>new_str)
{
    if(new_str.sz>3)
        sort(new_str.begin()+1,new_str.end()-1);
    ans *= mp[new_str];
}
