/* functions for python style string manipulation */
#if !defined(STRINGTOOLS_HPP)
#define STRINGTOOLS_HPP 1

#include <string> 
#include <vector>

using namespace std;

/* join string_vec into a long string using seperator "sep" */
string Join(const string sep, const vector<string>& string_vec)
{
    if (string_vec.size()==0) return "";
    string joined_str=string_vec[0];
    for (int s=1;s<string_vec.size();s++) joined_str+=sep+string_vec[s];
    return joined_str;
}


/* split a long string into vectors by whitespace 
 * line          - input string
 * line_vec      - output vector 
 * delimiter     - delimiter */
void Split(const string &line, vector<string> &line_vec,
    const char delimiter=' ')
{
    bool within_word = false;
    for (size_t pos=0;pos<line.size();pos++)
    {
        if (line[pos]==delimiter)
        {
            within_word = false;
            continue;
        }
        if (!within_word)
        {
            within_word = true;
            line_vec.push_back("");
        }
        line_vec.back()+=line[pos];
    }
}

/* strip white space at the begining or end of string */
string Trim(const string &inputString,const string &char_list=" \n\r\t")
{
    string result = inputString;
    int idxBegin = inputString.find_first_not_of(char_list);
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxBegin >= 0 && idxEnd >= 0)
        result = inputString.substr(idxBegin, idxEnd + 1 - idxBegin);
    return result;
}

string rstrip(const string &inputString,const string &char_list=" \n\r\t")
{
    string result=inputString;
    int idxEnd = inputString.find_last_not_of(char_list);
    if (idxEnd >= 0)
        result = inputString.substr(0, idxEnd + 1);
    return result;
}

inline bool StartsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(0,shortString.size())==shortString);
}

inline bool EndsWith(const string &longString, const string &shortString)
{
    return (longString.size()>=shortString.size() &&
            longString.substr(longString.size()-shortString.size(),
                shortString.size())==shortString);
}
#endif
