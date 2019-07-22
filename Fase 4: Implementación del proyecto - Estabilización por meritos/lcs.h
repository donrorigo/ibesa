#ifndef LCS_H
#define LCS_H
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>

struct CDSmaxLCS{
int first;
int second;
int index;
std::string maxsubstring;
};

struct LCP {
    int sa_pos;
    int lcp;
};

void suffixArray(int* s, int* SA, int n, int K);

struct CDSmaxLCS lcs(std::string str1, std::string str2);

#endif

