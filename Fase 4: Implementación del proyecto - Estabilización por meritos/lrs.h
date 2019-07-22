#ifndef LRS_H
#define LRS_H

#include <cstdio>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#define MAX_CHAR 256

struct CDSmaxLRS{
int first;
int index;
std::string maxsubstring;
};

struct SuffixTreeNode {
    struct SuffixTreeNode *children[MAX_CHAR];

    //pointer to other node via suffix link
    struct SuffixTreeNode *suffixLink;

    /*(start, end) interval specifies the edge, by which the
     node is connected to its parent node. Each edge will
     connect two nodes,  one parent and one child, and
     (start, end) interval of a given edge  will be stored
     in the child node. Lets say there are two nods A and B
     connected by an edge with indices (5, 8) then this
     indices (5, 8) will be stored in node B. */
    int start;
    int *end;

    /*for leaf nodes, it stores the index of suffix for
      the path  from root to leaf*/
    int suffixIndex;
};

struct CDSmaxLRS lrs (std::string);

#endif