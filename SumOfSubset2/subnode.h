#ifndef SUBNODE_H
#define SUBNODE_H

#include<vector>
#include "treeNode.h"
using namespace std;
class treeNode;
class subnode{
public:
	int data;
	vector <subsetNode> children;

	subnode(int d){
		data = d;
	}
};

#endif