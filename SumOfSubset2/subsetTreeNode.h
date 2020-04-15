#ifndef TREENODE_H
#define TREENODE_H

#include<vector>
using namespace std;
class SubsetTreeNode{
public:
	int data;
	SubsetTreeNode *parent;
	vector<SubsetTreeNode*> children;
	int leaf_number;

	SubsetTreeNode(int d, SubsetTreeNode* p){
		data = d;
		parent = p;
		leaf_number = -1;
	}

	SubsetTreeNode(int d, SubsetTreeNode* p, int l_num){
		data = d;
		parent = p;
		leaf_number = l_num;
	}

	~SubsetTreeNode(){
		for (int i = 0; i < children.size(); i++){
			delete(children[i]);// ->~SubsetTreeNode();
		}
	}


};


#endif