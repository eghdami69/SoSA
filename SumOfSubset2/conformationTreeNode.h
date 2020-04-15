#ifndef CONFORMATIONTREENODE_H
#define CONFORMATIONTREENODE_H
#include<string>
#include<vector>
#include<algorithm>
using namespace std;



class ConformationTreeNode{
public:
	int node_number_of_beta_strands;
	int possible_number_of_strands_for_selection;
	bool end_of_processing;
	bool start_processing;
	int *possible_strands;
	string combinations;
	int *permutation;
	bool *interaction_pattern;
	int* remaining_strands;
	double score;
	ConformationTreeNode *parent;
	ConformationTreeNode* child;
	int factorial_number;
	int total_number_of_permutations;

	

	ConformationTreeNode(int ns, int pns, int* possible_strands, ConformationTreeNode *p);
	

	ConformationTreeNode(ConformationTreeNode *p);
	~ConformationTreeNode();

	void compute_score(double **alignment_scores, int all_protein_number_of_beta_strands);
	void next_state(double **alignment_scores, int all_protein_number_of_beta_strands);
	bool pruning_is_needed(double **alignment_scores, int all_protein_number_of_beta_strands);

};

#endif