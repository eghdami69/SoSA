#ifndef PROTEIN_H
#define PROTEIN_H

#include<vector>
#include"sheet.h"
#include"conformationTreeNode.h"
using namespace std;

#define PARALLEL 1
#define ANTIPARALLEL 0
class Protein
{
private:

	//void set_setcontactmap(int* max_states_of_two_beta_strands, int m, int n, int s1, int s2, int type);
	void set_setcontactmap(int s1, int s2, int * max_states, int* beta_residue_positions, int* segment_starts, int* segment_ends, int number_of_beta_residues, int parallel);
public:
	int number_of_sheets;
	vector<Sheet *> sheets;
	int **contact_map;
	double protein_score;
	int number_of_beta_strands;
	int number_of_beta_residues;
	//Prob prob;
	Protein(int number_of_strands, int number_of_beta_residues);
	~Protein();
	void init(ConformationTreeNode* c, double score);
	void make_contact_map(int *start_points, int *** max_states, int *segment_starts, int * segment_ends);
	void print();

};

#endif