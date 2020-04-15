#include"protein.h"
#include "sheet.h"
#include"mem.h"
#include<string>
#include<stdio.h>
#include<iostream>
using namespace std;


#define max_number_of_strands_per_sheet 6
#define max_number_of_sheets 3

Protein::Protein(int number_of_strands, int number_of_beta_residues){
	this->number_of_beta_strands = number_of_strands;
	this->number_of_beta_residues = number_of_beta_residues;
	number_of_sheets = 0;
	protein_score = -1;
	allocate_2d_int(contact_map, number_of_beta_residues, number_of_beta_residues);

}

Protein::~Protein(){
	deallocate_2d_int(contact_map, number_of_beta_residues);
	for (int i = 0; i < sheets.size(); i++){
		delete(sheets[i]);
	}
}

void Protein::init(ConformationTreeNode *c,double score){
	if (number_of_sheets != 0){
		for (int i = 0; i < sheets.size(); i++){
			delete(sheets[i]);
		}
		sheets.clear();
		number_of_sheets = 0;
	}
	protein_score = score;
	ConformationTreeNode *curr=c;
	while (curr != NULL){
		sheets.push_back(new Sheet(curr->node_number_of_beta_strands, curr->permutation, curr->interaction_pattern));
		number_of_sheets++;
		curr = curr->child;
	}
}

void Protein::print(){
	allocate_2d_int(contact_map, number_of_beta_residues, number_of_beta_residues);
	for (int i = 0; i < number_of_sheets; i++){
		for (int j = 0; j < sheets[i]->number_of_strands - 1; j++){
			cout << sheets[i]->spatial_ordering[j];
			if (sheets[i]->interaction_pattern[j]){
				cout << "p";
			}
			else{
				cout << "ap";
			}
		}
		cout << sheets[i]->spatial_ordering[sheets[i]->number_of_strands - 1] << endl;
	}
	cout << "--------------------" << endl;
}



void Protein::make_contact_map(int *beta_residue_positions, int *** beta_strands_max_states, int *segment_starts, int * segment_ends){
	//allocate_2d_int(contact_map, number_of_beta_residues, number_of_beta_residues);
	for (int i = 0; i<number_of_beta_residues; i++){
		for (int j = 0; j<number_of_beta_residues; j++){
			contact_map[i][j] = 0;
		}
	}
	int strand1, strand2, temp;
	for (int i = 0; i<number_of_sheets; i++){
		for (int j = 0; j<sheets[i]->number_of_strands - 1; j++){
			strand1 = sheets[i]->spatial_ordering[j];/////////////////////////////////////////////
			strand2 = sheets[i]->spatial_ordering[j+1];///////////////////////////////////////////
			if (sheets[i]->interaction_pattern[j] == PARALLEL){
				set_setcontactmap(strand1, strand2, beta_strands_max_states[strand1][strand2], beta_residue_positions, segment_starts, segment_ends, number_of_beta_residues, PARALLEL);
			}
			else if (sheets[i]->interaction_pattern[j] == ANTIPARALLEL){
				set_setcontactmap(strand1, strand2, beta_strands_max_states[strand1][number_of_beta_strands+strand2], beta_residue_positions, segment_starts, segment_ends, number_of_beta_residues, ANTIPARALLEL);
			}
		}
	}
}

int get_matrix_indice2(int * beta_residue_positions, int number_of_beta_residues, int i){

	int j;

	//i is the position of the beta-residue

	for (j = 0; j < number_of_beta_residues; j++){

		if (beta_residue_positions[j] == i)
			return j;

	}

	cout << "Error. Incorrect beta-residue position" << endl << flush;

	return 0;



}

void Protein::set_setcontactmap(int s1, int s2, int *max_states, int * beta_residue_positions, int* segment_starts, int* segment_ends, int number_of_beta_residues, int parallel){

	int done;
	int nr, nc;
	int i, j, k, m, n;
	
	/*n = strand_x_length;
	m = strand_y_length;

	nr = m + 1;
	nc = n + 1;
*/
	
	int strand_x_start = segment_starts[s1];
	int strand_y_start = segment_starts[s2];

	int strand_x_end = segment_ends[s1];
	int strand_y_end = segment_ends[s2];

	int strand_x_start_ind = get_matrix_indice2(beta_residue_positions, number_of_beta_residues, strand_x_start);
	int strand_y_start_ind = get_matrix_indice2(beta_residue_positions, number_of_beta_residues, strand_y_start);

	int strand_x_end_ind = get_matrix_indice2(beta_residue_positions, number_of_beta_residues, strand_x_end);
	int strand_y_end_ind = get_matrix_indice2(beta_residue_positions, number_of_beta_residues, strand_y_end);

	int strand_x_length = strand_x_end - strand_x_start + 1;
	int strand_y_length = strand_y_end - strand_y_start + 1;
	
	i = 0; j = 0;
	for (k = 0; ; k++){
		if (i >= strand_x_length || j >= strand_y_length){
			break;
		}
		if (max_states[k] == 0){

			i++; j++;

			if (parallel){

				contact_map[strand_x_start_ind + i - 1][strand_y_start_ind + j - 1] = 1;
				contact_map[strand_y_start_ind + j - 1][strand_x_start_ind + i - 1] = 1;
			}

			else{

				contact_map[strand_x_start_ind + i - 1][strand_y_end_ind - (j - 1)] = 1;
				contact_map[strand_y_end_ind - (j - 1)][strand_x_start_ind + i - 1] = 1;
			}
		}

		else if (max_states[k] == 1){

			i++;
		}

		else if (max_states[k] == 2){

			j++;
		}
	}

}
//void Protein::set_setcontactmap(int* max_states_of_two_beta_strands, int m, int n, int s1, int s2, int type){
//	int done = 0;
//	int k = 0;
//	int index1, index2;
//	int factor;
//	int i, j;
//	if (type == PARALLEL){
//
//		index1 = s1;
//		index2 = s2;
//		i = j = 0;
//		while (!done){
//
//			if (max_states_of_two_beta_strands[k] == 0){ //M state
//				contact_map[index1 + i][index2 + j] = 1;
//				contact_map[index2 + j][index1 + i] = 1;
//				i++;
//				j++;
//
//				//index1--;
//				//index2=index2+factor;
//
//			}
//
//			else if (max_states_of_two_beta_strands[k] == 1){ //Ix state
//
//				j++;
//				//index1--;
//
//			}
//
//			else if (max_states_of_two_beta_strands[k] == 2){ //Iy state
//
//				i++;
//				//index2=index2+factor;
//
//			}
//
//			if ((i >= m) || (j >= n))
//				done = 1;
//
//			k++;
//
//		} //end of while (backtrack)
//	}
//	else{
//		index1 = s1;
//		index2 = s2 + n - 1;
//		//i = m - 1;
//		//j = 0;
//		while (!done){
//
//			if (max_states_of_two_beta_strands[k] == 0){ //M state
//				contact_map[index1][index2] = 1;
//				contact_map[index2][index1] = 1;
//				index1++;
//				index2--;
//
//				//index1--;
//				//index2=index2+factor;
//
//			}
//
//			else if (max_states_of_two_beta_strands[k] == 1){ //Ix state
//
//				index1++;
//				//index1--;
//
//			}
//
//			else if (max_states_of_two_beta_strands[k] == 2){ //Iy state
//
//				index2--;
//				//index2=index2+factor;
//
//			}
//
//			if ((index1 >= s1 + m) || (index2<s2))
//				done = 1;
//
//			k++;
//
//		} //end of while (backtrack)
//
//	}
//
//}
