#include<iostream>
#include"conformationTreeNode.h"
using namespace std;


int fact(int n){
	int res = 1;
	for (int i = 1; i <= n; i++){
		res *= i;
	}
	return res;
}

bool increment(bool *ip, int len){
	int i;
	for (i = 0; i < len; i++){
		if (!ip[i]){
			ip[i] = true;
			return false;
		}
		ip[i] = false;
	}
	if (i == len){
		return true;
	}
	else{
		return false;
	}
}

ConformationTreeNode:: ConformationTreeNode(int ns, int pns, int* possible_strands, ConformationTreeNode *p){
	parent = p;
	node_number_of_beta_strands = ns;
	possible_number_of_strands_for_selection = pns;
	factorial_number = 0;
	total_number_of_permutations = fact(node_number_of_beta_strands) / 2;
	end_of_processing = false;
	start_processing = true;
	this->possible_strands = possible_strands;
	combinations.resize(node_number_of_beta_strands, 1); // K leading 1's
	combinations.resize(possible_number_of_strands_for_selection, 0);// N-K trailing 0's
	remaining_strands = new int[possible_number_of_strands_for_selection - node_number_of_beta_strands];
	permutation = new int[node_number_of_beta_strands];
	interaction_pattern = new bool[node_number_of_beta_strands - 1];
	for (int i = 0; i < node_number_of_beta_strands - 1; i++){//is necessary???????
		interaction_pattern[i] = 0;
	}
	int k = 0;
	int j = 0;
	for (int i = 0; i < possible_number_of_strands_for_selection; i++){
		if (combinations[i] == 0){
			remaining_strands[k++] = possible_strands[i];
		}
		else{
			permutation[j++] = possible_strands[i];
		}
	}
}

ConformationTreeNode:: ConformationTreeNode(ConformationTreeNode *p){
	parent = p;
}
ConformationTreeNode:: ~ConformationTreeNode(){
	delete(remaining_strands);
	delete(interaction_pattern);
	delete(permutation);
}

void ConformationTreeNode::compute_score(double **alignment_scores, int all_protein_number_of_beta_strands){
	score = 0;
	for (int i = 0; i < node_number_of_beta_strands - 1; i++){
		/*if (i == 0){
			cout << permutation[i]+1;
		}*/
		if (interaction_pattern[i]){
			(score) += alignment_scores[permutation[i]][permutation[i + 1]];
			//cout << 'p' << permutation[i + 1]+1;
		}
		else{
			(score) += alignment_scores[permutation[i]][permutation[i + 1] + all_protein_number_of_beta_strands];
			//cout << "ap" << permutation[i + 1]+1;
		}
	}
	//cout << endl;
	//score /= (node_number_of_beta_strands - 1);
	
}

bool ConformationTreeNode::pruning_is_needed(double **alignment_scores, int all_protein_number_of_beta_strands){
	return false;
	double threshold = 0.01;
	for (int i = 0; i < node_number_of_beta_strands - 1; i++){
		if (interaction_pattern[i]){
			if (alignment_scores[permutation[i]][permutation[i + 1]] < threshold){
				//cout << "**" << endl;
				return true;
			}
			if (i + 1 < node_number_of_beta_strands - 1){
				if (interaction_pattern[i + 1]){
					if (alignment_scores[permutation[i]][permutation[i + 1]] + alignment_scores[permutation[i + 1]][permutation[i + 2]] < alignment_scores[permutation[i]][permutation[i + 2]]
						|| alignment_scores[permutation[i]][permutation[i + 1]] + alignment_scores[permutation[i + 1]][permutation[i + 2]] < alignment_scores[permutation[i]][permutation[i + 2]+all_protein_number_of_beta_strands]){
						//cout << "**" << endl;
						return true;
					}
				}
				else{
					if (alignment_scores[permutation[i]][permutation[i + 1]] + alignment_scores[permutation[i + 1]][permutation[i + 2]+node_number_of_beta_strands] < alignment_scores[permutation[i]][permutation[i + 2]]
						|| alignment_scores[permutation[i]][permutation[i + 1]] + alignment_scores[permutation[i + 1]][permutation[i + 2]+node_number_of_beta_strands] < alignment_scores[permutation[i]][permutation[i + 2] + all_protein_number_of_beta_strands]){
						//cout << "**" << endl;
						return true;
					}
				}
			}
		}
		else{
			if (alignment_scores[permutation[i]][permutation[i + 1]+all_protein_number_of_beta_strands] < threshold){
				//cout << "**" << endl;
				return true;
			}
			if (i + 1 < node_number_of_beta_strands - 1){
				if (interaction_pattern[i + 1]){
					if (alignment_scores[permutation[i]][permutation[i + 1]+node_number_of_beta_strands] + alignment_scores[permutation[i + 1]][permutation[i + 2]] < alignment_scores[permutation[i]][permutation[i + 2]]
						|| alignment_scores[permutation[i]][permutation[i + 1]+node_number_of_beta_strands] + alignment_scores[permutation[i + 1]][permutation[i + 2]] < alignment_scores[permutation[i]][permutation[i + 2] + all_protein_number_of_beta_strands]){
						//cout << "**" << endl;
						return true;
					}
				}
				else{
					if (alignment_scores[permutation[i]][permutation[i + 1]+node_number_of_beta_strands] + alignment_scores[permutation[i + 1]][permutation[i + 2] + node_number_of_beta_strands] < alignment_scores[permutation[i]][permutation[i + 2]]
						|| alignment_scores[permutation[i]][permutation[i + 1]+node_number_of_beta_strands] + alignment_scores[permutation[i + 1]][permutation[i + 2] + node_number_of_beta_strands] < alignment_scores[permutation[i]][permutation[i + 2] + all_protein_number_of_beta_strands]){
						//cout << "**" << endl;
						return true;
					}
				}
			}
		}
	}
	return false;

}

void ConformationTreeNode::next_state(double **alignment_scores, int all_protein_number_of_beta_strands){
	if (!start_processing){
		do{
			if (increment(interaction_pattern, node_number_of_beta_strands - 1)){
				do{
					std::next_permutation(permutation, permutation + node_number_of_beta_strands);
				} while (permutation[0] > permutation[node_number_of_beta_strands - 1]);
				factorial_number++;
				if (factorial_number >= total_number_of_permutations){
					factorial_number = 0;
					if (!std::prev_permutation(combinations.begin(), combinations.end())){
						end_of_processing = true;
						break;///////////////////////
					}
					else{
						int k = 0;
						int j = 0;
						for (int i = 0; i < possible_number_of_strands_for_selection; i++){
							if (combinations[i] == 0){
								remaining_strands[k++] = possible_strands[i];
							}
							else{
								permutation[j++] = possible_strands[i];
							}
						}
					}
				}

			}
		} while (pruning_is_needed(alignment_scores, all_protein_number_of_beta_strands));
	}
	else{
		start_processing = false;
		if (pruning_is_needed(alignment_scores, all_protein_number_of_beta_strands)){
			do{
				if (increment(interaction_pattern, node_number_of_beta_strands - 1)){
					do{
						std::next_permutation(permutation, permutation + node_number_of_beta_strands);
					} while (permutation[0] > permutation[node_number_of_beta_strands - 1]);
					factorial_number++;
					if (factorial_number >= total_number_of_permutations){
						factorial_number = 0;
						if (!std::prev_permutation(combinations.begin(), combinations.end())){
							end_of_processing = true;
							break;
						}
						else{
							int k = 0;
							int j = 0;
							for (int i = 0; i < possible_number_of_strands_for_selection; i++){
								if (combinations[i] == 0){
									remaining_strands[k++] = possible_strands[i];
								}
								else{
									permutation[j++] = possible_strands[i];
								}
							}
						}
					}

				}
			} while (pruning_is_needed(alignment_scores, all_protein_number_of_beta_strands));
		}
	}
}