#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include<stdlib.h>
#include<time.h>
#include<conio.h>
#include<Windows.h>
#include<vector>
#include<iostream>
#include<algorithm>
#include<fstream>
#include<string>
#include"subsetTreeNode.h"
#include"IO.h"
#include"mem.h"
#include"conformationTreeNode.h"
#include"protein.h"
using namespace std;


void discard_filename_extension(char * filename_input){

	int filename_length;
	int i;

	filename_length = (int)strlen(filename_input);

	for (i = 0; i < filename_length; i++){

		if (filename_input[i] == '.')
			break;
	}

	filename_input[i] = '\0';
}


void beta_strand_segment_borders(char * SS, int * segment_starts, int * segment_ends,
	int * number_of_beta_strands, int * number_of_beta_residues)
{

	int j = 0, n;

	n = (int)strlen(SS);

	*number_of_beta_strands = 0;
	*number_of_beta_residues = 0;

	// infinite loop!
	while (1) {

		while ((SS[j] != 'E') && (j < n)) {

			j++;
		}

		if (j >= n) {

			break;
		}

		segment_starts[*number_of_beta_strands] = j;

		while ((SS[j] == 'E') && (j < n)) {

			(*number_of_beta_residues)++;
			j++;
		}

		segment_ends[*number_of_beta_strands] = j - 1;

		if (j >= n) {

			break;
		}

		(*number_of_beta_strands)++;
	}

}

void C2L(char * s){

	int j, length;

	length = (int)strlen(s);

	for (j = 0; j < length; j++){

		if (s[j] == 'C')
			s[j] = 'L';
	}

}


void set_2d_int_to_a_value(int ** array1, int nr, int nc, int value){

	int j, k;

	for (j = 0; j < nr; j++){

		for (k = 0; k < nc; k++)

			array1[j][k] = value;
	}
}

void reverse_int_array(int * array1, int n){

	int max_it = n / 2;
	int k;
	int temp;

	for (k = 0; k < max_it; k++){

		temp = array1[n - k - 1];
		array1[n - k - 1] = array1[k];
		array1[k] = temp;
	}

}

int get_matrix_indice(int * beta_residue_positions, int number_of_beta_residues, int i){

	int j;

	//i is the position of the beta-residue

	for (j = 0; j < number_of_beta_residues; j++){

		if (beta_residue_positions[j] == i)
			return j;

	}

	cout << "Error. Incorrect beta-residue position" << endl << flush;

	return 0;



}
double compute_pairwise_alignment_of_two_beta_strands2(int strand_x_start_ind, int strand_y_start_ind, int strand_x_end_ind, int strand_y_end_ind, int strand_x_length, int strand_y_length, int parallel, double ** beta_residue_pairing_probs, int number_of_beta_residues, int & alignment_length, int * max_states){

	int done;
	int nr, nc;
	int i, j, k, m, n;
	double ** max_alignment_scores;
	int ** max_states_2d; //indices of the maximum value coming from previous states
	//int * max_states;


	double score_i, score_j, score_ij;
	double max_score;
	double max_alignment_score;

	n = strand_x_length;
	m = strand_y_length;

	nr = m + 1;
	nc = n + 1;

	allocate_2d_double(max_alignment_scores, nr, nc);
	allocate_2d_int(max_states_2d, nr, nc);

	//max_states = new int[100];

	max_alignment_scores[0][0] = 0.0;

	for (i = 1; i <= n; i++)
		max_states_2d[0][i] = 1;

	for (j = 1; j <= m; j++)
		max_states_2d[j][0] = 2;

	for (i = 1; i <= n; i++){

		for (j = 1; j <= m; j++){

			if (parallel)
				score_ij = max_alignment_scores[j - 1][i - 1] + beta_residue_pairing_probs[strand_y_start_ind + j - 1][strand_x_start_ind + i - 1];

			else score_ij = max_alignment_scores[j - 1][i - 1] + beta_residue_pairing_probs[strand_y_end_ind - (j - 1)][strand_x_start_ind + i - 1];


			max_states_2d[j][i] = 0;
			max_score = score_ij;

			score_i = max_alignment_scores[j][i - 1] + 0.0; //gap score is chosen as 0.0

			if (score_i > max_score){

				max_score = score_i;
				max_states_2d[j][i] = 1;
			}

			score_j = max_alignment_scores[j - 1][i] + 0.0; //gap score is chosen as 0.0

			if (score_j > max_score){

				max_score = score_j;
				max_states_2d[j][i] = 2;
			}

			max_alignment_scores[j][i] = max_score;

		}
	}

	//BACK TRACE THE STATES ON THE MAXIMUM SCORING PATH 
	max_alignment_score = max_alignment_scores[m][n];
	max_states[0] = max_states_2d[m][n];

	i = n;
	j = m;

	done = 0;
	alignment_length = 0;
	k = 1;

	while (!done){

		if (max_states[k - 1] == 0){ //M state

			i--; j--;
			max_states[k] = max_states_2d[j][i];
			(alignment_length)++;

		}

		else if (max_states[k - 1] == 1){ //Ix state

			i--;
			max_states[k] = max_states_2d[j][i];
			(alignment_length)++;

		}

		else if (max_states[k - 1] == 2){ //Iy state

			j--;
			max_states[k] = max_states_2d[j][i];
			(alignment_length)++;

		}

		if ((i == 0) && (j == 0))
			done = 1;

		k++;

	} //end of while (back trace)

	reverse_int_array(max_states, alignment_length);

	deallocate_2d_double(max_alignment_scores, nr);
	deallocate_2d_int(max_states_2d, nr);
	//delete [] max_states;

	return max_alignment_score;
}



void compute_pairwise_alignments_of_beta_strands2(int * segment_starts, int * segment_ends, int number_of_beta_strands, double ** beta_residue_pairing_probs, int * beta_residue_positions, int number_of_beta_residues, double** max_scores, int*** beta_strand_interactions_max_states){

	int i, j, k;
	int strand_x, strand_y;
	int strand_x_start, strand_y_start, strand_x_end, strand_y_end;
	int strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind;
	int strand_x_length, strand_y_length;
	int parallel;
	int alignment_length;

	double alignment_score_p, alignment_score_ap;

	int * max_states;

	for (strand_x = 0; strand_x < number_of_beta_strands; strand_x++){
		for (strand_y = 0; strand_y < number_of_beta_strands; strand_y++){


			strand_x_start = segment_starts[strand_x];
			strand_y_start = segment_starts[strand_y];

			strand_x_end = segment_ends[strand_x];
			strand_y_end = segment_ends[strand_y];

			strand_x_start_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_x_start);
			strand_y_start_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_y_start);

			strand_x_end_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_x_end);
			strand_y_end_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_y_end);

			strand_x_length = strand_x_end - strand_x_start + 1;
			strand_y_length = strand_y_end - strand_y_start + 1;

			max_states = new int[strand_x_length + strand_y_length + 1];

			//compute anti-parallel alignment
			parallel = 0;
			alignment_score_ap = compute_pairwise_alignment_of_two_beta_strands2(strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind, strand_x_length, strand_y_length, parallel, beta_residue_pairing_probs, number_of_beta_residues, alignment_length, max_states);
			max_scores[strand_x][strand_y + number_of_beta_strands] = alignment_score_ap;
			for (int k = 0; k < alignment_length; k++){
				beta_strand_interactions_max_states[strand_x][strand_y + number_of_beta_strands][k] = max_states[k];
			}
			//compute parallel alignment
			parallel = 1;
			alignment_score_p = compute_pairwise_alignment_of_two_beta_strands2(strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind, strand_x_length, strand_y_length, parallel, beta_residue_pairing_probs, number_of_beta_residues, alignment_length, max_states);
			max_scores[strand_x][strand_y] = alignment_score_p;
			for (int k = 0; k < alignment_length; k++){
				beta_strand_interactions_max_states[strand_x][strand_y][k] = max_states[k];
			}

			delete[] max_states;

		}
	}

}



void extract_positions_of_beta_residues(char * s, int * beta_residue_positions){

	int j, k, length;

	length = (int)strlen(s);

	k = 0;
	for (j = 0; j < length; j++){

		if (s[j] == 'E')
			beta_residue_positions[k++] = j;

	}

}

double compute_pairwise_alignment_of_two_beta_strands3(int strand_x_start_ind, int strand_y_start_ind, int strand_x_end_ind, int strand_y_end_ind, int strand_x_length, int strand_y_length, int parallel, double ** beta_residue_pairing_probs, int number_of_beta_residues, int ** temp_contact_map){

	int done;
	int nr, nc;
	int i, j, k, m, n;
	double ** max_alignment_scores;
	int ** max_states_2d; //indices of the maximum value coming from previous states
	//int * max_states;
	int alignment_length;

	double score_i, score_j, score_ij;
	double max_score;
	double max_alignment_score;

	n = strand_x_length;
	m = strand_y_length;

	nr = m + 1;
	nc = n + 1;

	allocate_2d_double(max_alignment_scores, nr, nc);
	allocate_2d_int(max_states_2d, nr, nc);

	//max_states = new int[100];
	int max_states[100];

	max_alignment_scores[0][0] = 0.0;

	for (i = 1; i <= n; i++)
		max_states_2d[0][i] = 1;

	for (j = 1; j <= m; j++)
		max_states_2d[j][0] = 2;

	for (i = 1; i <= n; i++){

		for (j = 1; j <= m; j++){

			if (parallel)
				score_ij = max_alignment_scores[j - 1][i - 1] + beta_residue_pairing_probs[strand_y_start_ind + j - 1][strand_x_start_ind + i - 1];

			else score_ij = max_alignment_scores[j - 1][i - 1] + beta_residue_pairing_probs[strand_y_end_ind - (j - 1)][strand_x_start_ind + i - 1];


			max_states_2d[j][i] = 0;
			max_score = score_ij;

			score_i = max_alignment_scores[j][i - 1] + 0.0; //gap score is chosen as 0.0

			if (score_i > max_score){

				max_score = score_i;
				max_states_2d[j][i] = 1;
			}

			score_j = max_alignment_scores[j - 1][i] + 0.0; //gap score is chosen as 0.0

			if (score_j > max_score){

				max_score = score_j;
				max_states_2d[j][i] = 2;
			}

			max_alignment_scores[j][i] = max_score;

		}
	}

	//BACK TRACE THE STATES ON THE MAXIMUM SCORING PATH 
	max_alignment_score = max_alignment_scores[m][n];
	max_states[0] = max_states_2d[m][n];

	i = n;
	j = m;

	done = 0;
	alignment_length = 0;
	k = 1;

	while (!done){

		if (max_states[k - 1] == 0){ //M state

			i--; j--;
			max_states[k] = max_states_2d[j][i];
			alignment_length++;

		}

		else if (max_states[k - 1] == 1){ //Ix state

			i--;
			max_states[k] = max_states_2d[j][i];
			alignment_length++;

		}

		else if (max_states[k - 1] == 2){ //Iy state

			j--;
			max_states[k] = max_states_2d[j][i];
			alignment_length++;

		}

		if ((i == 0) && (j == 0))
			done = 1;

		k++;

	} //end of while (back trace)

	reverse_int_array(max_states, alignment_length);

	i = 0; j = 0;
	for (k = 0; k < alignment_length; k++){

		if (max_states[k] == 0){

			i++; j++;

			if (parallel){

				temp_contact_map[strand_x_start_ind + i - 1][strand_y_start_ind + j - 1] = 1;
				temp_contact_map[strand_y_start_ind + j - 1][strand_x_start_ind + i - 1] = 1;
			}

			else{

				temp_contact_map[strand_x_start_ind + i - 1][strand_y_end_ind - (j - 1)] = 1;
				temp_contact_map[strand_y_end_ind - (j - 1)][strand_x_start_ind + i - 1] = 1;
			}
		}

		else if (max_states[k] == 1){

			i++;
		}

		else if (max_states[k] == 2){

			j++;
		}
	}

	deallocate_2d_double(max_alignment_scores, nr);
	deallocate_2d_int(max_states_2d, nr);
	//delete [] max_states;

	return max_alignment_score;
}


void compute_pairwise_alignments_of_beta_strands3(Protein* protein, int * segment_starts, int * segment_ends, int number_of_beta_strands, double ** beta_residue_pairing_probs, int * beta_residue_positions, int number_of_beta_residues){

	int i, j, k;
	int strand_x, strand_y;
	int strand_x_start, strand_y_start, strand_x_end, strand_y_end;
	int strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind;
	int strand_x_length, strand_y_length;
	int parallel;

	int ** temp_contact_map;
	int ** temp_contact_map2;

	double alignment_score_p, alignment_score_ap;

	allocate_2d_int(temp_contact_map, number_of_beta_residues, number_of_beta_residues);
	allocate_2d_int(temp_contact_map2, number_of_beta_residues, number_of_beta_residues);

	for (i = 0; i < protein->number_of_sheets; i++){
		for (int j = 0; j < protein->sheets[i]->number_of_strands - 1; j++){

			set_2d_int_to_a_value(temp_contact_map, number_of_beta_residues, number_of_beta_residues, 0);
			set_2d_int_to_a_value(temp_contact_map2, number_of_beta_residues, number_of_beta_residues, 0);

			strand_x = protein->sheets[i]->spatial_ordering[j];
			strand_y = protein->sheets[i]->spatial_ordering[j + 1];

			strand_x_start = segment_starts[strand_x];
			strand_y_start = segment_starts[strand_y];

			strand_x_end = segment_ends[strand_x];
			strand_y_end = segment_ends[strand_y];

			strand_x_start_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_x_start);
			strand_y_start_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_y_start);

			strand_x_end_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_x_end);
			strand_y_end_ind = get_matrix_indice(beta_residue_positions, number_of_beta_residues, strand_y_end);

			strand_x_length = strand_x_end - strand_x_start + 1;
			strand_y_length = strand_y_end - strand_y_start + 1;
			if (protein->sheets[i]->interaction_pattern[j] == ANTIPARALLEL){
				//compute anti-parallel alignment
				parallel = 0;
				alignment_score_ap = compute_pairwise_alignment_of_two_beta_strands3(strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind, strand_x_length, strand_y_length, parallel, beta_residue_pairing_probs, number_of_beta_residues, temp_contact_map);
				for (int j = 0; j < strand_x_length; j++){

					for (int k = 0; k < strand_y_length; k++){

						protein->contact_map[strand_x_start_ind + j][strand_y_start_ind + k] = temp_contact_map[strand_x_start_ind + j][strand_y_start_ind + k];
						protein->contact_map[strand_y_start_ind + k][strand_x_start_ind + j] = protein->contact_map[strand_x_start_ind + j][strand_y_start_ind + k];
					}
				}
			}
			else{
				//compute parallel alignment

				parallel = 1;
				alignment_score_p = compute_pairwise_alignment_of_two_beta_strands3(strand_x_start_ind, strand_y_start_ind, strand_x_end_ind, strand_y_end_ind, strand_x_length, strand_y_length, parallel, beta_residue_pairing_probs, number_of_beta_residues, temp_contact_map2);
				for (int j = 0; j < strand_x_length; j++){

					for (int k = 0; k < strand_y_length; k++){

						protein->contact_map[strand_x_start_ind + j][strand_y_start_ind + k] = temp_contact_map2[strand_x_start_ind + j][strand_y_start_ind + k];
						protein->contact_map[strand_y_start_ind + k][strand_x_start_ind + j] = protein->contact_map[strand_x_start_ind + j][strand_y_start_ind + k];
					}
				}
			}

		}
	}

	deallocate_2d_int(temp_contact_map, number_of_beta_residues);
	deallocate_2d_int(temp_contact_map2, number_of_beta_residues);

}



void main(void)
{

	//////////////////////////////////////*conformations*//////////////////////////////////////////
	//////////////////////////////////////definitions////////////////////////////////////////////
	double ** max_alignment_scores;

	char * filename_input = new char[200]; // input filename
	char * filename_prediction1 = new char[200];
	char * filename_prediction2 = new char[200];
	char * fold_path_str = new char[1000];
	char * protein_information_path_str = new char[1000];
	char * fold_list_str = new char[1000];
	char * pdb_id = new char[200];
	char* fold_number_str = new char[3];

	int beta_strand_threshold = 6; // maximum number of strands in each sheet, you can change the threshold	
	int top_selected;
	int number_of_beta_residues; // number of beta residues
	int number_of_beta_strands; // number of beta strands
	int number_of_residues;

	FILE * input_file; // betapro output file
	FILE *contactfile;
	char *contactfileStr = new char[1000];
	double ** pp;

	double **beta_residue_pairing_probs;

	int *** beta_strand_interactions_max_states;

	int max_length = 3000;

	char * R = new char[max_length];
	char * SS = new char[max_length];
	char * s1 = new char[max_length];
	char * s2 = new char[max_length];

	int * segment_starts = new int[100];
	int * segment_ends = new int[100];
	int done;

	FILE* ifp_list;
	FILE * ofp_conformation;
	FILE * ofp_contact_map;

	clock_t all_config_time = 0;
	clock_t start_config_time;
	clock_t start = clock();
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// +                                                                                                  +
	// +                         read protein information from file                                       + 
	// +                                                                                                  +
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	for (int fold_number = 1; fold_number <= 10; fold_number++){
		cout << "fold" << fold_number << "..." <<endl;
		_itoa(fold_number, fold_number_str, 10);
		strcpy(fold_path_str, "Proteins\\");
		strcpy(fold_list_str, fold_path_str);
		strcat(fold_list_str, "fold_lists\\fold");
		strcat(fold_list_str, fold_number_str);
		strcat(fold_list_str, "_list.txt");
		strcat(fold_path_str, "fold");
		strcat(fold_path_str, fold_number_str);
		strcat(fold_path_str, "\\");

		done = 0;
		ifp_list = fopen(fold_list_str, "r");
		int num = 0;
		while (!done){
			read_string_from_file(filename_input, 2000, ifp_list);

			if (filename_input[0] == 'e'){//end of file condition
				done = 1; continue;
			}

			discard_filename_extension(filename_input); //remove the extension (.fasta) to get the protein name 
			strcpy(protein_information_path_str, fold_path_str);
			//strcat(protein_information_path_str,"betapro_predictions\\");//change
			strcat(protein_information_path_str, filename_input);
			strcat(protein_information_path_str, ".out");
			input_file = fopen(protein_information_path_str, "r");
			read_string_from_file(pdb_id, 200, input_file);
			if (strstr(pdb_id, "1A15A")){
				int y = 0;
			}
			//read_string_from_file(s, max_length, input_file);
			//read_string_from_file(s, max_length, input_file);
			read_string_from_file(R, max_length, input_file);
			//read_string_from_file(s, max_length, input_file);
			//read_string_from_file(s, max_length, input_file);
			read_string_from_file(SS, max_length, input_file);
			read_string_from_file(s1, max_length, input_file); //skip 4th line of input file
			read_string_from_file(s1, max_length, input_file); //skip 5th line of input file
			// convert the C characters to L (loop state)
			C2L(SS);

			beta_strand_segment_borders(SS, segment_starts, segment_ends, &number_of_beta_strands,
				&number_of_beta_residues);///////////Protein class

			if (number_of_beta_strands > beta_strand_threshold){
				fclose(input_file);
				continue;
			}


			//cout << pdb_id << "  " << number_of_beta_strands << endl;
			////////read deepLearning pairing probs prediction////////////////
			number_of_residues = (int)strlen(SS);
			int *beta_residue_positions = new int[number_of_beta_residues];
			extract_positions_of_beta_residues(SS, beta_residue_positions);
			allocate_2d_double(pp, number_of_residues, number_of_residues);
			allocate_2d_double(beta_residue_pairing_probs, number_of_beta_residues, number_of_beta_residues);
			float db;
			bool sw_read = false;
			////////////////you can change the address of contact file
			/*if (strstr(pdb_id, "1J98A") || strstr(pdb_id, "1MFWA") || strstr(pdb_id, "1Q9UB") || strstr(pdb_id, "1QGWA")
			|| strstr(pdb_id, "1ESRA") || strstr(pdb_id, "1QYSA") || strstr(pdb_id, "1L8RA")){
			continue;
			}*/

			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// +                                                                                                  +
			// +                              read contact scores from file                                       + 
			// +                                                                                                  +
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			strcpy(contactfileStr, "ultra_deep_learning_contactmap\\");
			strcat(contactfileStr, pdb_id);
			strcat(contactfileStr, ".GCNN");
			contactfile = fopen(contactfileStr, "r");
			if (contactfile != NULL){
				sw_read = true;

				for (int i = 0; i < number_of_residues; i++){
					for (int j = 0; j < number_of_residues; j++){
						fscanf(contactfile, "%f", &db);
						pp[i][j] = (double)db;
					}
				}


				for (int i = 0; i < number_of_beta_residues; i++){
					for (int j = 0; j < number_of_beta_residues; j++){
						beta_residue_pairing_probs[i][j] = pp[beta_residue_positions[i]][beta_residue_positions[j]];
					}
				}
				fclose(contactfile);
			}

			if (sw_read == false){
				for (int i = 0; i < number_of_beta_residues; i++){
					for (int j = 0; j < number_of_beta_residues; j++){
						fscanf(input_file, "%f", &db);
						beta_residue_pairing_probs[i][j] = (double)db;
					}
				}
			}
			fclose(input_file);
			deallocate_2d_double(pp, number_of_residues);
			allocate_2d_double(max_alignment_scores, number_of_beta_strands, 2 * number_of_beta_strands);
			allocate_3d_int(beta_strand_interactions_max_states, number_of_beta_strands, 2 * number_of_beta_strands, 1000);
			////compute alignments/////
			compute_pairwise_alignments_of_beta_strands2(segment_starts, segment_ends, number_of_beta_strands, beta_residue_pairing_probs, beta_residue_positions, number_of_beta_residues, max_alignment_scores, beta_strand_interactions_max_states);
			//cout << pdb_id << "   " << number_of_beta_strands << endl;
			////////////////sum of subsets/////////////////////////////
			int n = number_of_beta_strands;
			int num = 0;
			for (int i = 2; i < n - 1; i++){
				num += n / i;
			}

			//////////////////*creating S array*///////////////////////////////////
			int *S = (int*)(calloc(num, sizeof(int)));
			int kk = 0, j;
			for (int i = 2; i < n - 1; i++){
				for (int j = 0; j < n / i; j++){
					S[kk++] = i;
				}

			}
			///////////////////////////////////////////////////////////////////////
			////////////////////////*building the tree of subsets*//////////////////////////////////////
			int *X = (int*)(calloc(n / 2, sizeof(int)));
			int k = 0; int i = 0;
			int sum = 0;
			int leaf_num = 0;
			SubsetTreeNode root(0, NULL);
			SubsetTreeNode *curr;
			curr = &root;
			while (i < num && X[0] <= n / 2){
				X[k] = S[i];
				if (sum + S[i] < n)
				{
					if (i + 1 < num && sum + S[i] + S[i + 1] <= n){
						curr->children.push_back(new SubsetTreeNode(S[i], curr));
						curr = (curr->children.back());
						sum += X[k];
						k++; i++;

					}
					else{
						i++;
					}
				}
				else if (sum + S[i] > n) {
					k--;
					curr = curr->parent;
					delete(curr->children.back());
					curr->children.pop_back();//////how to delete?
					sum -= X[k];
					for (int j = 0; j <num; j++){//find the minimum number in S that is greater than X[k](a formula can be replaced)
						if (S[j] > X[k]){
							i = j;
							break;
						}
					}
				}
				else {
					//for (int p = 0; p <= k; p++) printf("%d", X[p]); // to be extended for the next step
					//printf("\n");
					k--;
					//curr->children.push_back(new SubsetTreeNode(S[i], curr));//change for top topologies
					curr->children.push_back(new SubsetTreeNode(S[i], curr, leaf_num));
					leaf_num++;
					curr = curr->parent;
					/*if (curr->parent->parent == NULL){
					curr = &root;
					}
					else{
					curr = curr->parent->parent;
					curr->sNode.back().children.push_back(subsetNode(curr));
					curr = &(curr->sNode.back().children.back());
					}*/
					sum -= X[k];
					for (int j = 0; j <num; j++){//find the minimum number in S that is greater than X[k](a formula can be replaced)
						if (S[j] > X[k]){
							i = j;
							break;
						}
					}
				}
			}
			root.children.push_back(new SubsetTreeNode(n, &root, leaf_num));
			leaf_num++;
			top_selected = leaf_num;
			free(S);
			free(X);
			////////////////////*conformations*//////////////////////////////////////////////////////
			int* all_strands = new int[number_of_beta_strands];
			for (int i = 0; i < number_of_beta_strands; i++){
				all_strands[i] = i;
			}
			k = 0;
			int c = 0;
			//double max_score = -1;
			//Protein best_protein(number_of_beta_strands,number_of_beta_residues);
			vector<double> max_scores;
			vector<Protein*> best_proteins;
			for (int i = 0; i < top_selected; i++){
				max_scores.push_back(-1);
				best_proteins.push_back(new Protein(number_of_beta_strands, number_of_beta_residues));
			}
			double score = 0;
			double temp_score;
			X = new int[number_of_beta_strands / 2];
			ConformationTreeNode* conformation_root = NULL;
			ConformationTreeNode *curr_conformation = NULL;
			SubsetTreeNode *curr_subset_node;
			for (int i = 0; i < root.children.size(); i++){
				//cout << score << endl;////////////////////////////////////////
				curr_subset_node = (root.children[i]);
				curr_conformation = new ConformationTreeNode(curr_subset_node->data, number_of_beta_strands, all_strands, curr_conformation);
				c = 0;
				conformation_root = curr_conformation;
				while (curr_conformation != NULL){//??????
					curr_conformation->next_state(max_alignment_scores, number_of_beta_strands);
					if (!curr_conformation->end_of_processing){
						curr_conformation->compute_score(max_alignment_scores, number_of_beta_strands);
						score += curr_conformation->score;
						if (c < curr_subset_node->children.size()){
							X[k++] = c;
							curr_subset_node = (curr_subset_node->children[c]);
							curr_conformation->child = new ConformationTreeNode(curr_subset_node->data, curr_conformation->possible_number_of_strands_for_selection - curr_conformation->node_number_of_beta_strands, curr_conformation->remaining_strands, curr_conformation);
							curr_conformation = curr_conformation->child;
						}
						else if (curr_subset_node->children.empty()){
							if (score>max_scores[curr_subset_node->leaf_number]){
								best_proteins[(curr_subset_node->leaf_number)]->init(conformation_root, score);
								max_scores[(curr_subset_node->leaf_number)] = score;
							}

							score -= curr_conformation->score;//change
							/*Protein p(number_of_beta_strands, number_of_beta_residues);
							p.init(conformation_root, score);
							p.print();*/
						}
						else{
							int y = 0;
						}
					}
					else{
						////////////////////////////////////////////////////////////
						/*if (!curr_subset_node->children.empty()){
						int y = 0;
						score -= curr_conformation->score;
						}*/
						//////////////////////////////////////////////////////////////
						curr_subset_node = curr_subset_node->parent;
						curr_conformation = curr_conformation->parent;
						if (curr_conformation != NULL){
							delete(curr_conformation->child);// ->~ConformationTreeNode();
							k--;
							c = X[k] + 1;
							if (c < curr_subset_node->children.size()){
								X[k++] = c;
								curr_subset_node = (curr_subset_node->children[c]);
								curr_conformation->child = new ConformationTreeNode(curr_subset_node->data, curr_conformation->possible_number_of_strands_for_selection - curr_conformation->node_number_of_beta_strands, curr_conformation->remaining_strands, curr_conformation);
								curr_conformation = curr_conformation->child;
								c = 0;
							}
							else{
								score -= curr_conformation->score;/////change
								c = 0;
							}
						}
					}

				}
				delete(conformation_root);
			}
			///////////*write the best protein and free memory*////////////////////////
			/*int * start_point_for_strands = new int[best_protein.number_of_beta_strands];
			start_point_for_strands[0] = 0;


			for (int i = 1; i<best_protein.number_of_beta_strands; i++){
			start_point_for_strands[i] = start_point_for_strands[i - 1] + (segment_ends[i - 1] - segment_starts[i - 1] + 1);
			}*/

			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			// +                                                                                                  +
			// +                                  write output files                                              + 
			// +                                                                                                  +
			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			char iter_str[10];
			for (int iter = 0; iter < top_selected && iter<best_proteins.size(); iter++){

				best_proteins[iter]->make_contact_map(beta_residue_positions, beta_strand_interactions_max_states, segment_starts, segment_ends);

				strcpy(filename_prediction1, "Predictions\\fold");
				strcat(filename_prediction1, fold_number_str);
				strcpy(filename_prediction2, "Predictions\\fold");
				strcat(filename_prediction2, fold_number_str);
				strcat(filename_prediction1, "\\conformation_"); strcat(filename_prediction1, filename_input); strcat(filename_prediction1, "_");
				_itoa(iter, iter_str, 10);
				strcat(filename_prediction1, iter_str); strcat(filename_prediction1, ".out");
				strcat(filename_prediction2, "\\contact_map_"); strcat(filename_prediction2, filename_input); strcat(filename_prediction2, "_"); strcat(filename_prediction2, iter_str); strcat(filename_prediction2, ".out");

				ofp_conformation = fopen(filename_prediction1, "w");
				ofp_contact_map = fopen(filename_prediction2, "w");
				//print betaza predictions to file
				int strand1, strand2;
				char ch;
				for (int i = 0; i < best_proteins[iter]->number_of_sheets; i++){
					for (int j = 0; j < best_proteins[iter]->sheets[i]->number_of_strands - 1; j++){
						strand1 = best_proteins[iter]->sheets[i]->spatial_ordering[j];
						strand2 = best_proteins[iter]->sheets[i]->spatial_ordering[j + 1];
						fprintf(ofp_conformation, "%d ", strand1 + 1);
						fprintf(ofp_conformation, "%d\n", strand2 + 1);
						if (segment_starts[strand1] == segment_ends[strand1] || segment_starts[strand2] == segment_ends[strand2]){
							ch = 'B';
						}
						else if (best_proteins[iter]->sheets[i]->interaction_pattern[j] == PARALLEL){
							ch = 'P';
						}
						else{
							ch = 'A';
						}
						fprintf(ofp_conformation, "%c\n\n", ch);

					}
				}


				for (int j = 0; j < number_of_beta_residues; j++){

					for (int k = 0; k < number_of_beta_residues; k++){

						fprintf(ofp_contact_map, "%d ", best_proteins[iter]->contact_map[j][k]);

					}

					fputs("\n", ofp_contact_map);
				}
				//delete subset tree?????

				fclose(ofp_contact_map);
				fclose(ofp_conformation);
			}
			freemat(beta_residue_pairing_probs);
			freemat(max_alignment_scores);
			delete[] beta_residue_positions;
			delete[] all_strands;
			//delete[] start_point_for_strands;
			deallocate_3d_int(beta_strand_interactions_max_states, number_of_beta_strands, number_of_beta_strands);
			for (int i = 0; i < top_selected && i<best_proteins.size(); i++){
				delete best_proteins[i];
			}
			best_proteins.clear();
			max_scores.clear();
		}

	}

	delete[] filename_input;
	delete[] filename_prediction1;
	delete[] filename_prediction2;
	delete[] fold_list_str;
	delete[] fold_path_str;
	delete[] contactfileStr;
	delete[] protein_information_path_str;
	delete[] pdb_id;
	delete[] R;
	delete[] SS;
	delete[] s1;
	delete[] s2;
	delete[] segment_starts;
	delete[] segment_ends;
	getchar();

}
