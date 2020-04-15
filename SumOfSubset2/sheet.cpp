#include"sheet.h"
#include <stdio.h>


Sheet::Sheet(int ns){
	number_of_strands = ns;
	spatial_ordering = new int[ns];
	interaction_pattern = new bool[ns - 1];
}

Sheet::Sheet(int ns, int *spatial_ordering_arr, bool *interaction_pattern_arr){
	number_of_strands = ns;
	spatial_ordering = new int[ns];
	interaction_pattern = new bool[ns - 1];
	init(spatial_ordering_arr, interaction_pattern_arr);
}

void Sheet::init(int *spatial_ordering, bool *interaction_pattern){

	for (int i = 0; i<number_of_strands; i++){//mishe behtaresh kard
		if (i<number_of_strands - 1){
			this->interaction_pattern[i] = interaction_pattern[i];
		}
		this->spatial_ordering[i] = spatial_ordering[i];
	}
}



Sheet::~Sheet(){
	delete[] spatial_ordering;
	delete[] interaction_pattern;
}
