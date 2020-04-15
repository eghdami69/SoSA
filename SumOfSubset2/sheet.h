#ifndef SHEET_H
#define SHEET_H

class Sheet
{
private:

public:
	int number_of_strands;
	int* spatial_ordering;
	bool* interaction_pattern;

	Sheet(int ns);
	Sheet(int ns, int *spatial_ordering, bool *interaction_pattern);
	void init(int *spatial_ordering, bool *interaction_pattern);
	~Sheet();
};

#endif