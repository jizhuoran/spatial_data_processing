#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>

#define AGE_RANGE 81
#define BINS_NUMBER 8
#define THRESHOLD_NUMBER 7

using namespace std;

double estimate_eqWidth(int LEFT, int RIGHT, char * histogram_path)
{
	unsigned int person_of_bin[BINS_NUMBER];
	for (int i = 0; i < BINS_NUMBER; ++i) {
		person_of_bin[i] = 0;
	}
	//get the input
	ifstream input;
	input.open(histogram_path);
	if (input.is_open()) {
		for(int i = 0; i < BINS_NUMBER; ++i) {
			input >> person_of_bin[i];
		}
  	} else {
  		cerr << "Can not open the file: " << histogram_path << endl;
  		return 0.0;
  	}
  	input.close();

	double sum = 0;

	int thresholds[THRESHOLD_NUMBER + 2] = {0,10,20,30,40,50,60,70,81};

	int ll = -1;
	int rr = -1;

	//find the largest thresholds before LEFT
	for (int i = 0; i < THRESHOLD_NUMBER + 1; ++i) {
		if(thresholds[i+1] > LEFT) {
			ll = i;
			break;
		}
	}
	//find the smallest thresholds after RIGHT
	for (int i = 0; i < THRESHOLD_NUMBER + 2; ++i) {
		if(thresholds[i] > RIGHT) {
			rr = i;
			break;
		}
	}

	// if LEFT and RIGHT in the same interval
	if((ll + 1) == rr) {
		return 1.0 * person_of_bin[ll] * (RIGHT-LEFT+1)/(thresholds[rr] - thresholds[ll]);
	}

	// sum the interval which totally between LEFT and RIGHT
	for (int i = ll + 1; i < rr-1 ; ++i) {
		sum += person_of_bin[i];
	}

	//add sum with right side of LEFT and its right and RIGHT and its left.
	sum += 1.0 * person_of_bin[ll]*(thresholds[ll+1] - LEFT) / (thresholds[ll+1] - thresholds[ll]);
	sum += 1.0 * person_of_bin[rr - 1]*(RIGHT - thresholds[rr-1] + 1) / (thresholds[rr] - thresholds[rr-1]);

	return sum;
}

double estimate_eqDepth(int LEFT, int RIGHT, char * histogram_path)
{
	unsigned int thresholds[THRESHOLD_NUMBER+2];
	int total = 0;

	for (int i = 0; i < THRESHOLD_NUMBER+2; ++i) {
		thresholds[i] = 0;
	}
	thresholds[THRESHOLD_NUMBER + 1] = 81;
	// get the input
	ifstream input;
	input.open(histogram_path);
	if (input.is_open()) {
		input>>total;
		for(int i = 1; i < THRESHOLD_NUMBER+1; ++i) {
			input >> thresholds[i];
		}
  	} else {
  		cerr << "Can not open the file: " << histogram_path << endl;
  		return 0.0;
  	}
  	input.close();

	double record_per_bin = total / BINS_NUMBER;
	double sum = 0;

	int ll = -1;
	int rr = -1;
	//find the largest thresholds before LEFT
	for (int i = 0; i < THRESHOLD_NUMBER+1; ++i) {
		if(thresholds[i+1] > LEFT) {
			ll = i;
			break;
		}
	}
	//find the smallest thresholds after RIGHT
	for (int i = 0; i < THRESHOLD_NUMBER+2; ++i) {
		if(thresholds[i] > RIGHT) {
			rr = i;
			break;
		}
	}
	// if LEFT and RIGHT in the same interval
	if((ll + 1) == rr) {
		return record_per_bin * (RIGHT-LEFT+1)/(thresholds[rr] - thresholds[ll]);
	}

	// sum the interval which totally between LEFT and RIGHT
	sum += (rr-ll-2) * record_per_bin;

	//add sum with right side of LEFT and its right and RIGHT and its left.
	sum += record_per_bin*(thresholds[ll+1] - LEFT) / (thresholds[ll+1] - thresholds[ll]);
	sum += record_per_bin*(RIGHT - thresholds[rr-1] + 1) / (thresholds[rr] - thresholds[rr-1]);

	return sum;
}

int get_result(int LEFT, int RIGHT, char * dat_path)
{
	string line, age_in_string;

	//person_per_age is use to store the number of records in that age
	unsigned int person_per_age[AGE_RANGE];
	for (int i = 0; i < AGE_RANGE; ++i) {
		person_per_age[i] = 0;
	}


	ifstream input;
	input.open(dat_path);
	if (input.is_open()) {
    	while (getline (input,line)){
 			age_in_string = strtok(&(line[0])," ");
 			age_in_string = strtok(NULL," ");//get the second column of the data
 			++person_per_age[atoi(&(age_in_string[0]))];//add 1 to the corresponding elements of the array
  		}
  	} else {
  		cerr << "Can not open the file: " << dat_path << endl;
  		return 0;
  	}
  	input.close();
	unsigned int sum = 0;

	for(int i = LEFT; i <= RIGHT; ++ i) {
		sum += person_per_age[i];
	}
	return sum;
}

int main(int argc, char** argv)
{
	if (argc != 6){
		cerr << "Usage: " << argv[0] << " LEFT RIGHT WIDTH_HISTOGRAM_PATH DEPTH_HISTOGRAM_PATH DATA_PATH" << endl;
		/*
			LEFT(int): the lower bound of the interval
			RIGHT(int): the upper bound of the interval
			WIDTH_HISTOGRAM_PATH(char *): the file path of the equal-width histogram
			DEPTH_HISTOGRAM_PATH(char *): the file path of the equal-depth histogram
			DATA_PATH(char *): the file path of final_general.dat
  		*/
		return -1;
	}
	cout << estimate_eqWidth(atoi(argv[1]), atoi(argv[2]), argv[3]) << endl;
	cout << estimate_eqDepth(atoi(argv[1]), atoi(argv[2]), argv[4]) << endl;
	cout << get_result(atoi(argv[1]), atoi(argv[2]), argv[5]) << endl;
	return 0;
}
