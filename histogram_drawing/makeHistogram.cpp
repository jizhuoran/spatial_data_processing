#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

#define BINS_NUMBER 8
#define AGE_RANGE 81
#define THRESHOLD_NUMBER 7


using namespace std;

void create_histogram(char * data_path, char * eq_width_path, char * eq_depth_path)
{
	// To create an equi-width histogram and an equi-depth histogram, and save them to files.
	// Each output file should contain exactly eight lines, and each line should contain a single integer.

	string line, age_in_string;

	//person_per_age is use to store the number of records in that age
	unsigned int person_per_age[AGE_RANGE];
	for (int i = 0; i < AGE_RANGE; ++i) {
		person_per_age[i] = 0;
	}


	ifstream input;
	input.open(data_path);
	if (input.is_open()) {
    	while (getline (input,line)){
 			age_in_string = strtok(&(line[0])," ");
 			age_in_string = strtok(NULL," ");//get the second column of the data
 			++person_per_age[atoi(&(age_in_string[0]))];//add 1 to the corresponding elements of the array
  		}
  	}else {
  		cerr << "Can not open the file: " << data_path << endl;
  		return;
  	}
  	input.close();


  	unsigned int equi_width[BINS_NUMBER];

  	unsigned int sum = person_per_age[AGE_RANGE - 1];
  	for (int i = 0; i < BINS_NUMBER; ++i) {
		equi_width[i] = 0;
	}


	for (int i = 0; i < AGE_RANGE - 1; ++i) {//from age 0 to age 79
		equi_width[i/10] += person_per_age[i];//[10i, 10(i+1) )
		sum += person_per_age[i];
	}

	equi_width[BINS_NUMBER-1] += person_per_age[AGE_RANGE - 1]; // add age 80 to [70,80]

	//output the result
  	ofstream equi_width_output;
  	equi_width_output.open(eq_width_path);
	if (equi_width_output.is_open()) {
    	for (int i = 0; i < BINS_NUMBER; ++i) {
    		equi_width_output << equi_width[i] <<endl;
    	}
  	}else {
  		cerr << "Can not open the file: " << eq_width_path << endl;
  		return;
  	}
  	equi_width_output.close();



	int threshold[THRESHOLD_NUMBER];
	for (int i = 0; i < THRESHOLD_NUMBER; i++) {
		threshold[i] = 0;
	}

	unsigned int person_per_interval = sum / BINS_NUMBER;

  	int pointer = 1;// pointer to a threshold value, 1 to a1, 2 to a2...

  	int sum_until_now = 0;//sum the number of records until this age

  	for (int i = 0; i < AGE_RANGE; ++i) {
		sum_until_now += person_per_age[i];
		unsigned int pointer_range = person_per_interval * pointer;//expect total number of records from 0 to this records
  		if (sum_until_now >= pointer_range) {
				//if the pointer_range is in some age interval, then we need to determine
				//whether we should include this age
				//the method I use is that calculating the difference with this age and without this age
				//call them ld and rd
				//if ld is smaller, without this age, otherwise, include it.
  		 		unsigned int rd = abs(sum_until_now -pointer_range);
  		 		unsigned int ld = abs(pointer_range - (sum_until_now - person_per_age[i]));
  		 		if (ld >= rd) {//this reason why ld ==rd, we still include it is that
					//when calculating the person_per_interval, we drop the reminder
					//and the experience also show that it is a good choice
  		 			threshold[pointer - 1] = i + 1;
  		 		} else {
  		 			threshold[pointer - 1] = i;
					//this is question that if I do not include it,
					//whether I should minus it from the sum, which means that follow should be include
					/*
						sum_until_now -= person_per_age[i];
						i--;
					*/
					//this is not a good choice since if some age with too many records, there will be empty interval with is misleading
  		 		}

				if (pointer == THRESHOLD_NUMBER) {
					break;
				}

  		 		++pointer;

  		 	}
  	}


	//output the result
	ofstream eq_depth_output;
  	eq_depth_output.open(eq_depth_path);
	if (eq_depth_output.is_open()) {
		eq_depth_output << sum << endl;
    	for (int i = 0; i < THRESHOLD_NUMBER; ++i) {
    		eq_depth_output << threshold[i] <<endl;
    	}
  	} else {
  		cerr << "Can not open the file: " << eq_depth_path << endl;
  		return;
  	}
  	eq_depth_output.close();

}

int main(int argc, char ** argv)
{
	if (argc != 4){
		cerr << "Ussage: " << argv[0] << " DATA_PATH EQ_WIDTH_PATH EQ_DEPTH_PATH" << endl;
		/*
			DATA_PATH(char *): the file path of final_general.dat
			EQ_WIDTH_PATH(char *): the output file path of the equal-width histogram
  			EQ_DEPTH_PATH(char *): the output file path of the equal-depth histogram
  		*/
		return -1;
	}
	create_histogram(argv[1], argv[2], argv[3]);
	return 0;
}
