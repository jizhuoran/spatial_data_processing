#include "roadnetwork.cpp"
#include <ctime>
#include <cstdlib>
#define EXPERIMENT 0

int main(int argc, char ** argv)
{
	// You can edit this main function to test your roadnetwork class
	RoadNetwork rn;

	rn.readRoadNetwork(argv[1], argv[2]);
	
	rn.output(argv[3]);
#if EXPERIMENT
	rn.buildLandMark();

	ofstream output;
	output.open("output.txt");
	srand (time(NULL));
	
	unsigned int dij_average = 0;
	unsigned int a_star_average = 0;
	unsigned int a_star_landmark_average = 0;

	int dij_best = 0;
	int a_star_best = 0;
	int a_star_landmark_best = 0;

	for (int i = 0; i < 30000; ++ i) {
		cout << "--------------" << i << "------------"<<endl;
		int s = rand() % rn.numOfNodes_;
		int t = rand() % rn.numOfNodes_;
		vector<int> v1, v2, v3;
		int iteration1, iteration2, iteration3;
		double dij = rn.dijkstra(s, t, v1, iteration1);
		double a_star = rn.a_star(s, t, v2, iteration2);
		double a_star_landmark = rn.a_star_landmark(s, t, v3, iteration3);
		
		dij_average += iteration1;
		a_star_average += iteration2;
		a_star_landmark_average += iteration3;

		if(iteration1 < iteration2 && iteration1 < iteration3) dij_best++;
		if(iteration2 < iteration1 && iteration2 < iteration3) a_star_best++;
		if(iteration3 < iteration2 && iteration3 < iteration1) a_star_landmark_best++;


		output<< "a_star_landmark: " << iteration1 <<"    a_star: "<<iteration2 << "    a_star_landmark: " << iteration3 << endl;
	}
	

	cout << "dij_average: "<< dij_average / 30000 << endl;
	cout << "a_star_average: "<< a_star_average / 30000 << endl;
	cout << "a_star_landmark_average: "<< a_star_landmark_average / 30000 << endl;

	cout << "dij_best: "<< dij_best<< endl;
	cout << "a_star_best: "<< a_star_best<< endl;
	cout << "a_star_landmark_best: "<< a_star_landmark_best<< endl;
#endif

	return 0;
}
