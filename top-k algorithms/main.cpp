#include "roadnetwork.cpp"
#include <vector>
#include <set>
#define DEBUG 0

#if DEBUG
#include <ctime>
#include <cstdlib>
#endif
using namespace std;


#define TA 1
#define NRA 2

#define SUM 1
#define MAX 2


int main(int argc, char ** argv)
{


	// You can edit this main function to test your roadnetwork class
	RoadNetwork rn;
	rn.readRoadNetwork(argv[1], argv[2]);
	rn.output(argv[3]);
#if DEBUG
	vector<int> topk1;
	vector<int> topk2;
	srandom(3323);
	std::clock_t start;
    double duration;
    double tatime = 0;
    double nratime = 0;


    for (int i = 0; i < 50; ++i) {
    	int base = random() % 10000;
    	std::set<int> s;
    	while(s.size() < 4) {
    		int tmp = (random() % 10)*1000 + base;
    		s.insert(tmp);
    	}
    	std::vector<int> v;
    	for (int i = 0; i < 3; ++i) {
    		v.push_back(*(s.begin()));
    		s.erase(s.begin());
    	}
    	start = std::clock();
    	topk1 = rn.topK(v, TA, SUM);
    	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    	tatime += duration;
    	cout << "TA" << duration << endl;

    	start = std::clock();
    	topk2 = rn.topK(v, NRA, SUM);
    	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    	nratime += duration;
    	cout << "NRA" << duration << endl;
    }

    cout << "average of ta is " << tatime / 50 << endl;
	cout << "average of nra is " << nratime / 50 << endl;    
#endif	
	vector<int> inputs;
	inputs.push_back(1);
	inputs.push_back(22);
	inputs.push_back(55);
	vector<int> topk;
	topk = rn.topK(inputs, TA, SUM);
#if DEBUG
	cout << topk.size() << endl;
	cout << topk[0] <<endl;
#endif
	return 0;
}