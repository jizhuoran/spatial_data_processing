#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>

#define LANDMARK_NUM 10
#define DEBUG 1
using namespace std;

class Edge
{
	
public:
	int id;
	int endNode;
	double distance;
	Edge(int _id, int _end, double _distance)
		: id(_id), endNode(_end), distance(_distance) {}
	~Edge();
	
};

class Node
{
public:
	/*
		According to the structure of the data, use the id of the node as the hash map.
		the hash function is f(x) = x, in other word, In the node list, for example, named
		node_list, node_list[i] is the node with identifier i.
	*/
	int id; 

	double longitude;
	double latitude;

	std::vector<Edge *> edgeList;

	bool visited;
	int predecessor;
	double SDP;
	double heuristic;
	
	Node(int _id, double _longitude, double _latitude) 
		: id(_id), longitude(_longitude), latitude(_latitude), visited(false), predecessor(-1), SDP(-1.0), heuristic(-1.0) {}
	~Node();

	
};

class Comparator_dijkstra {
public:
    bool operator()(const Node* a, const Node* b)
    {
		if(a->SDP == b->SDP) {
			return (a->id < b->id);
		}
		return (a->SDP < b->SDP);
    }
};


class Comparator_a_star {
public:
    bool operator()(const Node* a, const Node* b)
    {
		if((a->SDP + a->heuristic) == (b->SDP + b->heuristic)) {
			return (a->id < b->id);
		}
		return ((a->SDP + a->heuristic) < (b->SDP + b->heuristic));
    }
};

class RoadNetwork {
	public:
	int numOfNodes_;
	double **s;

	double **spd_set;


	std::vector<Node *> node;

	RoadNetwork() {
		numOfNodes_ = 0;
		spd_set = new double*[LANDMARK_NUM];
	}
	void readRoadNetwork(char * nodepath, char * edgepath);
	void output(char * outpath);
	Node* hashmap(int i) {
		return node[i];
	}

	double dijkstra(int src, int dst, vector<int> & path, int & iteration);
	double a_star(int src, int dst, vector<int> & path, int & iteration);

	
	void dijkstra_tree(int src, double spd_set[]);
	void buildLandMark();
	double a_star_landmark(int src, int dst, vector<int> & path, int & iteration);


};


double heuristic_a_star(Node  *a, double dst_longtitude, double dst_latitude) {
	return sqrt((a->longitude - dst_longtitude) * (a->longitude - dst_longtitude) + (a->latitude - dst_latitude) * (a->latitude - dst_latitude));
}

double heuristic_landmark(int src, int dst, double** spd_set) {
	double max = 0.0;
	double tmp = 0.0;
	for(int i = 0; i < 10; ++i) {
		tmp = abs(spd_set[i][src] - spd_set[i][dst]);
		if(max < tmp) max = tmp;
	}
	return max;
}

bool notLandmark(int node, int number, int landmark[]) {
	for(int i = 0; i < number; ++i) {
		if(node == landmark[i]) {
			return false;
		}
	}
	return true;
}


void RoadNetwork::readRoadNetwork(char * nodepath, char * edgepath) {
	int __id;
	double __longitude, __latitude;
	ifstream node_input;
	node_input.open(nodepath);
	if (node_input.is_open()) {
    	while (node_input >> __id >> __longitude >> __latitude){ // read the node
			Node *tmp = new Node(__id, __longitude, __latitude);
			node.push_back(tmp);
			++numOfNodes_;
  		}
  	}
  	node_input.close();


	int __startID, __endID;
	double __distance;
	ifstream edge_input;
	edge_input.open(edgepath);
	if (edge_input.is_open()) {
    	while (edge_input >> __id >> __startID >> __endID >> __distance){ // read the edge
			/*
				As the edge(id, startNode, endNode, distance) is bidirection, we need 
				to add this edge to the edge list of startNode as well as the endNode.
				The difference is that the another side of the edge is each other.
			*/
 			Edge *tmp = new Edge(__id, __endID, __distance);
			node[__startID]->edgeList.push_back(tmp);
			tmp = new Edge(__id, __startID, __distance);
			node[__endID]->edgeList.push_back(tmp);
  		}
  	}
  	edge_input.close();	
}

void RoadNetwork::output(char * outpath) {
	ofstream output;
  	output.open(outpath);
	if (output.is_open()) {
		output << numOfNodes_ << endl;
    	for (int i = 0; i < numOfNodes_; ++i) {
    		output << node[i]->id << "\t" << node[i]-> longitude << "\t" << node[i]-> latitude;
			for(int j = 0; j < node[i]->edgeList.size(); ++j) {
				output << "\t" << node[i]->edgeList[j]->id << "\t" << node[i]->edgeList[j]->endNode << "\t" << node[i]->edgeList[j]->distance;
			}
			output << endl;
    	}
  	}
  	output.close();
		
}

double RoadNetwork::dijkstra(int src, int dst, vector<int> & path, int & iteration) {
	/*
		initiate the varible.
		The reason why store "visited", "predecessor", "SDP" in node is that
		they only occupy 12 bytes, which is much smaller then other information, 
		and it can save the memory allocation time during search.
	*/
	iteration = 0; 
	for(int i = 0; i < node.size(); ++ i) {
		node[i]->visited = false;
		node[i]->predecessor = -1;
		node[i]->SDP = -1.0;
	}
	/*
		Using ordered set but not hash nor priority queue can make the implatation easier, and 
		do not increase the time-complexity. But we still call it a q to make it easy to understand.
		There is a trick method, if two node with same SDP, we compare them
		by their id.
	*/
	std::set<Node*, Comparator_dijkstra > q; 
	node[src]->SDP = 0;
	q.insert(node[src]);// insert the src node.

	while(!q.empty()) {
		Node *top = *(q.begin());
		++iteration; // we assume that we need to "de-heap" it before knowing the information of the node.

		if(top->id == dst) { // the top node is the dst node.
			Node *tmp = top;
			path.insert(path.begin(),tmp->id); //insert dst to the path.
			while(tmp->predecessor != -1) {
				path.insert(path.begin(),tmp->predecessor); //insert the predecessor to the front of the path
				tmp = node[tmp->predecessor];
			}
			return top->SDP;
		}
			
		q.erase(q.begin());
		top->visited = true;
		for(int i = 0; i < top->edgeList.size(); ++i) { // for each neighbor u of v
			Node *u = node[top->edgeList[i]->endNode];
			if(!(u->visited)) { // if u is unvisited
				if(u->SDP < 0 || u->SDP > top->SDP + top->edgeList[i]->distance) {
					/*
						erase the old u and insert it again to update the node information and position.
						if the node is non-existent, q.erase do nothing, in this situation, the result is
						just inserting this node.
					*/
					q.erase(u);
					node[top->edgeList[i]->endNode]->SDP = top->SDP + top->edgeList[i]->distance;
					u->predecessor = top->id;
					q.insert(u);
				}
			}
				
		}
	}
	return -1.0; // not find this node.
}

double RoadNetwork::a_star(int src, int dst, vector<int> & path, int & iteration) {
	double dst_longtitude = node[dst]->longitude;
	double dst_latitude = node[dst]->latitude;

	/*
		initiate the varible.
		The reason why store "visited", "predecessor", "SDP" and "heuristic" in node is that
		they only occupy 12 bytes, which is much smaller then other information, 
		and it can save the memory allocation time during search.
	*/

	iteration = 0;
	for(int i = 0; i < node.size(); ++ i) {
		node[i]->visited = false;
		node[i]->predecessor = -1;
		node[i]->heuristic = 0.0;
		node[i]->SDP = -1.0;
	}
	/*
		Using ordered set but not hash nor priority queue can make the implatation easier, and 
		do not increase the time-complexity. But we still call it a q to make it easy to understand.
		There is a trick method, if two node with same SDP, we compare them
		by their id.
	*/
	std::set<Node*, Comparator_a_star > q;
	node[src]->SDP = 0;
	node[src]->heuristic = heuristic_a_star(node[src],dst_longtitude ,dst_latitude);
	q.insert(node[src]);

	while(!q.empty()) {
		Node *top = *(q.begin());//Find the node with least SDP+heuristic
		++iteration;// we assume that we need to "de-heap" it before knowing the information of the node.
		if(top->id == dst) {// the top node is the dst node.
			Node *tmp = top;
			path.insert(path.begin(),tmp->id);
			while(tmp->predecessor != -1) {
				path.insert(path.begin(),tmp->predecessor);
				tmp = node[tmp->predecessor];
			}
			return top->SDP;
		}
		q.erase(q.begin());
		top->visited = true;
		for(int i = 0; i < top->edgeList.size(); ++i) {
			Node *u = node[top->edgeList[i]->endNode];
			if(!(u->visited)) {
				// we check whether the SDP but not the SDP+heuristic
				if(u->SDP < 0 || u->SDP > top->SDP + top->edgeList[i]->distance) {
					q.erase(u);
					if(u -> SDP < 0) { // the node is first add into the queue.
						u->heuristic = heuristic_a_star(u,dst_longtitude ,dst_latitude);
					}
					node[top->edgeList[i]->endNode]->SDP = top->SDP + top->edgeList[i]->distance;
					u->predecessor = top->id;
					q.insert(u);
				}
			}
				
		}
	}
	return -1.0;
}

void  RoadNetwork::dijkstra_tree(int src, double spd_set[]) {
	/*
		This function is same as the dijkstra function, the different is that it not find a 
		particular destination node but give the tree. 
	*/
	for(int i = 0; i < node.size(); ++ i) {
		node[i]->visited = false;
		node[i]->predecessor = -1;
		node[i]->SDP = -1.0;
	}	
	std::set<Node*, Comparator_dijkstra > q;
	node[src]->SDP = 0;
	q.insert(node[src]);

	while(!q.empty()) {
		Node *top = *(q.begin());	
		spd_set[top->id] = top->SDP;// when the node pop out, we can determine the SDP	
		q.erase(q.begin());
		top->visited = true;
		for(int i = 0; i < top->edgeList.size(); ++i) {
			Node *u = node[top->edgeList[i]->endNode];
			if(!(u->visited)) {
				if(u->SDP < 0 || u->SDP > top->SDP + top->edgeList[i]->distance) {
					q.erase(u);
					node[top->edgeList[i]->endNode]->SDP = top->SDP + top->edgeList[i]->distance;
					u->predecessor = top->id;
					q.insert(u);
				}
			}
				
		}
	}
}

void RoadNetwork::buildLandMark() {
	int landmark[LANDMARK_NUM];
	double max = 0.0;
	double tmp = 0.0;
	for (int i = 0; i < LANDMARK_NUM; ++i) {
		spd_set[i] = new double[numOfNodes_];
		for (int j = 0; j < numOfNodes_; ++j) {
			spd_set[i][j] = -1.0;
		}
	}

	dijkstra_tree(0, spd_set[0]);// choose the v0

	max = 0.0;
	//choose the furthest node to v0, which is same as choose the furthest from v0.
	for(int i = 1; i < numOfNodes_; ++i) {
		if(spd_set[0][i] > max) {
			max = tmp;
			landmark[0] = i;
		}	
	}
	for(int i = 1; i < LANDMARK_NUM; ++i) {
		max = 0.0;
		dijkstra_tree(landmark[i-1], spd_set[i-1]);//calculate the v_i-1, since v1, v2 ...v_i-2 have been stored already.
	
for(int j = 0; j < numOfNodes_; ++j) {
			tmp = 0.0;

			/*
				calculate SPD(v1, vi) + SPD(v2, vi) + ... + SPD(v_i-1, vi)
			*/
			for (int k = 0; k < i; ++k) {
				tmp += spd_set[k][j];
				if(spd_set[k][j] < 0.0) {
					tmp = -1;
					break;
				}
			}

			//make sure the new node is not alread a landmark
			if(tmp > max && notLandmark(j, i, landmark)) {
				max = tmp;
				landmark[i] = j;
			}
		}
	}
	dijkstra_tree(landmark[LANDMARK_NUM-1], spd_set[LANDMARK_NUM-1]);//calculate the last spd_set
		
}

double RoadNetwork::a_star_landmark(int src, int dst, vector<int> & path, int & iteration) {
	/*
	initiate the varible.
	The reason why store "visited", "predecessor", "SDP" and "heuristic" in node is that
	they only occupy 12 bytes, which is much smaller then other information, 
	and it can save the memory allocation time during search.
	*/
	iteration = 0;
	for(int i = 0; i < node.size(); ++ i) {
		node[i]->visited = false;
		node[i]->predecessor = -1;
		node[i]->heuristic = 0.0;
		node[i]->SDP = -1.0;
	}
	/*
	Using ordered set but not hash nor priority queue can make the implatation easier, and 
	do not increase the time-complexity. But we still call it a q to make it easy to understand.
	There is a trick method, if two node with same SDP, we compare them
	by their id.
	*/
	std::set<Node*, Comparator_a_star > q;
	node[src]->SDP = 0;
	node[src]->heuristic = heuristic_landmark(src, dst, spd_set);
	q.insert(node[src]);

	while(!q.empty()) {
		Node *top = *(q.begin());//Find the node with least SDP+heuristic
		++iteration;// we assume that we need to "de-heap" it before knowing the information of the node.
		if(top->id == dst) {
			Node *tmp = top;
			path.insert(path.begin(),tmp->id);
			while(tmp->predecessor != -1) {
				path.insert(path.begin(),tmp->predecessor);
				tmp = node[tmp->predecessor];
			}
			return top->SDP;
		}
		q.erase(q.begin());
		top->visited = true;
		for(int i = 0; i < top->edgeList.size(); ++i) {
			Node *u = node[top->edgeList[i]->endNode];
			if(!(u->visited)) {
				if(u->SDP < 0 || u->SDP > top->SDP + top->edgeList[i]->distance) {
					q.erase(u);
					if(u -> SDP < 0) {
						u->heuristic = heuristic_landmark(u->id, dst, spd_set);
					}
					node[top->edgeList[i]->endNode]->SDP = top->SDP + top->edgeList[i]->distance;
					u->predecessor = top->id;
					q.insert(u);
				}
			}
				
		}
	}
	return -1.0;
}
