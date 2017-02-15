#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <string.h>
#include <cmath>
#include <set>

using namespace std;



class Bound_node {
public:
    int id;
    int count;
    double bound;
    Bound_node(int i, double l) : count(1), id(i), bound(l) {};
};

class Edge
{
    
public:
    int id;
    int endNode;
    double distance;
    Edge(int _id, int _end, double _distance)
    : id(_id), endNode(_end), distance(_distance) {}
    ~Edge() {};
    
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
    ~Node() {
        for(int i = 0; i < edgeList.size(); ++i) {
            delete edgeList[i];
        }
    }
    
    
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



class Comparator_bound {
public:
    bool operator()(const Bound_node* a, const Bound_node* b)
    {
        if(a->bound == b->bound) {
            return (a->id < b->id);
        }
        return (a->bound < b->bound);
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

double heuristic_a_star(Node  *a, double dst_longtitude, double dst_latitude) {
    return sqrt((a->longitude - dst_longtitude) * (a->longitude - dst_longtitude) + (a->latitude - dst_latitude) * (a->latitude - dst_latitude));
}


class Dijkstra_tree {
public:
    int id;
    std::set<Node*, Comparator_dijkstra > q;
    size_t numOfNodes_;
    std::vector<Node *> node;
    
    Dijkstra_tree(int i, std::vector<Node *> __node, int src) {
        numOfNodes_ = __node.size();
        id = i;
        for(int i = 0; i < numOfNodes_; ++i) {
            Node *tmp = new Node(__node[i]->id, __node[i]->longitude, __node[i]->latitude);
            for(int j = 0; j < __node[i]->edgeList.size(); ++j) {
                Edge *tmp_edge = new Edge(__node[i]->edgeList[j]->id, __node[i]->edgeList[j]->endNode, __node[i]->edgeList[j]->distance);
                tmp->edgeList.push_back(tmp_edge);
            }
            tmp->visited = false;
            tmp->predecessor = -1;
            tmp->SDP = -1.0;
            node.push_back(tmp);
        }
        node[src]->SDP = 0;
        q.insert(node[src]);
    }
    
    ~Dijkstra_tree() {
        for(int i = 0; i < numOfNodes_; ++i) {
            delete node[i];
        }
    }
    
    
    
    pair<int, double> getNext() {
        if(!q.empty()) {
            Node *top = *(q.begin());
            q.erase(q.begin());
            top->visited = true;
            for(int i = 0; i < top->edgeList.size(); ++i) {
                Node *u = node[top->edgeList[i]->endNode];
                if(!(u->visited)) {
                    if(u->SDP < 0 || u->SDP > top->SDP + top->edgeList[i]->distance) {
                        q.erase(u);
                        node[top->edgeList[i]->endNode]->SDP = top->SDP + top->edgeList[i]->distance;
                        q.insert(u);
                    }
                }
            }
            return make_pair(top->id,top->SDP);
        }else {
            cerr << "No more node" << endl;
            return make_pair(-1, -1);
        }
    }
};



/*class Dijkstra_node {
public:
    int id;
    Dijkstra_tree* tree;
    Dijkstra_node(int i, Dijkstra_tree* t) : id(i), tree(t) {}
};*/



class RoadNetwork{
public:
    int numOfNodes_;
    std::vector<Node *> node;
    std::vector<Dijkstra_tree *> tree_list;
    RoadNetwork() {
        numOfNodes_ = 0;
    }
    void readRoadNetwork(char * nodepath, char * edgepath);//same as a2
    
    void output(char * outpath);//same as a2
    
    pair<int, double> getNext(int u);
    
    double getSPD(int u, int v);//same as a2

    vector<int> TA(vector<int> inputs, int Aggregation);

    vector<int> NRA(vector<int> inputs, int Aggregation);
    
    vector<int> topK(vector<int> inputs, int Algorithm, int Aggregation)
    {
        if(Algorithm == 1) {
            return TA(inputs, Aggregation);
        }else if(Algorithm == 2) {
            return NRA(inputs, Aggregation);
        }else {
            cerr << "Unexpected Aggregation" << endl;
            return vector<int> ();
        }
        return vector<int> ();
    }
    
};



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

pair<int, double> RoadNetwork::getNext(int u) {
	/*
		As the number of nodes in the query set is usually small
		and the number of total nodes is uauslly large,
		it is better use a vector to store the node we want to 
		preserve the state. 
	*/
    for (int i = 0; i < tree_list.size(); ++i) {
        if (u == tree_list[i]->id) {
            return tree_list[i]->getNext();
        }
    }
    //if it is not exist, create a new one.
    Dijkstra_tree* tmp = new Dijkstra_tree(u, node, u);
    tree_list.push_back(tmp);
    return tmp->getNext();
}

double RoadNetwork::getSPD(int u, int v) {
    double dst_longtitude = node[v]->longitude;
    double dst_latitude = node[v]->latitude;
    
    /*
     initiate the varible.
     The reason why store "visited", "predecessor", "SDP" and "heuristic" in node is that
     they only occupy 12 bytes, which is much smaller then other information,
     and it can save the memory allocation time during search.
     */
    
    for(int i = 0; i < node.size(); ++ i) {
        node[i]->visited = false;
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
    node[u]->SDP = 0;
    node[u]->heuristic = heuristic_a_star(node[u],dst_longtitude ,dst_latitude);
    q.insert(node[u]);
    
    while(!q.empty()) {
        Node *top = *(q.begin());//Find the node with least SDP+heuristic
        if(top->id == v) {// the top node is the dst node.
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
                    q.insert(u);
                }
            }
            
        }
    }
    return -1.0;
}

vector<int> RoadNetwork::TA(vector<int> inputs, int Aggregation)
{
    std::vector<int> hash;
    /*
     hash: like a bit table, 0 means not found, 1 means found, if v has already
     been found, we do not need ti calculated it again.
     */
    for(int i = 0; i < numOfNodes_; ++i) {
        hash.push_back(0);
    }
    double threshold;
    
    std::set<Bound_node*, Comparator_bound > node_list;
    
    //current is the nodes found in this round.
    std::vector<std::pair<int, double> > current;
    for(int i = 0; i < inputs.size(); ++i) {
        current.push_back(make_pair(0,0));
    }
    
    for(int i = 0; i < numOfNodes_; ++i) {
        threshold = 0;
        for(int j = 0; j < inputs.size(); ++j) {
            current[j] = getNext(inputs[j]);
            //update the threshold according to the aggregation.
            if(Aggregation == 1) {
                threshold += current[j].second;
            } else if(Aggregation == 2) {
                if(threshold < current[j].second) {
                    threshold = current[j].second;
                }
            } else {
                cerr << "Unexpected Aggregation" << endl;
                return vector<int>();
            }
        }
        
        for(int j = 0; j < inputs.size(); ++j) {
            if(hash[current[j].first] == 0) {
                double tmp = current[j].second;
                for(int k = 0; k < inputs.size(); ++k) {
                    if(j != k) {
                        //get the aggregation
                        if(Aggregation == 1) {
                            tmp += getSPD(current[j].first, inputs[k]);
                        } else if(Aggregation == 2) {
                            double dis = getSPD(current[j].first, inputs[k]);
                            if(tmp < dis) {
                                tmp = dis;
                            }
                        } else {
                            cerr << "Unexpected Aggregation" << endl;
                            return vector<int>();
                        }
                        
                    }
                }
                hash[current[j].first] = 1;
                Bound_node *tmpnode = new Bound_node(current[j].first, tmp);
                node_list.insert(tmpnode);
            }
        }
        //if the best result is less than the threshold, then we can terminated.
        //the reason why it is not <= is that we want to get all node in the best set.
        //the new found node in the next round can still in this set.
        if ((*node_list.begin())->bound < threshold || i == (numOfNodes_ - 1)) {
            vector<int> result;
            double cutline = (*node_list.begin())->bound;
            while(!node_list.empty() && (*node_list.begin())->bound == cutline){
                result.push_back((*node_list.begin())->id);
                node_list.erase(node_list.begin());
            }

            for(int i = 0; i < tree_list.size(); ++i) {
                delete tree_list[i];
            }
            tree_list.clear();
            return result;
        }
    }
    //should never be there
    for(int i = 0; i < tree_list.size(); ++i) {
        delete tree_list[i];
    }
    tree_list.clear();
    return vector<int>();
}

vector<int> RoadNetwork::NRA(vector<int> inputs, int Aggregation)
{

    std::set<Bound_node*, Comparator_bound > up;
    std::set<Bound_node*, Comparator_bound > lb;
    
    std::vector<Bound_node *> hash;
    
    for(int i = 0; i < numOfNodes_; ++i) {
        hash.push_back(NULL);
    }
    
    for(int i = 0; i < numOfNodes_; ++i) {
        
        for(int j = 0; j < inputs.size(); ++j) {
            pair<int, double> tmp = getNext(inputs[j]);
            if(hash[tmp.first] == NULL) {
                hash[tmp.first] = new Bound_node(tmp.first, tmp.second);
            } else {
                //update the lower bound of the node v
                lb.erase(hash[tmp.first]);
                if(Aggregation == 1) {
                    hash[tmp.first]->bound += tmp.second;
                } else if(Aggregation == 2) {
                    if(hash[tmp.first]->bound < tmp.second) {
                        hash[tmp.first]->bound = tmp.second;
                    }
                } else {
                    cerr << "Unexpected Aggregation" << endl;
                    return vector<int>();
                }
                lb.insert(hash[tmp.first]);
                
                /*
                 count the times that v is found, if the count == the size of the
                 inputs, then we get its upper bound.
                 */
                hash[tmp.first]->count++;
                if(hash[tmp.first]->count == inputs.size()) {
                    up.insert(hash[tmp.first]);
                }
            }
        }
        
        if(up.size() > 0) {
            if ((*up.begin())->bound < (*lb.begin())->bound || i == (numOfNodes_ - 1)) {
                double cutline = (*up.begin())->bound;
                vector<int> result;
                while(!up.empty() && (*up.begin())->bound == cutline){
                    result.push_back((*up.begin())->id);
                    up.erase(up.begin());
                }
                for(int i = 0; i < tree_list.size(); ++i) {
                    delete tree_list[i];
                }
                tree_list.clear();
                return result;
            }
        }
    }
    //should never be there
    for(int i = 0; i < tree_list.size(); ++i) {
        delete tree_list[i];
    }
    tree_list.clear();
    return vector<int>();
}
