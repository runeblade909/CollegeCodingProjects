#ifndef GRAPH_H
#define GRAPH_H
#include<vector>
#include<iostream>

struct vertex;

struct adjVertex{
    vertex *v;
};

struct vertex{
    int key;                    //Data
    bool visited = false;       //Visited  
    int distance = 0;           //Weight
    std::vector<adjVertex> adj; //Pointer list to adjacent items
};

class Graph
{
    public:
        void addEdge(int v1, int v2);
        void addVertex(int v);
        bool isBridge(int key1, int key2);
        void removeEdge(int key1, int key2);
        void DFTraversal(vertex *n);
        void setAllVerticesUnvisited();
        void printGraph();


    private:
        std::vector<vertex*> vertices; // vector of vertex pointers called verticies

};

#endif
