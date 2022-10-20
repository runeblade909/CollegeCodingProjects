#include "Part2_Graph.hpp"
#include <vector>
#include <queue>
#include <iostream>

using namespace std;

// function to add edge between two vertices
void Graph::addEdge(int v1, int v2){

    for(int i = 0; i < vertices.size(); i++){
        if(vertices[i]->key == v1){
            for(int j = 0; j < vertices.size(); j++){
                if(vertices[j]->key == v2 && i != j){
                    adjVertex av;
                    av.v = vertices[j];
                    vertices[i]->adj.push_back(av);
                    //another vertex for edge in other direction
                    adjVertex av2;
                    av2.v = vertices[i];
                    vertices[j]->adj.push_back(av2);
                }
            }
        }
    }
}


// function to add a vertex to the graph
void Graph::addVertex(int n){
    bool found = false;
    for(int i = 0; i < vertices.size(); i++){
        if(vertices[i]->key == n){
            found = true;
        }
    }
    if(found == false){
        vertex * v = new vertex;
        v->key = n;
        vertices.push_back(v);
    }
}


// GOLD : function to print the path of the shortest path from src to dest vertex
// First modify the findShortestPath function to update the predecessor(pred) when updating the distance of the adjacent vertices,
// then complete the following function

void Graph::printPath(int src, int dest)
{
    // path vector stores the vertices of the shortest path from the vertices vector
    vector<vertex*> path;
    vertex *node;
        for (int i = 0; i < vertices.size(); i++) 
        {
            if (vertices[i]->key == dest) 
            {
                node = vertices[i];
            }
        }
        while (true) 
        {   
         // insert(place to insert, data to insert)
         path.insert(path.begin(), node);
            if (node->key == src) 
            {
                break;
            }
         node = node->pred;
        }
    for (int i = 0; i < path.size(); i++) 
    {
        cout << path[i]->key << " ";
    }
}

// SILVER: Complete the following function and return the length of the shortest path (return -1 if you can't find a path)
int Graph::findShortestPath(int src, int dest){
    // Need to find the source vertex since we only have it's key 'src'
    // A pointer for head vertex
    vertex *head;
    for(int i = 0; i < vertices.size(); i++)
    {
        if(vertices[i]->key == src)
        {
            head = vertices[i];
        }
    }
    head->visited = true;
    // Use the queue to keep track of visited vertices
    queue<vertex*> q;
    // Enqueue the head
    q.push(head);
    vertex *n;

    while(!q.empty())
    {
        n = q.front();
        q.pop();
        if(n->key == dest)
        {
        return n->distance;
        }
        // go to all the adjacent vertices of n
        for( int x = 0; x < n->adj.size(); x++ )
        {
            vertex *r = (n->adj[x]).v;
            if(r->visited == false)
            {
                r->visited = true;
                r->distance = (n->distance + 1);
                r->pred = n;
                q.push(r);
            }
        }
    }
    return -1;
}

