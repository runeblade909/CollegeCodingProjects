#include "Graph.hpp"
#include <vector>
using namespace std;


/*
 * adds a vertex to the graph
 */
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

/*
 * adds an edge between two vertices (directed graph)
 */
void Graph::addEdge(int src, int dest){
    for(int i = 0; i < vertices.size(); i++) {
        if(vertices[i]->key == src) {
            for(int j = 0; j < vertices.size(); j++) {
                if(vertices[j]->key == dest && i != j) {
                    adjVertex av;
                    av.v = vertices[j];
                    vertices[i]->adj.push_back(av);
                }
            }
        }
    }
}

/*
 * Complete the following function which checks if vert is a boss vertex
 */
bool Graph::isVertexABoss(vertex *vert) 
{

        int visitedNode = 0;
      vert->visited = true;


    // Use DFS to look and see if you can traverse the whole table in one path
    for(int x = 0; x < vert->adj.size(); x++ )
    {
        if (!vert->adj[x].v->visited)
        {
            isVertexABoss(vert->adj[x].v);
        }
    }
        // Count to make sure all nodes are visited
        for(int i = 0; i < vertices.size() ; i++)
        {
            if(vertices[i]->visited == true)
            {
                visitedNode = visitedNode + 1;

            }

        }
        //If all nodes visited, that vertex is a boss vertex.
        if(visitedNode == vertices.size())
        {

            return true;

        }
        //Or its not.......
        else
        return false;

}

/*
 * Complete the following function which checks if the graph is a Boss
 */
bool Graph::isGraphABoss() 
{
    // If isVertexABoss returns true ever, then it is a boss graph.
        int counter = 0;

    for(int i = 0; i < vertices.size(); i++)
    {
        if(isVertexABoss(vertices[i]) == true)
        {
            counter = counter +1;
        
        }
    }
    
    if( counter > 0 ){
    return true;
    }
    else
    return false;
    
}
