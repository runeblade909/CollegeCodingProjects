#include <iostream>
#include "bst.hpp"

using namespace std;


BST::BST() {
	root = NULL;
}

BST::~BST() {
	deleteTree(root);
	root = NULL;
}
void BST::deleteTree(Node* node) {
	if (node != NULL) {
		deleteTree(node->leftChild);
		deleteTree(node->rightChild);
		delete node;
	}
}

/*
** Implement the following function to return the count of comparisons, 
**   you may implement helper functions.
*/
int BST::searchCounter(int target) {
	// your code here!

	Node *next = root;
	int count = 1;
	//look through the bst until the next node in line is NULL
	//We are just going to use a basic, greater than or less than code analysis to select which child is needed to be visited

	while(next != NULL)
	{

		//If not already there, count
		if(target != next ->key)
		{
			count = count+1;
		}

		//Right Tree
		if (target > next->key && next->rightChild == NULL)
		{
			count = count - 1;
			return count;

		}
		else if(target > next->key)
		{
			next = next->rightChild;
		}

		// Left Tree
		if (target < next->key && next->leftChild ==NULL)
		{
			count = count - 1;
			return count;

		}
		else if(target < next->key)
		{
			next = next->leftChild;
		}

		// Double NULL
		if(next->leftChild ==NULL && next->rightChild ==NULL)
		{
			count = count +1;
				if(next->key == target)
				{
					count = count - 1;
				}


			return count;
		}

		// Return when found
		if(target == next->key)
		{			
			return count;
		}

	}

	return count;

}
