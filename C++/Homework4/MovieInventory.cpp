#include "MovieInventory.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

MovieInventory::MovieInventory() {
    root = NULL;

}

MovieInventory::~MovieInventory() {
  delete root;
  
}

void inOrderPrint(MovieItem*node){

    if(node == NULL){
    return;
    }
    inOrderPrint(node->left);
    cout << "Movie: " << node->title << " " << node->rating << endl;    //Prints movies in alphabetical order
    inOrderPrint(node->right);
}
void MovieInventory::printMovieInventory()
{

    if(root == NULL){
        cout << "Tree is Empty. Cannot print" << endl;    //no movies
    }
    else{
    MovieItem*node = root;          //calls function directly above
    inOrderPrint(node);
    }
}

void MovieInventory::addMovieItem(int rating, string title, int year, float ranking) {

  MovieItem *temp = new MovieItem(0, "", 0, 0);
  MovieItem *compare = new MovieItem(0, "", 0, 0);
  temp->ranking = rating;     //sets all stats for the inputed movie 
  temp->title = title;
  temp-> year = year;
  temp->rating = ranking;
  temp -> left = NULL;
  temp -> right = NULL;
  if(root == NULL){
    root = temp;        //if tree is empty
  }
  else{
    compare = root;
    while(compare-> left != NULL || compare -> right != NULL){
      if(title.compare(compare->title) < 0 && compare -> left != NULL){       //if compare is higher in alphabet move left
      compare = compare -> left;

      }
      if(title.compare(compare->title) > 0 && compare -> right != NULL){   //if compare is lower in alphabet move right
        compare = compare -> right;

      }

      
      if(title.compare(compare->title) > 0 && compare -> right == NULL){      //found open slot and its lower in alphabet
        compare -> right = temp;
        break;

      }
      if(title.compare(compare->title) < 0 && compare -> left == NULL){       //open slot higher in alphabet
      compare -> left = temp;
      break;

      
      }

    }
  
    if(compare ->left == NULL && compare -> right == NULL){     //if there is only a root
        if(title.compare(compare->title) < 0){
            compare -> left = temp;
        }
        else{
            compare->right = temp;
        }
    }
  }

}

void postOrderSearch(MovieItem*&temp, string &s, bool &found)
{
  if(temp == NULL){

    return;
  }
  

  postOrderSearch(temp->left, s, found);
  postOrderSearch(temp->right, s, found);         //traverses through and prints movie stats below when found
  if(temp->title == s){

    cout << "Movie Info:" << endl;
    cout << "==================" << endl;
    cout << "Ranking:" << temp->ranking << endl; 
    cout << "Title :" << temp->title << endl; 
    cout << "Year :" << temp->year << endl; 
    cout << "rating :" << temp->rating << endl;
        found = true;
        return;
    }


}


void MovieInventory::getMovie(string title)
{
    // Your code here
    MovieItem *temp = root;
    string s = title;
    bool found = false;
    postOrderSearch(temp, s, found);        //call function above
    if(found == false){
      cout<<"Movie not found."<<endl;       //if the input title doesn't match any titles in the tree
    }
        
    
}


void inOrderSearch(MovieItem*node, float rating, int year){

    if(node == NULL){
    return;
    }
    inOrderSearch(node->left, rating, year);                    //recursive in order
    if(node->year >= year && node->rating > rating){                            //if the year is after or equal and rating is better
    cout << node->title << "(" << node->year << ") " << node->rating << endl;   //this is printed
    }
    inOrderSearch(node->right, rating, year);
}

void MovieInventory::searchMovies(float rating, int year) {
  //write your code
  if(root != NULL){
    cout << "Movies that came out after " << year << " with rating at least " << rating << ":" << endl;   //header
    MovieItem *temp = root;
    inOrderSearch(temp, rating, year);      //function directly above

  

  }
  else{                                                       //no movies in tree
    cout << "Tree is Empty. Cannot query Movies" << endl;
  }



}

void countAndTotal(MovieItem *node,double &runningtotal,double &totalcount){    //adds up all rating and counts
  if(node == NULL){
    return;
  }
  else{
  countAndTotal(node->left, runningtotal, totalcount);      //recursive
  countAndTotal(node->right, runningtotal, totalcount);
  runningtotal += node-> rating;                        //adds running total of ratings
  totalcount++;                                         //adds total count
  }
  
}

void MovieInventory::averageRating() {
  //write your code
  MovieItem *temp = root;
  double runningtotal = 0.0;
  double totalcount = 0.0;
  double average = 0.0;
  countAndTotal(temp, runningtotal, totalcount);
    average = runningtotal / totalcount;            //running total of ratings is then divided by total count to get the average
    if(temp == NULL){
       cout<<"Average rating:0.0"<<endl;          //no movies in  tree
    }
    else{
        cout<<"Average rating:"<<average<<endl;
  }
}


void postOrderFind(MovieItem*&temp, string &s, bool &found)
{
while(temp->title.compare(s) != 0){
  if(temp->title.compare(s) > 0){     //compares the alphabetical until the exact word is found
    temp = temp->left;
  }
  else{
    temp = temp->right;
  }

}
found = false;
return;
}


void MovieInventory::deleteMovieItem(string title)
{
  //write your code
  MovieItem *temp = new MovieItem(0, "", 0, 0);
  temp = root;
  bool found = false;
  MovieItem *holdingnode = new MovieItem(0, "", 0, 0);
  postOrderFind(temp, title, found);                //call function above

  
if(temp->left == NULL && temp->right == NULL)       //temp is the root
    {
        delete temp;
        temp = NULL;
    }
    //Only right child
    else if(temp->left == NULL)
    {
        holdingnode = temp->right;          //search for smallest on right side to replace
        while(holdingnode->left != NULL){
          holdingnode = holdingnode->left;
        }
        delete temp;
        temp = holdingnode;
        
    }
    //Only left child
    else if(temp->right == NULL)
    {
        delete temp;
        temp = temp->left;
        temp->left = NULL;
    }
    //Both left and right child
    else
    {
      ///Replace with Minimum from right subtree
    holdingnode = temp->right;
    while(holdingnode->left != NULL){
        holdingnode = holdingnode->left;
    }
    delete temp;
    temp->title = holdingnode->title;   //transfer info
    temp->ranking = holdingnode->ranking;
    temp->year = holdingnode->year;
    temp->rating = holdingnode->rating;

    }






}



void MovieInventory::leftRotate(string title)
{
    MovieItem*X = root;
    MovieItem*Y = new MovieItem(0, "", 0, 0);
    bool xroot = false;
    bool rotate = false;
    bool found = false;
    postOrderFind(X, title, found);

    Y = X->right;
    if(root == X){
      xroot = true;
    }
    if(X->right != NULL){
    X->right = Y->left;     //left subtree of Y becomes right of X
    Y->left = X;            //left side of Y becomes x
    rotate = true;
    }

    if(xroot == true && rotate == true){
      root = Y;
    }

}