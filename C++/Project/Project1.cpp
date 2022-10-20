#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;
int const keyPhraseLength = 18;
int const populationSize = 200;
char letterArr[] = {32,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122};
int const length = 53;
string keyPhrase = "To be or not to be";
string Child = ""; 
char keyPhrase1[keyPhraseLength] = {}; 
char DNA[keyPhraseLength] = {}; 
char DNA2[keyPhraseLength] = {};
char childDNA[keyPhraseLength] = {};
string populationArray[populationSize] = {}; 
string TempPopulationArray[populationSize] = {};
int fitness[keyPhraseLength] = {}; 
int fitness1[keyPhraseLength] = {};
int Scores = {}; 
float ScoresArray[populationSize] ={};
float ScoresPercentArray[populationSize] = {};
float FitnessScores[populationSize] = {};
int FitnessTickets[populationSize] = {};
float max_ele;
const float matingFactor = 10;
float indiciesCounter = 0;
float indiciesArray[populationSize] = {};
int const matingPoolSize = populationSize*matingFactor;
float matingPool[matingPoolSize] = {};
int ticketCounter = 0;
int p1 = 0; 
int p2 = 0; 
float randParentArr[populationSize] = {};
int midpoint = 0; 
int mutationRate = 4; 
int generations = 1;
int numOfMatches = 0;
int THESPOT = 0;


void buildPopulation();
void calculateFitness();
void buildMatingPool();
void practicuum (float max_ele);
void breed();
void causeMutation();
int randDNA[keyPhraseLength-1] = {0};

int main()
{
    //Trying to do the whole graph thing
    ofstream outputStream;
    ofstream outputYEET;

    outputYEET.open("fitness.csv");

    outputStream.open("max_ele.csv");


    
    //Seed timer

    srand(time(NULL)); 

    // show our desired phrase

    cout << endl << "We are looking to generate the phrase: " << keyPhrase << endl;
    
    //function to build our population of 200 people with 23 char long arrays

    buildPopulation();

    //Fitness Calculation
while ( numOfMatches != 1)
 {
        
    //Fitness Calculation
    calculateFitness();

    //Build Mating build mating pool

    buildMatingPool();
    
    practicuum (max_ele);

    // for (int i = 0; i<=populationSize-1;i++)
    // {
    //     cout << populationArray[i] << " ";
    // }
    // cout << endl;


     
    for ( int j = 0; j<=populationSize-1; j++)
    {
        //Now time for the parent strands to make a kid
        breed();
        //cout << ticketCounter << " ZOOWEEMAMA" << " " << p1<< " " << p2 << " " << Child << " " << endl;
        //Now to make that kid just a little different
        causeMutation();
        TempPopulationArray[j] = Child;
        
    }
    // cout << endl;
    //Make a temp array to store the new children better and swap them
    for (int k = 0; k<= populationSize-1;k++)
    {
        populationArray[k] = TempPopulationArray[k];
    }
     p1 = 0;
     p2 = 0;
    // Just make sure we are working with all of the right parts
    calculateFitness();
    
    // Now, we will show the fittest member of the population.
    for (int i =0; i<=populationSize-1;i++)
    {
        if (max_ele <= ScoresPercentArray[i])
        {
            max_ele = ScoresPercentArray[i];
            THESPOT = i;

        }
        outputYEET << ScoresPercentArray[i] << ",";

    }
    outputYEET << "\n";
    cout << endl;
    numOfMatches = 0;
    //This is what will stop the loop, it scans for the finished string 
    for (int k = 0; k<=populationSize -1; k++)
    {
        if (ScoresPercentArray[k] >= .99)
        {
           numOfMatches = numOfMatches+1;
        }


    }


    // Here is the final output per generation
     cout << endl << "Generation: " << generations << " The fittest member is: " << populationArray[THESPOT] << " Fitness: " << max_ele << endl;
     generations++;
    //Grabing them most fit peeps
    outputStream << max_ele;
    outputStream << "\n";
 } // END OF WHILE LOOP
    outputYEET.close();
    outputStream.close();
    
    return 0;
}

void buildPopulation()
{
    
    for (int j=0; j<= populationSize-1; j++)
    {
        for (int i=0; i<=keyPhraseLength-1; i++)
        {  
            DNA[i] = letterArr[rand()%53];
        }
        populationArray[j] = DNA;
     int count;
     count = j +1;
    }
}

void calculateFitness()
{
  for (int y = 0; y <= populationSize-1; y++)
 {
     //Put element of population array (string) back into a char array
    strcpy(DNA , populationArray[y].c_str());

    //Put keyphrase into a char array    
    strcpy(keyPhrase1, keyPhrase.c_str());
           
    Scores = 0;
    //Calculate Fitness
    for (int k = 0; k <= keyPhraseLength-1; k++)
    {
     fitness[k] = DNA[k] == keyPhrase1[k];

     Scores += fitness[k];

     //cout << fitness[k];

    }

    //cout << " " << Scores << " ";
 // spit out number of matching values and also the percent
    ScoresArray[y] = Scores;
    Scores = 0;
    ScoresPercentArray[y] = ScoresArray[y]/keyPhraseLength;
    //cout << ScoresPercentArray[y] << " ";

 }
}

void buildMatingPool() 
//Puts everything that was calculated into arrays for orginization
{
    //loop for all parts of the array
    for (int y = 0; y <= populationSize-1; y++)
    {
        // rounds to second decimal place
        FitnessScores[y] = (round(ScoresPercentArray[y]*100)/100);

        // Helps in the beginning
        if (FitnessScores[y]<=.099 && FitnessScores[y] > 0)
        {
           FitnessScores[y] = .1;
        }

        // Create Number of tickets
        FitnessScores[y] = round(FitnessScores[y]*10);
        //cout << FitnessScores[y] << " ";

    }
    //cout << endl;
}

void practicuum (float max_ele) // taking in fitness tick, pop size and max_ele
// Second part of building the mating pool, index matrix along with the tickets
{    
    max_ele = 0;
 // Find max element
    for (int i =0; i<=populationSize-1;i++)
    {
        if (max_ele <= FitnessScores[i])
        {
            max_ele = FitnessScores[i];

        }   
    }
    //cout << endl << max_ele << endl;
 //Give people their tickets
        FitnessTickets[populationSize] = {0};
    for (int i = 0; i<=populationSize-1; i++)
    {
    
        FitnessScores[i] = FitnessScores[i]/max_ele;
        FitnessTickets[i] = FitnessScores[i] * matingFactor;
        //cout << FitnessTickets[i] << " ";
    }
    

    //cout << endl;
 //Create an indicies array
    indiciesArray[populationSize] = {0};
    indiciesCounter = 0;
    for (int i = 0; i<= populationSize-1; i++)
    {
        if(FitnessTickets[i] > 0)
        {
        indiciesArray[i] = indiciesCounter;
        }
        else
        {
            indiciesArray[i] = {};
        }
        indiciesCounter += 1;    
        //cout << indiciesArray[i] << " ";    
    }
    //cout << endl;
 // TRYING to create the mating pool (RIGHT NOW ITS ONLY DISPLAYING WHAT I WANT TO SEE, NOT GIVING ME AN ACTUAL ARRAY TO WORK with)
    matingPool[2000] = {0};
    ticketCounter = 0;
    
    for(int i = 0; i <= populationSize-1; i++)
    {
        for(int k = 0; k <= FitnessTickets[i] -1; k++)
        {
            if (indiciesArray[i] > 0)
            {
                matingPool[ticketCounter] = indiciesArray[i];
                ticketCounter++;
            }

        }
    } 
    // for (int k = 0; k<= 1999 ; k++)
    // {
    //     cout << matingPool[k];
    // }
    // cout << endl << ticketCounter << endl;
}

void breed()
{
    
 // Select a random person from the indiciesArray
    p1 = rand() % ticketCounter;
    
    
 // Now for p2
    p2 = rand() % ticketCounter;
    

 //Thank CHRIST it works
    p1 = matingPool[p1];
    p2 = matingPool[p2];

     //Put element of population array (string) back into a char array
    strcpy(DNA , populationArray[p1].c_str());
     //Then for p2
    strcpy(DNA2 , populationArray[p2].c_str());

 // MIDPOINT HERE, we want to make a random number from 0-23 to choose as our midpoint, parent 1 gets the left side of the midpoint
 // and parent 2 gets the right side

    midpoint = 1 + rand() % keyPhraseLength;

   for (int i = 0; i<midpoint-1;i++)
   {
       childDNA[i] = DNA[i];

   }

   for (int k = midpoint; k<=keyPhraseLength -1  ; k++)
    {
        childDNA[k] = DNA2[k];
    }


    //  // RANDOM BREEDING HERE
    //  //Im going to select random indicies 
    //     for (int i = 0; i<= keyPhraseLength-1; i++)
    //     {
    //      randDNA[i] = rand() %2;
    //      //cout << randDNA[i] <<" ";
    //     }
    //     //cout << endl;
        
    //     //for (int i = 0; i<= keyPhraseLength-1; i++)
    //     //{
         
    //      //cout << DNA[i] <<" ";
    //     //}
    //     //cout << endl;
    //     //for (int i = 0; i<= keyPhraseLength-1; i++)
    //     //{
    //      //cout << DNA2[i] <<"T";
    //     //}

    //     for (int j = 0; j<= keyPhraseLength-1; j++)
    //     {
    //         if(randDNA[j]=1)
    //         {
    //             childDNA[j] = DNA[j];
    //         }
    //         else if (randDNA[j]=0)
    //         {
    //             childDNA[j] = DNA2[j];
    //         }
    //         //cout << childDNA[j];
            
    //     }
        //cout << endl;
}   


void causeMutation()
{
    // We are going to make a random number pool from 0-100 where only 4 spots are more than 0, therefore giving a 4 percent chance that those 4 numbers are picked
 
 int randChoice = 0;
    for (int i = 0; i<= keyPhraseLength -1; i++)
    {
        randChoice = rand() % 101;
        if (randChoice <= mutationRate)
        {
            childDNA[i] = letterArr[rand()%53];
            break;
        }


    }

    Child = childDNA;
}