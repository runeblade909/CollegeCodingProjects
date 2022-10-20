#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <bits/stdc++.h>

using namespace std;
char letterArr[] = {32,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122};
int const length = 53;
string keyPhrase = "May your mountains rise";
char keyPhrase1[23] = {};
float const keyPhraseLength = 23;
char DNA[23] = {};
char DNA2[23] = {};
char childDNA[23] = {};
string populationArray[200] = {};
int const populationSize = 200;
int fitness[23] = {};
int fitness1[23] = {};
int Scores = {};
float ScoresArray[200] ={};
float ScoresPercentArray[200] = {};
float FitnessScores[200] = {};
float max_ele;
float matingFactor = 10;
float indiciesCounter = 0;
float indiciesArray[200] = {};
float matingPool[2000] = {};
int matingPoolSize = 2000;
int ticketCounter = 0;
int p1 = 0;
int p2 = 0;
float randParentArr[200] = {};
int midpoint = 0;
int mutationRate = 4;

void buildPopulation();
void calculateFitness();
void buildMatingPool();
void practicuum (float arr[], int size, float max_ele);
void breed();
void causeMutation();

int main()
{
    //Seed timer

    srand(time(NULL)); 

    // show our desired phrase

    cout << endl << "We are looking to generate the phrase: " << keyPhrase << endl;
    
    //function to build our population of 200 people with 23 char long arrays

    buildPopulation();

    //Fitness Calculation

    calculateFitness();

    //Build Mating build mating pool

    buildMatingPool();
    practicuum (FitnessScores, populationSize, max_ele);


    //Now time for the parent strands to make a kid
    breed();

    //Now to make that kid just a little different
    causeMutation();


    return 0;
}

void buildPopulation()
{
    for (int j=0; j<= 199; j++)
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
           

    //Calculate Fitness
    for (int k = 0; k <= keyPhraseLength-1; k++)
    {
     fitness[k] = DNA[k] == keyPhrase1[k];

     Scores += fitness[k];

    }
    
 // spit out number of matching values and also the percent
    ScoresArray[y] = Scores;
    Scores = 0;
    ScoresPercentArray[y] = ScoresArray[y]/keyPhraseLength;

 }
}

void buildMatingPool() 
//Puts everything that was calculated into arrays for orginization
{
    //loop for all paprts of the array
    for (int y = 0; y <= populationSize-1; y++)
    {
        // rounds to second decimal place
        FitnessScores[y] = (ceil(ScoresPercentArray[y]*100)/100);

        // Helps in the beginning
        if (FitnessScores[y]<=.099 && FitnessScores[y] > 0)
        {
           FitnessScores[y] = .1;
        }

        // Create Number of tickets
        FitnessScores[y] = round(FitnessScores[y]*10);

        
    }
}

void practicuum (float arr[], int size, float max_ele) 
// Second part of building the mating pool, index matrix along with the tickets
{    
    max_ele = 0;
 // Find max element
    for (int i =0; i<=size;i++)
    {

        if (max_ele <= arr[i])
        {
            max_ele = arr[i];

        }

        
    }
 //Give people their tickets
    for (int i = 0; i<=size; i++)
    {
        arr[i] = arr[i]/max_ele;
        arr[i] = arr[i] * matingFactor;
    }

 //Create an indicies array
    for (int i = 0; i<= size; i++)
    {
        if(arr[i] > 0)
        {
        indiciesArray[i] = indiciesCounter;
        }
        else
        {
            indiciesArray[i] = {};
        }
        indiciesCounter += 1;        
    }
 
 // TRYING to create the mating pool (RIGHT NOW ITS ONLY DISPLAYING WHAT I WANT TO SEE, NOT GIVING ME AN ACTUAL ARRAY TO WORK with)
    
    ticketCounter = 0;
    
    for(int i = 0; i <= size; i++)
    {
        for(int k = 0; k <= arr[i]-1; k++)
        {
            if (indiciesArray[i] > 0)
            {
                matingPool[ticketCounter] = indiciesArray[i];
                ticketCounter++;
            }
        }
    } 
}

void breed()
{
    cout << endl;
    
 // Select a random person from the indiciesArray
    p1 = rand() % ticketCounter;
    
 // Now for p2
    p2 = rand() % ticketCounter;

    cout <<"Our ticket number is: " << p1 << " Followed by the ever so lucky: " << p2 << endl;
 //Thank CHRIST it works
    p1 = matingPool[p1];
    p2 = matingPool[p2];
    cout << keyPhrase << endl << populationArray[p1] << endl << populationArray[p2] << endl << endl;

 // MIDPOINT HERE, we want to make a random number from 0-23 to choose as our midpoint, parent 1 gets the left side of the midpoint
 // and parent 2 gets the right side

    midpoint = 1 + rand() % 21;


 //Put element of population array (string) back into a char array
    strcpy(DNA , populationArray[p1].c_str());
 //Then for p2
    strcpy(DNA2 , populationArray[p2].c_str());

   for (int i = 0; i<midpoint;i++)
   {
       childDNA[i] = DNA[i];

   }

   for (int k = midpoint; k<=keyPhraseLength ; k++)
    {
        childDNA[k] = DNA2[k];
    }

   for (int i = 0; i<=keyPhraseLength;i++)
   {
       cout << childDNA[i];
   }
 cout << endl;

     // RANDOM BREEDING HERE
     //Im going to select random indicies 
    //     int randDNA[23] = {};
    //     for (int i = 0; i<= keyPhraseLength; i++)
    //     {
    //      randDNA[i] = rand() %2;
    //     }
    //     for (int j = 0; j<= keyPhraseLength; j++)
    //     {
    //         if(randDNA[j]>0)
    //         {
    //             childDNA[j] = DNA[j];
    //         }
    //         else
    //         {
    //             childDNA[j] = DNA2[j];
    //         }
            
    //         cout << childDNA[j];
    //     }
    //  cout << endl;
}   


void causeMutation()
{
    // We are going to make a random number pool from 0-100 where only 4 spots are more than 0, therefore giving a 4 percent chance that those 4 numbers are picked
 int randNumChance[100] = {0};
 int counter = 0;
 for (int i = 0; i <= 99 ; i++)
 {
    randNumChance[i] = rand() %2;

    if (randNumChance[i] > 0)
    {
        counter++;
    }

    if (counter == mutationRate)
    {
     break;
    }

 }
 int randChoice = 0;
    for (int i = 0; i<= keyPhraseLength; i++)
    {
        randChoice = rand() % 100;
        if (randNumChance[randChoice] > 0)
        {
            childDNA[rand()%23] = letterArr[rand()%53];
        }

        cout << childDNA[i];
    }
    cout << endl;
}





