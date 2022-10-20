#include <iostream>
#include<iomanip>
#include<string>
#include<fstream>
#include<sstream>

using namespace std;

struct wordRecord {
   string word;
   int count;
};

//Function prototypes
void getIgnoreWords(const char* ignoreWordFileName, string ignoreWords[]);

bool isIgnoreWord(string word, string ignoreWords[]);

int getTotalNumberNonIgnoreWords(wordRecord distinctWords[], int length);

void sortArray(wordRecord distinctWords[], int length);

void printTenFromN(wordRecord distinctWords[], int N, int totalNumWords);



int main(int argc,char* argv[])
{
   //Command line error check
   if (argc != 4) {
       cout << "Usage: Assignment2Solution <number of words> <inputfilename.txt> <ignorewordsfilename.txt>" << endl;
       exit(0);
   }

   //Read values from command line
   int index = stoi(argv[1]);

   ifstream in(argv[2]);

   //Array for words
   wordRecord* records;

   string ignoreWords[50];

   records = new wordRecord[100];   //Pointer to dynamic vector

   int size = 100; // Starting Size

   int doublingCnt = 0; // Double Count

   //Read ignore words into array
   getIgnoreWords("ignoreWords.txt", ignoreWords);

   //Error check in open file of words file
   if (!in) {
       cout << "Failed to open " << argv[2] << endl;
       exit(0);
   }

   string line; // for getline and ss
    string word;

   int i = 0;

   //Read until the end of file line by line

    // get rid of getline and do word<< in

    // Need to break up the lines by spaces and paragraph breaks "\n"


   while (!in.eof()) {
       getline(in, line);

       stringstream ss(line);

       string word;

       //Split line into words
       while (getline(ss, word, ' ')) {
           //Check each word present in ignore array
           if (isIgnoreWord(word, ignoreWords) == false) {

               if (i < size) {
                   // Preset check
                   bool check = false;

                    // Check the arrays
                   for (int j = 0; j < i; j++) {

                       if (records[j].word == word) {

                           records[j].count += 1;

                           check = true;

                           break;
                       }
                   }
                   if (!check) {

                       records[i].word = word;

                       records[i].count = 1;

                       i++;
                   }
               }

               else {
                   bool check = false;

                   for (int j = 0; j < i; j++) {

                       if (records[j].word == word) {

                           records[j].count += 1;

                           check = true;
                           
                           break;
                       }
                   }

                    // Array Doubling
                   if (!check) {

                       doublingCnt++;

                       wordRecord *temp;

                       temp = new wordRecord[size];

                       for (int k = 0; k < size; k++) {

                           temp[k] = records[k];
                       }

                       records = new wordRecord[size * 2];

                       for (int k = 0; k < size; k++) {

                           records[k] = temp[k];
                       }

                       records[size].word = word;

                       records[size].count = 1;

                       size = size * 2;

                       i++;
                   }
               }
           }
       }
   }
   in.close();

   sortArray(records, i);



  //  cout << "Size of records:" << sizeof(records) << endl;

    //This gets rid of the beginning mess up thing

    if(records[0].word.empty()==1){


                       wordRecord *temp;

                       temp = new wordRecord[size];

                       for (int k = 0; k < size; k++) {

                           temp[k] = records[k];
                       }

                       records = new wordRecord[size];

                       for (int k = 0; k < size; k++) {

                           records[k] = temp[k+1];
                       }
                       i=i-1;

    }
    
    /*


    cout << "Size of records:" << sizeof(records) << endl;

   //bool t = records[0].word.empty();

    //  cout << "THIS IS" << t << endl;




   cout << "This is some word:" << records[0].word << "Yeah" << endl;
   
    */
   cout << "Array doubled: " << doublingCnt << endl;

   cout << "Distinct non-common words: " << i << endl;

   cout << "Total non-common words: " << getTotalNumberNonIgnoreWords(records, i) << endl;

   cout << "Probability of next 10 words from rank " << index << endl;

   cout << "---------------------------------------\n";

   printTenFromN(records, index, getTotalNumberNonIgnoreWords(records, i));

   return 0;
}




//Function read ignore words from file
void getIgnoreWords(const char* ignoreWordFileName, string ignoreWords[]) {
   int i = 0;

   ifstream in(ignoreWordFileName);

   if (!in) {
       cout << "Failed to open "<<ignoreWordFileName << endl;
       exit(0);
   }

   while (in >> ignoreWords[i]) {
       i++;
   }

   in.close();
}

//Function return true if the word found in array else false
bool isIgnoreWord(string word, string ignoreWords[]) {

   for (int i = 0; i < 50; i++) {

       if (ignoreWords[i] == word) {
           return true;
       }
   }
   
   return false;
}

//Return the sum of all words
int getTotalNumberNonIgnoreWords(wordRecord distinctWords[], int length) {
   int sum = 0;

   for (int i = 0; i < length; i++) {
       sum += distinctWords[i].count;
   }

   return sum;
}

//Sort array according to frequency
void sortArray(wordRecord distinctWords[], int length) {

   for (int i = 0; i < length; i++) {

       for (int j = i + 1; j < length; j++) {

           if (distinctWords[i].count < distinctWords[j].count) {
               wordRecord wr = distinctWords[i];

               distinctWords[i] = distinctWords[j];

               distinctWords[j] = wr;
           }
       }
   }
}

//Display 10 words their probability after given index
void printTenFromN(wordRecord distinctWords[], int N, int totalNumWords) {

   for (int i = N; i < N + 10; i++) {
       float prob = (float)distinctWords[i].count / totalNumWords;
       
       cout << fixed << setprecision(4) << prob << " - " << distinctWords[i].word << endl;
   }
}