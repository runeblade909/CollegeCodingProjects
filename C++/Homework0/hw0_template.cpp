/*
** Name: Zak Reichenbach
** Assighment: Week 1 - HW0
** Description: Student Grades Array of Structures
**
** Command-line example:
** hw0 students.csv output.csv C A
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

struct studentData {
   string studentName;
   int homework;
   int recitation;
   int quiz;
   int exam;
   double average;
};

void addStudentData(studentData students[], string studentName, int homework, int
recitation, int quiz, int exam, int length, int i)
{
   // insantiate studentData struct
   // calculate average

      // I am not sure why, but this is the only way I could figure out how to actually output the structs needed.

      studentData weh;
   
      students[i].studentName = studentName;
      students[i].homework = homework;
      students[i].recitation = recitation;
      students[i].quiz = quiz;
      students[i].exam = exam;
      students[i].average = (weh.homework+weh.recitation+weh.quiz+weh.exam)/4;
         
         
         /*
         cout << weh.studentName << endl;
         cout << weh.homework << endl;
         cout << weh.recitation << endl;
         cout << weh.quiz << endl;
         cout << weh.exam << endl;
         */


      students[i] = weh;

      i++;
      // cout << i << endl;
     
   // populate all struct variables
   // add struct to the students array
}

char calcLetter(double avg)
{
   // calculate grade letter
   // >90 = A
   // 80-89.9 = B
   // 70-79.9 = C
   // 60-69.9 = D
   // <60 = F
   char grade;
   

   if (avg >= 90)
   {
      grade = 'A';
   }

   else if (avg >= 80 && avg < 89.9)
   {
      grade = 'B';
   }
   else if (avg >= 70 && avg < 79.9)
   {
      grade = 'C';
   }
   else if (avg >= 60 && avg < 69.9)
   {
      grade = 'D';
   }
   else if (avg < 60)
   {
      grade = 'F';
   }


   return grade;
   // return grade
}

void printList(const studentData students[], int length,char grade[] , int i)
{
   // Loops through vector and calls calcLetter to give people their grade.
   while(i<length)
   {
    grade[i] = calcLetter(students[i].average);
    cout << students[i].studentName << " earned a " << students[i].average << "%. Which is a " << grade[i] <<"." << endl;


    i++;
   }
   // loop through students and display output
}



int main(int argc, char* argv[])
{
   // declare local variables
   
   ifstream in_file("students.csv");

      if(!in_file.is_open())
      {

         cout <<"It didn't open it broken." << endl;
      }
 
   // read command-line arguments

   // Grade Bounds
   char UpperBound, LowerBound;




   UpperBound = argv[1][0];
   LowerBound = argv[2][0];
   



   // Change length to the input arguements for the given csv file useing the argv and argc stuff.
 
   studentData students[10];
 

 string line;
 string name;
 int homework;
 int recitation;
 int quiz;
 int exam;
 int length = 10;

     



     
   //getline(in_file, line,",");

   int i = 0;
   int j = 0;
   studentData TIM;
   while(getline(in_file, line) && i < length)
   {


      // use getline() to parse students.csv
      stringstream ss(line);
      getline(ss,name,',');
      getline(ss,line,',');
      homework = stoi(line);
      getline(ss,line,',');
      recitation = stoi(line);
      getline(ss,line,',');
      quiz = stoi(line);
      getline(ss,line,',');
      exam = stoi(line);
      length = 10;
   
      /*
      cout << name << endl;
      cout << homework << endl;
      cout << recitation << endl;
      cout << quiz << endl;
      cout << exam << endl;
      */


     

      addStudentData(students, name, homework, recitation, quiz, exam, length,0);
     
         
     
         students[j] = students[i];    //THIS LINE DOES NOT GET COMMENTED OUT OK YOU SMALL BIG BRAINED YOU
         
           // cout << students[j].studentName << endl;

         j++;


         // The first variable is constantly overwritten, so this catches the first person on the list, who is named Tim Thomas
         if (students->studentName == "Tim Thomas")
         {
            TIM = students[0];
         }



        // cout << TIM.studentName<< endl;
   }
      // This puts Tim back on top, woooooooo
      students[0] = TIM;


      // Sanity Checker
      /*
      for(int k = 0; k <10; k++)
      {
         cout << students[k].studentName << endl;
         cout << students[k].homework << endl;
         cout << students[k].recitation << endl;
         cout << students[k].quiz << endl;
         cout << students[k].exam << endl;
         cout << students[k].average << endl;
      }
      */
     

     char Grades[10];

     

     printList(students,10,Grades,0);

     


   // write to output.csv based on lower bound and upper bound
   // which are retrieved from command-line args
   
   //= printList();

 // now to output our grade range we want.

     ofstream out_file;
     out_file.open("Grades.csv");



      // char array with indexing for user input grades
   



   for(int i = 0; i< length; i++)
   {
      if(Grades[i] == 'A')
      {
      out_file<<students[i].studentName<<","<<students[i].average<<","<<Grades[i]<<endl;
      }

      else if(Grades[i] == 'B')
      {
      out_file<<students[i].studentName<<","<<students[i].average<<","<<Grades[i]<<endl;
      }

      else if(Grades[i] == 'C')
      {
      out_file<<students[i].studentName<<","<<students[i].average<<","<<Grades[i]<<endl;
      }


   }


      out_file.close();

   return 0;
}