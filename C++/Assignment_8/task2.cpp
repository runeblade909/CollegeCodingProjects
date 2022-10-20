#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<vector>


using namespace std;

int roll,choiceMode,choiceRoll,turn,turnTotal;

int main(){

srand(time(NULL));
cout << "Hey there! Your job is to get the highest number you can!\nI'll roll for how many turns you got, make the most of it! Try not to roll a 3!" << endl;
cout << "Are we playing with the computer or a human? (1=Human,0=Computer)" <<endl;
cin >> choiceMode;


if (choiceMode == 1){
  cout <<"Time to decide the number of turns."<<endl;
    roll =(1+(rand()%10));
    cout << "The set number of turns is: " << roll << endl;
    turn = roll;

    int j;
    for (j = 1; j<=turn ; j++){
        
        cout << "Roll(1) or Hold(0)?" << endl;
        cin >> choiceRoll;
       
        if (choiceRoll == 0){
            j =turn;
            turnTotal = turnTotal-roll;
        }

        else if (choiceRoll ==1){
    
          
          roll =  1+(rand()%6);
          cout << "Roll Outcome:\n"<< roll <<endl;

        }
        
        turnTotal+=roll;

                if (roll == 3){
            j = turn+1;
         turnTotal = 3;
                }
                
        cout << "Total: " << turnTotal << endl;
        }

}


    else if (choiceMode == 0){
     cout <<"Time to decide the number of turns."<<endl;
     roll =(1+(rand()%10));
     cout << "The set number of turns is: " << roll << endl;
     turn = roll;

    
        
        int i;
        for (i = 1; i<=turn; i++) {
            
            cout << "Roll(1) or Hold(0)?" << endl;

            choiceRoll = (rand()%6);
                if (choiceRoll >=1){
                    choiceRoll = 1;
                }

            cout << choiceRoll <<endl;
        
            if (choiceRoll == 0){
                i =turn;
                turnTotal = turnTotal-roll;
            }

            else if (choiceRoll ==1){
        
            
            roll =  1+(rand()%6);
            cout << "Roll Outcome:\n"<< roll <<endl;

            }
            
            turnTotal+=roll;

                    if (roll == 3){
                i = turn+1;
            turnTotal = 3;
                    }
                    
            cout << "Total: " << turnTotal << endl;



    }
}


return 0;
}