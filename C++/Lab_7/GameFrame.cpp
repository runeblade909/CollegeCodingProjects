#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

// Zak Reichenbach
// 11/5/2019
// I tried to work a little on things I think I know how to do reasonablly


// Function to roll dice
int roll()
{
return rand() %6 + 1;
}

// function to play the game
int oneTurn(bool turn, int playerPts, int computerPts)
{
int rollD, turnTotal, ch;

//Players turn
if(turn == 1)
{
cout << "\n\n Player Turn: " << endl;
// Roll dice
rollD = roll();

turnTotal=0;
ch=0;

while(ch==0)
{
cout << "\n Roll Outcome:" << rollD << endl;
// Check dice for 3's
if(rollD == 3)
{
// Stealing points
playerPts = playerPts + 3;
computerPts -= 3;
break;
}
else
{
do
{
// Adding to turn turnTotal
turnTotal += rollD;
//Roll dice
rollD = roll();
cout << "\n Roll Outcome: " << rollD << endl;
}while(rollD != 3);
}

cout << "\n 0-roll 1-Pass: ";
cin >> ch;
}

//Stealing points
playerPts = playerPts + (turnTotal);
computerPts -= turnTotal;

return playerPts;
}
else
{
cout << "\n\n Computer Turn: " << endl;    

//Roll dice
rollD = roll();

turnTotal = 0;
ch = 0;

while(ch == 0)
{
cout << "\n Roll Outcome: " << rollD << endl;
// Check roll
if (rollD == 3)
{
//Stealing points
computerPts = computerPts + 3;
playerPts -= 3;
break;
}
else
{
do
{
//Add to turn turnTotal
turnTotal += rollD;
//roll dice
rollD = roll();
cout << "\n Roll Outcome: " << rollD << endl;    
}while(rollD != 3);
}

ch = rand()%2;
}

//Steal points
computerPts = computerPts + turnTotal;
playerPts -= turnTotal;

return computerPts;
}
}

//Loop the game
bool loopGame (int winningPoints)
{
int playerPts = 100, computerPts = 100;

//Decide Turn
int turn =rand()%2;

//Loop till someone wins
while(playerPts < winningPoints && computerPts < winningPoints)
{
if(turn ==1)    
{
playerPts = oneTurn(turn, playerPts, computerPts);    
//Give turn to player
turn = 0;
cout << "\n Player Points: " << playerPts << endl;
}
else
{
computerPts = oneTurn(turn, playerPts, computerPts);
//Give turn to computer
turn = 1;
cout << "\n Computer Points: " << computerPts << endl;
}
}

if(playerPts > computerPts)
{
return 0;
}
else
{
return 1;
}
}

//Main game frame
int main()
{
int winningPoints;

//Seed rand
srand(time(NULL));

//Read winning score
cout << "\n Enter the number of points for victory (has to be over 100): ";
cin >> winningPoints;

// Victory
if(!loopGame(winningPoints))
{
cout << "\n Player won the game!" << endl;
}
else
{
cout << "\n Computer won the game!" << endl;
}

return 0;
}