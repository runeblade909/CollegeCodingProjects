#include <iostream>
#include <stdlib.h>
#include <time.h>

using namespace std;

//Function that rolls a dice
int roll()
{
return rand() % 6 + 1;
}

//Function that plays a game
int oneTurn(bool turn, int playerPts, int computerPoints)
{
int rollD, turnTotal, ch;

//Player turn
if(turn==1)
{
cout << "\n\n Player turn: \n";

//Rolling a dice
rollD = roll();

turnTotal = 0;
ch=0;

while(ch==0)
{
cout << "\n Roll Outcome: " << rollD << "\n";
//Checking roll
if(rollD == 3)
{
//Stealing points
playerPts = playerPts + 3;
computerPoints -= 3;
break;
}
else
{
do
{
//Adding to turn total
turnTotal += rollD;
//Rolling a dice
rollD = roll();
cout << "\n Roll Outcome: " << rollD << "\n";
}while(rollD != 3);
}

cout << "\n 0-roll 1-hold: ";
cin >> ch;
}

//Stealing points
playerPts = playerPts + (turnTotal);
computerPoints -= turnTotal;


return playerPts;
}
else
{
cout << "\n\n Computer turn: \n";

//Rolling a dice
rollD = roll();

turnTotal = 0;
ch=0;

while(ch==0)
{
cout << "\n Roll Outcome: " << rollD << "\n";
//Checking roll
if(rollD == 3)
{
//Stealing points
computerPoints = computerPoints + 3;
playerPts -= 3;
break;
}
else
{
do
{
//Adding to turn total
turnTotal += rollD;
//Rolling a dice
rollD = roll();
cout << "\n Roll Outcome: " << rollD << "\n";
}while(rollD != 3);
}

ch = rand()%2;
}

//Stealing points
computerPoints = computerPoints + turnTotal;
playerPts -= turnTotal;

return computerPoints;
}
}

//Function that loops the game
bool loopGame(int winningPoints)
{
int playerPts = 100, computerPoints = 100;

//Deciding turn
int turn = rand()%2;

//Loop till some player won the game
while(playerPts < winningPoints && computerPoints < winningPoints)
{
if(turn == 1)
{
playerPts = oneTurn(turn, playerPts, computerPoints);
//Giving turn to opponent
turn = 0;
cout << "\n Player Points: " << playerPts << " \n";
}
else
{
computerPoints = oneTurn(turn, playerPts, computerPoints);
//Giving turn to opponent
turn = 1;
cout << "\n Computer Points: " << computerPoints << " \n";
}
}

if(playerPts > computerPoints)
{
return 0;
}
else
{
return 1;
}
}

//Main function
int main()
{
int winningPoints;

//Seeding random number generator
srand(time(NULL));

//Reading winning points
cout << "\n Enter the number of points required to win the game: ";
cin >> winningPoints;

//Playing game
if(!loopGame(winningPoints))
{
cout << "\n Player won the game... \n";
}
else
{
cout << "\n Computer won the game... \n";
}

return 0;
}