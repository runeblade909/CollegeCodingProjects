/****************************************************************/
/*                  Stack Calculator Driver File                  */
/****************************************************************/
/*      TODO: Implement driver as described in the writeup      */
/****************************************************************/

#include "StackCalculator.hpp"
#include "StackCalculator.cpp"
#include <iostream>


using namespace std;

/*
 * Purpose: Determine whether some user input string is a
 *          valid floating point number
 * @param none
 * @return true if the string s is a number
 */
bool isNumber(string s)
{
    if(s.size() == 1 && s == "-") return false;
    else if(s.size() > 1 && s[0] == '-') s = s.substr(1);

    bool point = false;
    for(int i = 0; i < s.size(); i++)
    {
      if(!isdigit(s[i]) && s[i] != '.') return false;
      if(s[i]=='.' and !point) point = true;
      else if(s[i]=='.' and point) return false;
    }

    return true;
}

int main()
{
  // TODO - Declare a stack to hold the operands
  StackCalculator calculate;
  string input;

  cout << "Enter the operators and operands ('+', '*') in a postfix format" << endl;

  // End the math with an equals sign

  while(input != "=")
  {
       /* TODO
       1. Read input (operators and operands) until you encounter a "="
       2. Use the isNumber function to check if the input is an operand
       3. push all operands to the stack
       4. If input is not an operand, call the compute function to evaluate
          the partial expression
       */



    cout << "#> ";


    getline(cin, input);
    if(input != "="){
            if(isNumber(input)){
                calculate.push(stof(input));
            }
            else if(!isNumber(input)){
                calculate.calculate(input);

                if(calculate.isEmpty()){
                    cout << "No operands: Nothing to evaluate" << endl;
                    return 0;
                }
            }
    }
  }

  /* TODO - If the stack is empty then print "No operands: Nothing to evaluate" */

  /* TODO - Validate the expression
      1. If valid then print the result
      2. Else, print "Invalid expression"*/


     if(calculate.isEmpty()){
        cout << "No operands: Nothing to evaluate" << endl;
        return 0;
      }
      
    Op *temp = calculate.peek();
    float num = temp->num;
    calculate.pop();

    if(!calculate.isEmpty())
        cout << "Invalid expression" << endl;
    else
        cout << num << endl;

  return 0;

}
