#include "StackCalculator.hpp"
#include <iostream>
#include <string>

using namespace std;


StackCalculator::StackCalculator()
{
    head = NULL;
}

StackCalculator::~StackCalculator()
{
    while(!isEmpty())
    {
        pop();
    }
}

bool StackCalculator::isEmpty()
{
    if(head == NULL)
        return true;
    else
        return false;
}

void StackCalculator::push(float number)
{
    Op *temp = new Op();
    temp->num = number;
    temp->next = head;
    head = temp;
}

void StackCalculator::pop()
{
    if(!isEmpty()){
        Op *temp = new Op();
        temp = head;
        head = head->next;
        delete temp;
    }
    
    else{
        cout << "Stack empty, cannot pop an item." << endl;
    }
}

Op* StackCalculator::peek()
{
    if(!isEmpty()){
        return head;
    }
    
    else{
        cout << "Stack empty, cannot peek." << endl;
        return NULL;
    }
}

bool StackCalculator::calculate(string symbol)
{
    float num1;
    float num2;
    float ans;

    if(symbol == "+"){
        if(isEmpty()){
            cout << "err: not enough operands" << endl;
            return false;
        }

        num1 = head->num;
        pop();

        if(isEmpty()){
            cout << "err: not enough operands" << endl;
            push(num1);
            return false;
        }

        num2 = head->num;
        pop();
        ans = num1+num2;
        push(ans);
        return true;

    }
    else if(symbol == "*"){
        if(isEmpty()){
            cout << "err: not enough operands" << endl;
            return false;
        }
        num1 = head->num;
        pop();

        if(isEmpty()){
            cout << "err: not enough operands" << endl;
            push(num1);
            return false;
        }

        num2 = head->num;
        pop();
        ans = num1*num2;
        push(ans);
        return true;
    }
    else if(symbol == "="){
        return true;
    }
    else{
        cout << "err: invalid operation" << endl;
        return false;
    }
}