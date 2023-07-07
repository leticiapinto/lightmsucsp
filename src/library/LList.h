#pragma once
#include<iostream>
using namespace std;
struct LNode
{
    LNode(int v){
        this->value = v;
    }
    int value;
    LNode* next;
    LNode* prev;
};

struct LList
{
    LNode* first;
    LNode* last;
    
    void insert(LNode* item)
    {
        if(this->last==NULL)
        {
            this->first = item;
            this->last = item;
            return;            
        }

        LNode* last_copy = this->last;
        last_copy->next = item;
        item->prev = last_copy;    
        this->last = item;

    }

    void remove(LNode* item)
    {
        //caso que este solito
        LNode* temp_prev = item->prev;
        LNode* temp_next = item->next;
        if(item == first && item == last)
        {
            first= NULL;
            last = NULL;
        }else if(item == first && item != last)//caso que sea el primero
        {
            temp_next->prev = NULL;
            first = first->next;
        }else if(item == last && item != first) //caso que sea el ultimo
        {
            temp_prev->next = NULL;
            last = last->prev;
            
        }else{
             temp_next->prev = temp_prev;
             temp_prev->next = temp_next;
        }
        
        
                
        delete item;
    }


    void print()
    {
        LNode* iterator = first;
        while(iterator!=NULL)
        {
            cout<< iterator->value<<", ";
            iterator = iterator->next;
        }
        cout<<endl;
    }

};

