#pragma once
#include<iostream>
using namespace std;

struct Element
{
    int index;
    int rnk;
    bool flag = false;
    int in_tree = 1;
    Element* parent = NULL;
    vector<Element* > children;
    //list<Element* > children;

    Element( int index )
    {
        this->index = index;
    }

    void addChild( Element* child )
    {
        if(child->parent != this)
        {
            child->parent = this;
            children.push_back(child);
        }else{
            children.push_back(child);
        }   
        
        //cout<<"child->index : "<<child->index<<endl; //cout<<"child->parent->index"<<child->parent->index<<endl;
    }
    void removeChild(Element* child )
    {
        int value = child->index;
        
        for(int i = 0; i < children.size(); i++)
        {
            if(children[i]->index == value)
            {
                children.erase (children.begin()+i);
            }
        }
    }
    void removeAllChildren()
    {
        children.clear();
    }

    
};


void print(vector< Element* > parents)
{
    //cout<<"qt_parents"<<endl;
    //for (int i = 0; i < parents.size(); i++) 
    for (int i = 0; i < parents.size(); i++) 
    { 
        if(parents[i]->parent!=NULL)
        {
            cout <<((parents[i]->parent)?"Has parent":"Not has parent");
            cout <<" ("<<parents[i]->index<<"-> " <<parents[i]->parent->index<<")";            
        }else{
            cout<<parents[i]->index<<" has not  parent";
        }
        cout <<", children: ";
        
        for(int j =0; j < parents[i]->children.size(); j++)
        {
            cout <<parents[i]->children[j]->index<<",";
        }
        cout<<endl;
        //cout<<"fin"<<endl;
    }
}