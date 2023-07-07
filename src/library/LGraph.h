#pragma once
#include<iostream>
using namespace std;

struct Graph
{
    int V, E;
    vector< pair<int, iPair> > edges;
 
    // Constructor
    Graph(int V, int E)
    {
        this->V = V;
        this->E = E;
    }
    // Utility function to add an edge
    void addEdge(int u, int v, int w)
    {
        edges.push_back({w, {u, v}});
    }
    // Function to find MST using Kruskal's MST algorithm
    int kruskalMST();
};

// To represent Disjoint Sets
struct DisjointSets
{

    int n;
    // Constructor.

    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
    }
    void qbt_makeSet(int value)
    {
        Element* parent1 = new Element(value);
        qbt_parents.push_back(parent1);
        QBT_size++;
    }
    int qbt_findCanonical(int value)
    {
        while(qbt_parents[value]->parent!=NULL)
        {
            value = qbt_parents[value]->parent->index;
        }
        return value;
    }
    
    int qbt_union(int cx, int cy)
    {
        Element* parent1 = new Element(QBT_size);
        parent1->addChild(qbt_parents[cx]);
        parent1->addChild(qbt_parents[cy]);
        qbt_parents.push_back(parent1);

        QBT_size++; 
    }

    //QT  ********************************************************
    void qt_makeSet(int value)
    {
        Element* parent1 = new Element(QT_size);
        parent1->rnk = 0;
        qt_parents.push_back(parent1);
        QT_size++;
    }

    int qt_findCanonical(int q)
    {
        int r = q;
        while(qt_parents[r]->parent!=NULL)
        {
            r = qt_parents[r]->parent->index;
        }
        while(qt_parents[q]->parent!=NULL)
        {
            int tmp = q;
            q = qt_parents[q]->parent->index;
            qt_parents[tmp]->parent = qt_parents[r];
        }
        return r;
    }

    int qt_union(int cx, int cy)
    {
        if(qt_parents[cx]->rnk > qt_parents[cy]->rnk)
        {
            int tmp = cx;
            cx = cy;
            cy= tmp;
        }
        if(qt_parents[cx]->rnk == qt_parents[cy]->rnk)
        {
            qt_parents[cy]->rnk +=1;
        }
        qt_parents[cx]->parent = qt_parents[cy];
        return cy;
    }
    //by altitute ordering

    void qebt_makeSet(int q)
    {
        qebt_root.push_back(q);
        qbt_makeSet(q);
        qt_makeSet(q);
        altitude.push_back(-1);
    }

    int qebt_union(int cx, int cy, int weight)
    {
        int tu = qebt_root[cx];
        int tv = qebt_root[cy];
        //union in QBT (without compression)
        qbt_parents[tu]->parent= qbt_parents[QBT_size];
        qbt_parents[tv]->parent= qbt_parents[QBT_size];

        // If children are needed, add them to the root
        //QBT .children[QBT .size].add({tu}); QBT .children[QBT .size].add({tv});
        Element* parent1 = new Element(QBT_size);
        parent1->addChild(qbt_parents[tu]);
        parent1->addChild(qbt_parents[tv]);
        qbt_parents.push_back(parent1);

        altitude.push_back(weight);
        int c = qt_union(cx, cy); //union in qt (with compression)
    
        qebt_root[c] = QBT_size;
        QBT_size++;

        return QBT_size -1;
    }

    int qebt_findCanonical(int q)
    {
        return qt_findCanonical(q);
    }

};

 /* Functions returns weight of the MST*/ 
int Graph::kruskalMST()
{
    int mst_wt = 0; // Initialize result
    QBT_size = 0;
    QT_size = 0;
    QCT_size = 0;
    // Create disjoint sets
    DisjointSets ds(V);
    cnt = 0;
    
    for (int i = 0; i < V; i++) 
    { 
        ds.qebt_makeSet(i);
        //cout <<((qt_parents[i]->parent!=NULL)?"true":"false")<<endl;
    }

    // Sort edges in increasing order on basis of cost
    sort(edges.begin(), edges.end());
    
    // Iterate through all sorted edges
    vector< pair<int, iPair> >::iterator it;

    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        int w = it->first;
        //cout <<"u : "<< u << ", v : "<<v;
        
        int set_u = ds.qebt_findCanonical(u);
        int set_v = ds.qebt_findCanonical(v);
        
        //cout << ", set_u :" << set_u << ", set_v : "<< set_v ;
        
        if (set_u != set_v)
        {
            //cout <<" --------- sera agregado : "<< u << " - " << v<< ", w : "<< w;
            MST.push_back({w, {u, v}});
            cnt++;

            // Update MST weight
            mst_wt += w;

            // Merge two sets
            //cout<<",  w : "<<w<<endl;
            ds.qebt_union(set_u, set_v,w);
            //print(qbt_parents);            
        }
        //cout<< endl;
    }
    //print(qbt_parents);     

    return mst_wt;
}
