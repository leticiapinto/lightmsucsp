#pragma once
#include<iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>



using namespace cv;
using namespace std;

typedef  pair<int, int> iPair;
int last_altitude=-1, cnt_altitude=-1;

// Structure to represent a graph
struct Graph
{
    int V, E;
    vector< Element* > parents_graph;
    vector< pair<double, iPair> > edges;
    vector< pair<double, iPair> > MST; 
    vector <double>  altitude;
    int cnt; 
    // Constructor
    Graph(int V, int E)
    {
        this->V = V;
        this->E = E;
    }
    // Utility function to add an edge
    void addEdge(int u, int v, double w)
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
    int QBT_size;
    int QT_size;
    int QCT_size;
    vector< Element* > parents;
    vector <double>  altitude;
    
    
    // Constructor.

    DisjointSets(int n, int QBT_size)
    {
        // Allocate memory
        this->n = n;
        this->QBT_size = QBT_size;
        //this->QT_size = QT_size;
        //this->QCT_size = QCT_size;
    }
    void qbt_makeSet(int q)
    {
        Element* parent1 = new Element(q);
        parents.push_back(parent1);
        altitude.push_back(-1);
        QBT_size++;
    }
    int qbt_findCanonical(int value)
    {
        while(parents[value]->parent!=NULL)
        {
            value = parents[value]->parent->index;
        }
        return value;
    }
    
    int qbt_union(int cx, int cy, double weight)
    {
        Element* parent1 = new Element(QBT_size);
        parent1->addChild(parents[cx]);
        parent1->addChild(parents[cy]);
        parents.push_back(parent1);
        /*
        if(weight!=last_altitude)
        {
            cnt_altitude++;
            last_altitude= weight;
        }
        altitude.push_back(cnt_altitude);
        */
        altitude.push_back(weight);
        QBT_size++; 

        return QBT_size -1;
    }
};


 /* Functions returns weight of the MST*/ 
int Graph::kruskalMST()
{
    int mst_wt = 0; // Initialize result
    int QBT_size = 0;
    int QT_size = 0;
    int QCT_size = 0;
    
    // Create disjoint sets

    
    DisjointSets ds(V, QBT_size);
    
    cnt = 0;
     
    for (int i = 0; i < V; i++) { 
        ds.qbt_makeSet(i);
    }

    //vector< pair<int, iPair> > edges;
    //edges = g.edges;

    sort(edges.begin(), edges.end());
    vector< pair<double, iPair> >::iterator it;
    for (it=edges.begin(); it!=edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        double w = it->first;
        
        int set_u = ds.qbt_findCanonical(u);
        int set_v = ds.qbt_findCanonical(v);
        
        if (set_u != set_v)
        {
            //cout<<endl<<"u : "<<u<<", v : "<<v<<", w : "<<w;
            MST.push_back({w, {u, v}});
            cnt++;
            mst_wt += w;
            ds.qbt_union(set_u, set_v, w);
        
        }
        
    }   
    
    parents_graph = ds.parents;
    altitude = ds.altitude;
    return mst_wt;
}


Graph readGraph(int &V, int &E, int &height, int &width,  string mygraphfile)
{
    
    Graph _g(V, E);
    string line;
    //ifstream myfile ("_1.graph");
    ifstream myfile (mygraphfile);

    if (myfile.is_open())
    {
        getline (myfile,line); 
        //cout << line << '\n';
        istringstream iss(line);
        vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
        height = stoi(results[1]); //rows
        width = stoi(results[3]); //cols
        //cout<< "rows  : " << height << endl;
        //cout<< "cols  : " << width << endl;

        getline (myfile,line);
        //cout << line << endl;
        
        istringstream iss2(line);
        vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
        V = stoi(results2[0]);
        E = stoi(results2[1]);
        //cout<< "cnt_vertices : "<<V<< ", cnt_edges : "<< E<<endl;
        getline (myfile,line);
        //cout << line << endl;

        for(int i  = 0; i < V; i++)
        {
            getline (myfile,line); // lectura de cada linea de vertices
            istringstream iss2(line);
            vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
            //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1])<<endl;  
        }
        //cout<<"***************************************************" <<endl;
        getline (myfile,line); // lectura de textos edges arcs
        //cout<<"line "<< line << endl;
        //cout<<"***************************************************" <<endl;
        Graph g(V, E);
        for(int i  = 0; i < E; i++)
        {
            getline (myfile,line); // lectura de cada linea de edges
            istringstream iss2(line);
            vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
            //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1]) << " " << stoi(results2[2]) << endl;
            
            g.addEdge(stoi(results2[0]), stoi(results2[1]), stoi(results2[2]));
        }
        _g = g;
        myfile.close();

    }else cout << "Unable to open file"; 
    
    return _g;
}


Graph saveGraph(Mat image, int &V, int &E, int &height, int &width)
{
    width = image.cols;
    height = image.rows;
    cout<<"width: "<<width <<", height : "<<height<<endl;
    int cant= (width-1)*(height-1)+2 + (width-1)*2 + (height-1); 
     V = width*height, E = cant;

    Graph g(V, E);
    for(int i = 0; i < height-1; i++)
    {
        for(int j = 0; j < width-1; j++)
        {
            int pixel_center =i*width + j;
            Vec3f pixel = image.at<Vec3b>(i, j);
            int b_ , g_, r_;
            b_ = pixel[0];
            g_ = pixel[1];
            r_ = pixel[2];

            //vecindad de 4
            int pixel_actual;
            Vec3f pixel_rgb_actual;

            pixel_actual = i*width+j+1;
            pixel_rgb_actual = image.at<Vec3b>(i,j+1);
            //cout<<"b : "<<pixel_rgb_actual[0]<<", g : "<<pixel_rgb_actual[1]<<", r : "<<pixel_rgb_actual[2]<<endl;
            int b_actual = abs(b_ - pixel_rgb_actual[0]);
            int g_actual = abs(g_ - pixel_rgb_actual[1]);
            int r_actual = abs(r_ - pixel_rgb_actual[2]);
            int w = int(sqrt(b_actual*b_actual + g_actual*g_actual + r_actual*r_actual ));
            //cout<<"pc : "<<pixel_center<<", pa : "<<pixel_actual<<", w : [ "<<w<<" ] \t";
            g.addEdge(pixel_center, pixel_actual, w);


            pixel_actual = (i+1)*width+j;
            pixel_rgb_actual = image.at<Vec3b>(i+1,j);
            //cout<<"b : "<<pixel_rgb_actual[0]<<", g : "<<pixel_rgb_actual[1]<<", r : "<<pixel_rgb_actual[2]<<endl;
            b_actual = abs(b_ - pixel_rgb_actual[0]);
            g_actual = abs(g_ - pixel_rgb_actual[1]);
            r_actual = abs(r_ - pixel_rgb_actual[2]);
            w = int(sqrt(b_actual*b_actual + g_actual*g_actual + r_actual*r_actual ));
            //cout<<"pc : "<<pixel_center<<", pa : "<<pixel_actual<<", w : [ "<<w<<" ] \t";
            g.addEdge(pixel_center, pixel_actual, w);
        }
        int j = width-1;
        int pixel_center =i*width + j;
        Vec3f pixel = image.at<Vec3b>(i, j);
        int b_ , g_, r_;
        b_ = pixel[0];
        g_ = pixel[1];
        r_ = pixel[2];

        //vecindad de 4
        int pixel_actual;
        Vec3f pixel_rgb_actual;

        pixel_actual = (i+1)*width+j;
        pixel_rgb_actual = image.at<Vec3b>(i+1,j);
        //cout<<"b : "<<pixel_rgb_actual[0]<<", g : "<<pixel_rgb_actual[1]<<", r : "<<pixel_rgb_actual[2]<<endl;
        int b_actual = abs(b_ - pixel_rgb_actual[0]);
        int g_actual = abs(g_ - pixel_rgb_actual[1]);
        int r_actual = abs(r_ - pixel_rgb_actual[2]);
        int w = int(sqrt(b_actual*b_actual + g_actual*g_actual + r_actual*r_actual ));
        //cout<<"pc : "<<pixel_center<<", pa : "<<pixel_actual<<", w : [ "<<w<<" ] \t";
        g.addEdge(pixel_center, pixel_actual, w);
    }

    int i = height-1;
    for(int j = 0; j < width-1; j++)
        {
            int pixel_center =i*width + j;
            Vec3f pixel = image.at<Vec3b>(i, j);
            int b_ , g_, r_;
            b_ = pixel[0];
            g_ = pixel[1];
            r_ = pixel[2];

            //vecindad de 4
            int pixel_actual;
            Vec3f pixel_rgb_actual;

            pixel_actual = i*width+j+1;
            pixel_rgb_actual = image.at<Vec3b>(i,j+1);
            //cout<<"b : "<<pixel_rgb_actual[0]<<", g : "<<pixel_rgb_actual[1]<<", r : "<<pixel_rgb_actual[2]<<endl;
            int b_actual = abs(b_ - pixel_rgb_actual[0]);
            int g_actual = abs(g_ - pixel_rgb_actual[1]);
            int r_actual = abs(r_ - pixel_rgb_actual[2]);
            int w = int(sqrt(b_actual*b_actual + g_actual*g_actual + r_actual*r_actual ));
            //cout<<"pc : "<<pixel_center<<", pa : "<<pixel_actual<<", w : [ "<<w<<" ] \t";
            g.addEdge(pixel_center, pixel_actual, w);
        }



    return g;
}

void saveGraphFile(string myGraphFile, Graph g, int rs, int cs)
{

    //imprimir la primera linea 
    int n = g.V, m = g.E;
    ofstream myfile;
    myfile.open (myGraphFile);
    if (myfile.is_open())
    {
        myfile << "#rs "<< rs <<" cs "<< cs<<endl;
        myfile << n <<" "<< m << endl;
        myfile << "val sommets"<< endl;
        for (int i = 0; i < n; i++)
        myfile << i << " " << "1" << endl;
        myfile << "arcs values" << endl;
        vector< pair<double, iPair> >::iterator it;
        for (it=g.edges.begin(); it!=g.edges.end(); it++) //iterar sobre graph
        {
            int u = it->second.first;
            int v = it->second.second;
            double w = it->first;
            myfile << u << " " << v << " " << w << endl;   
        }

        myfile.close();

    }else cout << "Unable to open file"; 


}