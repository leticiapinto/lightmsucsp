
#pragma once

/* imagen a grafo en vecindad de 8 */
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "library/LElement.h"
#include "library/Graph.h"
#include "library/ImageOperations.h"

#include "HGB_Edward/HgbSegInterval.h"

using namespace cv;
using namespace std;

typedef  pair<int, int> iPair;

template <class T>
void print_vector(vector <T> myvec)
{
	for(int i = 0; i < myvec.size(); i++)
    {
        cout<< myvec[i]<<" ";
    }
    cout<<endl;
}


int getEdge(int node);

int weight_node(vector <double>  &altitude, int node);

void canonizeQBT(vector< Element* > &qbt_parents, vector< Element* > &qct_parents, vector <double>  &altitude, int vertex_number, int QCT_size);

int level_node(vector< Element* > &parents, int q );

void iterate_pre(vector< Element* > &parents, vector<int> &euler, vector<int> &level, int i);

int LCA(int u, int v, vector<int > &R, vector<int > euler, vector<int > level );

bool isleaf(const vector< Element* > &parents, vector<int> &level_down, int node);

void fathertobottom(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node);

void fathertotop(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int father, int root);

bool reducenode(vector< Element* > &parents, int node);

void removefromfather(vector< Element* > &parents, int node, int father );

void updatenodeinformation(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node, int root);

void count_leaves(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int leaves_pixels, int root);

void childrenFilter(vector< Element* > &parents, vector<int> &leaves, int leaves_pixels, int parameter_number);

void iterate_cut(vector< Element* > &parents, int i, int evaluate_node, int color, vector<int> &myflag, vector<int> &vector_color, int leaves_pixels);

void iterate_coloring(vector< Element* > &parents, vector<int> level_down, vector<int> vector_color, int node, int color);

void ColorRegions(vector< Element* > &parents, vector<int> &R, vector<int> &vector_color, vector<int> &level, int max, int cut_tree, int leaves_pixels, bool option);

void visititeratead( const vector< Element* > &parents, vector<int> &newparents, vector<int> &vector_color, vector<int> &myflag,  vector<int> &level_down, int color, int node);

int cutleveltree(const vector< Element* > &parents, vector< vector<int> > &newregions,  vector<int> &level_down,  vector<int> &vector_color, int max, int cut_tree, int leaves_pixels);

int ColorRegionsHGB(const vector< Element* > &parents, vector<int> &vector_color, vector<int> &level_down, int max, int cut_tree, int leaves_pixels, bool option, int width, int height);

int ColorRegionsHGB2(string filename, int op, int level_tree, int typeComputeD, const vector< Element* > &parents, vector<int> &vector_color, vector<int> &level_down, int max, int cut_tree, int leaves_pixels, bool option, int width, int height);

void updateTree(vector< Element* > &parents, vector <int> &euler, vector <int> &level, vector <int> &R, vector <int> &leaves, vector <int> &level_down, int leaves_pixels, int root );

Graph saliencymapQCT(Graph g, vector<int> &vector_color, vector< Element* > &parents, int leaves_pixels, int level_tree ); // leaves_pixels : cantidad de hojas 

Mat ProcessGraph(int type, Mat image, string mygraphfile, string finalgraphfile, int level_tree);

// Depth-first preprocessing
int32_t LCApreprocessDepthFirst(JCctree *CT, int32_t node, int32_t depth, int32_t *nbr, int32_t *rep, int32_t *Euler, int32_t *Represent, int32_t *Depth, int32_t *Number);

int32_t ** LCApreprocess(JCctree *CT, int32_t *Euler, int32_t *Depth, int32_t *Represent, int32_t *Number, int32_t *nbR, int32_t *lognR);

int32_t LowComAncFast(int32_t n1, int32_t n2, int32_t *Euler, int32_t *Number, int32_t *Depth, int32_t **Minim);

void autoremovenodes(vector< Element* > &parents, vector<int> &level_down, Element* q );

void sendtoRoot(vector< Element* > &parents, int node, int root);

void iterate_remove(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node, int root);

void removePixels(vector< Element* > &parents, vector<int> &level_down, vector<int> &leaves, int node, int root);

void regionFilter(vector< Element* > &parents, vector<int> &leaves, int leaves_pixels, vector<int> &level_down, int width, int height);

Mat sm_code(string mygraphfile, int level_tree);

void usage(char *arg);

inline int mini(int a, int b);

inline int maxi(int a, int b);

inline float minf(float a, float b);

inline float maxf(float a, float b);

inline double mind(double a, double b);

inline double maxd(double a, double b);

inline void * myCalloc(size_t num, size_t size);

int * computeNodeArea(JCctree *CT);

void saliency(JCctree *ST, graphe *g, double * STAltitudes);

JCctree * UpperBoundBenjamin(graphe *g, double ratio);

Mat HGB_Edward(string myfile, string imagefile, int op, string finalgraphfile, int level_tree, int typeComputeD);

Mat generateboundbox2(string filename, int op, int level_tree, int typeComputeD );

void generateboundbox3(string filename, int op, int level_tree, int typeComputeD );

