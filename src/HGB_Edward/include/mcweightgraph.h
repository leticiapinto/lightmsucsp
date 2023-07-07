/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/* authors : J. Cousty - L. Najman and M. Couprie */

/* $Id: mcgraphe.h,v 1.3 2006/02/28 07:49:12 michel Exp $ */

/* ================================================ */
/* types publics */
/* ================================================ */
#pragma once
#define MAXEDGES 2000000
#define MAXDOUBLE FLOAT_MAX

#define MINDOUBLE (1e-999)

#include "MST.h"

//#include "cv.h"
//#include "highgui.h"
//#include "mcunionfind.h"


typedef struct cell {
  //  u_int32_t val;
  u_int32_t index;
  struct cell * next;

  /*
  void initCell(struct cell *head, u_int32_t n){
    head->index = n;
    head->next = NULL;
  }

  void addCell(struct cell *head, u_int32_t n){
    struct cell *NewCell = new cell;
    NewCell-> index = n;
    NewCell -> next = head;
    head = NewCell;
  }*/

} cell;
typedef cell * pcell;

typedef struct {
  u_int32_t v_sommets;
  u_int32_t nbsom;
  u_int32_t nbmaxar;
  u_int32_t nbar; 
  u_int32_t ind;  /*number of edges*/
  double * weight; /*edge weight indexed by cell index*/
  double * originalWeight; /*edge weight indexed by cell index*/
  u_int32_t *tete;
  u_int32_t *queue;
  u_int32_t cs;
  u_int32_t rs;
  cell * tasar;    /* tableau des cellules (tas) */
  cell * librear;  /* liste des cellules libres geree en pile lifo */
  pcell * listar;  /* tableau des listes d'aretes indexe par les sommets */
} graphe;

/* ================================================ */
/* prototypes */
/* ================================================ */

/*
extern graphe * initgraphe(u_int32_t nbsom, u_int32_t nbmaxar);
extern void terminegraphe(graphe * g);
extern pcell allouecell(graphe * g);
extern void liberecell(graphe * g, pcell p);
extern void retiretete(graphe * g, pcell * pliste);
extern void retirearete(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t estarete(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t indexarete(graphe * g, u_int32_t som, u_int32_t a, u_int32_t ind);
extern void ajoutearete(graphe * g, u_int32_t som, u_int32_t a);
extern void addarete(graphe * g, u_int32_t som, u_int32_t a, double weight );
extern void setweight(graphe * g, u_int32_t index, double value);
extern void maille4graphe(graphe * g, u_int32_t rs, u_int32_t cs);
//extern void nettoiegraphe(graphe * g);
extern void aretesommets(u_int32_t a, u_int32_t N, u_int32_t rs, u_int32_t * s1, u_int32_t * s2);
extern int32_t estsuccesseur(graphe * g, u_int32_t som, u_int32_t a);
extern int32_t estsymetrique(graphe * g);
extern double getweight(graphe * g, u_int32_t index);
extern u_int32_t voisin(graphe *g, u_int32_t x, u_int32_t u);
extern graphe * ReadGraphe(char * filename,double **Fv);
extern void SaveGraphe(graphe * g, char *filename,double *Fv );
extern void setSize(graphe *g, int32_t rs , int32_t cs);

*/

/* ========================================== */
static void erreur(string mess)
/* ========================================== */
{
  //fprintf(stderr, "%s\n", mess);
  exit(0);
} /* erreur() */


/* ====================================================================== */
graphe * initgraphe(u_int32_t nbsom, u_int32_t nbmaxar)
/* ====================================================================== */
{
  graphe * g;
  u_int32_t i;
  
  g = (graphe *)calloc(1,sizeof(graphe));
  g->cs = -1;
  g->rs = -1;
  if (g == NULL)
  {   //fprintf(stderr, "initgraph : malloc failed for g\n");
      return(0);
  }
  g->tasar = (cell *)calloc(1,nbmaxar * sizeof(cell));
  if (g->tasar == NULL)
  {   //fprintf(stderr, "initgraph : malloc failed for g->tasar\n");
      return(0);
  }
  g->listar = (pcell *)calloc(nbsom, sizeof(pcell));
  if (g->listar == NULL)
  {   //fprintf(stderr, "initgraph : calloc failed for g->listar\n");
      return(0);
  }

  g->weight = (double *)calloc(1,nbmaxar * sizeof(double));
  g->originalWeight = (double *)calloc(1,nbmaxar * sizeof(double));
  g->nbsom = nbsom;
  g->nbmaxar = nbmaxar;
  g->nbar = 0;
  g->ind = 0;
  //printf("\n nbsom : %d, nbmaxar : %d \n", nbsom, nbmaxar);
  for (i = 0; i < nbmaxar - 1; i++)
    (g->tasar+i)->next = g->tasar+i+1;
  (g->tasar+i)->next = NULL;
  g->librear = g->tasar;  

  //printf("%d \n", (g->tasar)->index);

  g->tete = (u_int32_t*)calloc(g->nbmaxar, sizeof(u_int32_t));
  if(g->tete == NULL){
    //fprintf(stderr,"initgrap: calloc failed for g->tete\n");
    exit(1);
  }
  g->queue = (u_int32_t*)calloc(g->nbmaxar, sizeof(u_int32_t));
  if(g->queue == NULL){
    //fprintf(stderr,"initgrap: calloc failed for g->queue\n");
    exit(1);
  }
  return g;
} /* initgraphe() */


/* ====================================================================== */
pcell allouecell(graphe * g)
/* ====================================================================== */
{
  pcell p;
  if (g->librear == NULL) erreur("allouecell : plus de cellules libres");
  p = g->librear;
  g->librear = g->librear->next;
  return p;
} /* allouecell() */

/* ====================================================================== */
void setweight(graphe * g, u_int32_t index, double value)
/* ====================================================================== */
{
  g->weight[index]= value;
  
} /* setweight */


/* ====================================================================== */
void setweightOriginal(graphe * g, u_int32_t index, double value)
/* ====================================================================== */
{
    g->originalWeight[index]= value;
    
}



/* ====================================================================== */
double getweight(graphe * g, u_int32_t index)
/* ====================================================================== */
{
  return g->weight[index];
 
} /* setweight */

/* ====================================================================== */
double getweightOriginal(graphe * g, u_int32_t index)
/* ====================================================================== */
{
    return g->originalWeight[index];
    
} /* setweight */



/* ====================================================================== */
void ajoutearete2(graphe * g, u_int32_t som, u_int32_t a, u_int32_t ind)
/* ====================================================================== */
{
  pcell p;
  p = allouecell(g);
  p->next = *(g->listar + som);
  //  p->val = a;
  p->index = ind;
  *(g->listar + som) = p;
  g->nbar++;
} /* ajoutearete() */

/* ====================================================================== */
void addarete(graphe * g, u_int32_t som, u_int32_t a, double weight)
/* ====================================================================== */
{
u_int32_t i;
  
  i = g->ind;
  ajoutearete2(g, som, a,i);
/*  indexarete(g, som, a, i); */
  ajoutearete2(g, a, som,i);
/*  indexarete(g, a, som, i); */
  g->tete[i] = a;
  g->queue[i] = som;
  setweight(g, i, weight);
  g->ind= i+1;
  
} /* addarete() */

/* ====================================================================== */
void addEdgeWithOriginalW(graphe * g, u_int32_t som, u_int32_t a, double weight, double ori_weight)
/* ====================================================================== */
{
    u_int32_t i;
    
    i = g->ind;
    ajoutearete2(g, som, a,i);
    /* 	indexarete(g, som, a, i); */
    ajoutearete2(g, a, som,i);
    /* 	indexarete(g, a, som, i); */
    g->tete[i] = a;
    g->queue[i] = som;
    setweight(g, i, weight);
    
    setweightOriginal(g, i, ori_weight);
    
    g->ind= i+1;
} /* addarete() */





void setSize(graphe *g, int32_t rs , int32_t cs)
{
  g->rs = rs;
  g->cs = cs;
}


graphe * ReadGraphe(string mygraphfile, double **Fv)
{

    graphe *g;
    int i, n, m, t, q, rs = -1, cs = -1;

    string line;
    //ifstream myfile ("_1.graph");
    ifstream myfile (mygraphfile);

    if (myfile.is_open())
    {
        getline (myfile,line); 
        //cout << line << '\n';
        istringstream iss(line);
        vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
        rs = stoi(results[1]); //rows
        cs = stoi(results[3]); //cols
        //cout<< "rows  : " << rs << endl;
        //cout<< "cols  : " << cs << endl;

        getline (myfile,line);
        //cout << line << endl;
        
        istringstream iss2(line);
        vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
        n = stoi(results2[0]);
        m = stoi(results2[1]);
        //cout<< "cnt_vertices : "<<n<< ", cnt_edges : "<< m<<endl;

        g = initgraphe(n, 2*m);
        setSize(g,rs,cs);

        getline (myfile,line);  //linea que lee  "val sommets"
        //cout << line << endl;
        double *F = (double *)calloc(1,n * sizeof(double));

        for(int i  = 0; i < n; i++) //reading weight of vertices 
        {
            getline (myfile,line); // lectura de cada linea de vertices
            istringstream iss2(line);
            vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
            //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1])<<endl; 
            F[i] = stoi(results2[1]);  //F[i] = v

        }
        //cout<<"***************************************************" <<endl;
        getline (myfile,line); // lectura de textos edges arcs
        //cout<<"line "<< line << endl;
        //cout<<"***************************************************" <<endl;

        for(int i  = 0; i < m; i++) // read weights of arcs
        {
            getline (myfile,line); // lectura de cada linea de edges
            istringstream iss2(line);
            vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
            //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1]) << " " << stoi(results2[2]) << endl;
            
            addarete(g , stoi(results2[0]), stoi(results2[1]), stod(results2[2]));  //addarete(g , t, q, e);
        }
        myfile.close();
        (*Fv)=F;
        
    }else cout << "Unable to open file"; 

    return g;
}

/* ====================================================================== */
/*! \fn void SaveGraphe(graphe * g, char *filename) 
    \param g (entr�) : un graphe.
    \param filename (entr�) : nom du fichier �g��er.
    \brief sauve le graphe g dans le fichier filename. 
*/
void SaveGraphe(graphe * g, char *filename, double *Fv ) 
/* ====================================================================== */
#undef F_NAME
#define F_NAME "SaveGraphe"
{
  int i, j, n = g->nbsom, m = g->ind,a,b;
  pcell p;
  FILE * fd = NULL;
  double v;

  fd = fopen(filename,"w");
  if (!fd)
  {
    fprintf(stderr, "%s: cannot open file: %s\n", F_NAME, filename);
    return;
  }
  
  if( (g->cs != -1) && (g->rs != -1))
    fprintf(fd,"#rs %d cs %d\n",g->rs, g->cs);

  fprintf(fd, "%d %d\n", n, m);

  
  fprintf(fd, "val sommets\n");
  for (i = 0; i < n; i++) 
    fprintf(fd, "%d %g\n", i, Fv[i]);
  
  fprintf(fd, "arcs values\n");
  for (i = 0; i < m; i++) 
   {
     a = g->tete[i];
     b = g->queue[i];
     v = getweight (g,i);
            
     fprintf(fd, "%d %d %g\n", a, b, v);
   }
  
  fclose(fd);
} /* SaveGraphe() */


void mysavegraphe(graphe * g, string finalgraph, double *Fv )
{
  int i, j, n = g->nbsom, m = g->ind,a,b;
  pcell p;
  double v;
  ofstream myfile;
  myfile.open (finalgraph);
  if (myfile.is_open())
  {
    if( (g->cs != -1) && (g->rs != -1))
    {
      myfile << "#rs "<< g->rs <<" cs "<< g->cs<<endl;
    }

    myfile << n <<" "<< m << endl;
    myfile << "val sommets"<< endl;
    for (int i = 0; i < n; i++)
      myfile << i << " " << Fv[i]<<endl;
    myfile << "arcs values" << endl;
    for (i = 0; i < m; i++) 
    {
      a = g->tete[i];
      b = g->queue[i];
      v = getweight (g,i);
      myfile << a << " " << b << " " << v << endl;       
    }  
    myfile.close();

  }else cout << "Unable to open file"; 


}


/* ====================================================================== */
void printliste(pcell p)
/* ====================================================================== */
{
  //for (; p != NULL; p = p->next)
  //  printf("%d", p->index);
    
    
  //printf("\n");
} /* printliste() */

void printindex(graphe * g)
/* ====================================================================== */
{

u_int32_t i;
pcell p;
  
  for (i = 0; i < g->nbsom; i++)
  {
    printf("[%d] ", i);
    p = *(g->listar+i);
    //for (; p != NULL; p = p->next){
    //  printf("%d ", p->index);
    //printf("w=%d :", g->weight[p->index]);
    //}
  //printf("\n");
  } 
        
  //printf("\n");
} /* printindex() */

/* ====================================================================== */
void printgraphe(graphe * g)
/* ====================================================================== */
{
  u_int32_t i;
  
  // printf("[#id_vertex]: list of adjacent vertex \n", i);
  for (i = 0; i < g->nbsom; i++)
  {
    //printf("[%d]: ", i);
    //printliste(*(g->listar+i));
  }
  //printf("\n");

  // printf("[#id_edge]:  vertex conforming the edge \n", i);  
  for (i = 0; i < g->ind; i++)
  { 
     //printf("[%d] ", i);
     //printf("%d ", g->tete[i]); 
     //printf("%d ", g->queue[i]);
     //printf("\n");  
   }
} /* printgraphe() */


/* ====================================================================== */
void terminegraphe(graphe * g)
/* ====================================================================== */
{
  free(g->tasar);
  free(g->listar);
  free(g->weight);
  free(g->tete);
  free(g->queue);
  free(g);
} /* terminegraphe() */



/* ====================================================================== */
void liberecell(graphe * g, pcell p)
/* ====================================================================== */
{
  p->next = g->librear;
  g->librear = p;
} /* liberecell() */

/* ====================================================================== */
void retiretete(graphe * g, pcell * pliste)
/* ====================================================================== */
{
  pcell p;
  p = (*pliste)->next;
  liberecell(g, *pliste);
  *pliste = p;
} /* retiretete() */


  /////////////////////////////////////////
// lca : Nearest (Lowest) Common Ancestor
// 
// From: The LCA Problem Revisited 
// M.A. Bender - M. Farach-Colton
//
// from lwshedtopo.c

graphe * graphFromMST(graphe **outg, MST * mst, int nVertex, int nEdge, long int value)
/* ====================================================================== */
#undef F_NAME
#define F_NAME "graphFromMST"
{
#define TAILLEBUF 4096
    //graphe * g;
    //, n, m, t, q, rs = -1, cs = -1;
    //cout << "nVertex : " <<nVertex<<  ", nEdge : " <<  nEdge<< endl;
   
    *outg = initgraphe(nVertex, nEdge*2);
     
    //cout<<"llegue"<<endl;
    
    setSize(*outg,nVertex, nEdge);
    int outset = 10;
    for(int i=0;i<nEdge;i++){
        //addarete(g , mst->MSTedges[i].x, mst->MSTedges[i].y, mst->MSTedges[i].weight);
        //addEdgeWithOriginalW(*outg , mst->MSTedges[i].x, mst->MSTedges[i].y, INFINITO, mst->MSTedges[i].weight);
        
        addEdgeWithOriginalW(*outg , mst->MSTedges[i].x, mst->MSTedges[i].y, value, mst->MSTedges[i].weight+outset);
        mst->MSTedges[i].weight +=outset;
        //printf("MST edges  %c,%c, %f \n", 97+mst->MSTedges[i].y, 97+mst->MSTedges[i].x , mst->MSTedges[i].weight);
       
       // printf(" MST edges  %d,%d,\t %f \n", mst->MSTedges[i].x, mst->MSTedges[i].y , mst->MSTedges[i].weight);
    
    }
    //*outg = g;
    return *outg;
}



// leti 

graphe * ImagetoGraph(string myimagefile, double **Fv)
{

    //cout << "Starting the image to graph process <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< "<< endl;
    graphe *g;
    int i, n, m, t, q, rs = -1, cs = -1;

    Mat image = imread(myimagefile.c_str(), IMREAD_GRAYSCALE );

    //cout << "image.cols : "<< image.cols << ", image.rows "<< image.rows<< endl;


    rs = image.cols; //rows s[0]
    cs = image.rows; //cols s[1]
    //cout<< "rows  : " << rs << endl;
    //cout<< "cols  : " << cs << endl;

    n = rs*cs;
    m = (rs-1)*( (cs-1)*2 + 1) + cs -1;

    //cout<< "cnt_vertices : "<<n<< ", cnt_edges : "<< m<<endl;

    g = initgraphe(n, 2*m);
    setSize(g,rs,cs);

    double *F = (double *)calloc(1,n * sizeof(double));

    //cout<< "Reading weights of the vertices "<<endl;
    for(int i  = 0; i < n; i++) //reading weight of vertices 
    {
        F[i] = 1;  //F[i] = v  // It is 1 for the first step
    }

    //cout << "arcs values"<< endl;

    for(int i = 0; i < cs; i++)
    {
      for(int j = 0; j < rs; j++)
      {
        //cout << i*rs + j+1 << " " << i*rs +j << ", " << int (image.at<uchar>(i, j)) << " " << int (image.at<uchar>(i,j+1)) << endl; 
        //cout << (i+1)*rs + j << " " << i*rs +j << ", " << int (image.at<uchar>(i,j)) << " " << int (image.at<uchar>(i+1,j)) << endl<< endl; 

        //cout << int (image.at<uchar>(j,i)) << " ";
      }
      //cout << endl;
    }

    /*
    for()
    {

    }
    */

    /*
    for(int i  = 0; i < m; i++) // read weights of arcs
    {
        getline (myfile,line); // lectura de cada linea de edges
        istringstream iss2(line);
        vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
        //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1]) << " " << stoi(results2[2]) << endl;
        
        addarete(g , stoi(results2[0]), stoi(results2[1]), stod(results2[2]));  //addarete(g , t, q, e);
    }
    myfile.close();
    (*Fv)=F;
    */


    return g;

}

graphe * TexttoGraphe(string filename, double **Fv)
{

  //Need to read
  //Need to pass from pixeles to the weight of the edge --- go check the generator.py file to generate graph :) love you Kristopher! 
  //line 1 216	256
  //line 2 slice 	 rs 	 cs
  graphe *g;
  int i, n, m, t, q, rs = -1, cs = -1;

  string line;
  //ifstream myfile ("_1.graph");
  ifstream myfile (filename);
  // Read file that have the information from the pixeles 
  if (myfile.is_open())
  {
    getline (myfile,line); 
    //cout << line << '\n';
    istringstream iss(line);
    vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());
    rs = stoi(results[0]); //rows
    cs = stoi(results[1]); //cols
    //cout << rs << "\t" << cs << endl;

    getline (myfile,line); 
    //cout << line << '\n';
    istringstream iss2(line);
    vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());

    //cout << results2[0] << "\t" << results2[1] << "\t" << results2[2] << endl;
    // Saving the pixeles from the file to a pixeles matrix
    int pixeles[rs+1][cs+1];
    for(int i  = 0; i < rs; i++) //reading weight of vertices 
    {
      for(int j = 0; j < cs; j++)
      {
        getline (myfile,line); // lectura de cada linea de vertices
        istringstream iss2(line);
        vector<string> results2(istream_iterator<string>{iss2}, istream_iterator<string>());
        //cout<<" "<< stoi(results2[0]) << " " << stoi(results2[1])<<endl; 
        //F[i] = stoi(results2[1]);  //F[i] = v
        //cout << results2[1] << "\t" << results2[2] << "\t" << results2[3] << endl;
        pixeles[stoi(results2[1])][stoi(results2[2])] = stoi(results2[3]);
      }
    }
    //cout << "I finished the pixeles"<< endl;
    //I saved the pixeles now build the graph:
    //vertices = s[0]*s[1]
    //edges = (s[0]-1)*( (s[1]-1)*2 + 1) + s[1] -1

    //I need to obtain the values and save it to the graph immediately 
    //working with vertices
    rs = stoi(results[1]);
    cs = stoi(results[0]);
    int vertices = rs*cs;
    int edges = (rs-1)*( (cs-1)*2 + 1) + cs -1;
    int n = rs*cs;
    int m = (rs-1)*( (cs-1)*2 + 1) + cs -1;
    //cout<< "cnt_vertices : "<<n<< ", cnt_edges : "<< m<<endl;
    g = initgraphe(n, 2*m);
    setSize(g,rs,cs);

    //getline (myfile,line);  //linea que lee  "val sommets"
    //cout << line << endl;
    double *F = (double *)calloc(1,vertices * sizeof(double));
    //cout << "before the F matrix" << endl;
    for(int i = 0; i < vertices; i++ )
    {
      F[i] = 1;
    }
    //cout << "I passed the F matrix" << endl;
    //cout << "rs: " << rs << ", cs: " << cs << endl;
    //Now using the pixeles matrix I need to load the pixeles value to the graph

    rs = stoi(results[0]);
    cs = stoi(results[1]);
    for(int i  = 0; i < rs-1; i++) //reading weight of vertices 
    {
      for(int j = 0; j < cs-1; j++)
      {
        //cout << "i: " << i << ", j: " << j << endl;
        int indexh_i = i*cs + j+1;
        int indexh_j = i*cs + j;
        int weighth = abs(pixeles[i][j+1] - pixeles[i][j]);
        //cout << " "<< indexh_i << " " << indexh_j << " " << weighth << endl;
        addarete(g ,indexh_i, indexh_j, weighth);  //addarete(g , t, q, e);

        int indexv_i = (i+1)*cs + j;
        int indexv_j = i*cs + j;
        int weightv =  abs(pixeles[i+1][j] - pixeles[i][j]);
        addarete(g , indexv_i, indexv_j, weightv);  //addarete(g , t, q, e);
      }
      //cout << "finished the first loop" << endl;
      int j2 = cs -1;
      int indexv_i = (i+1)*cs + j2;
      int indexv_j = i*cs + j2;
      int weightv =  abs(pixeles[i+1][j2] - pixeles[i][j2]);
      addarete(g , indexv_i, indexv_j, weightv);  //addarete(g , t, q, e);
    }
    //cout << "finished the two loops"<< endl;
    int i = rs - 1;
    for(int j = 0; j < cs-1; j++)
    {
      int indexh_i = i*cs + j+1;
      int indexh_j = i*cs + j;
      int weighth = abs(pixeles[i][j+1] - pixeles[i][j]);
      addarete(g ,indexh_i, indexh_j, weighth);  //addarete(g , t, q, e);
    }

    //cout<< "I finished the for loop " << endl;
    
    //start reading the next 3 elements to finally 
    myfile.close();
  }else cout << "Unable to open file"; 

  
  return g;
}
