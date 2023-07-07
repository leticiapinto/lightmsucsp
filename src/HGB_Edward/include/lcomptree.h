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


/*
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mcweightgraph.h>
#include <mcunionfind.h>
#include <mcsort.h>
#include <lcomptree.h>
*/


#include "JCousty.h"
#include "mcsort.h"
#include "mcunionfind.h"
#include "mcweightgraph.h"



/* ==================================== */
JCctree * componentTreeAlloc(int32_t N)
/* ==================================== */
#undef F_NAME
#define F_NAME "componentTreeAlloc"
{
  JCctree *CT;
  CT = (JCctree *)malloc(sizeof(JCctree)); 
  CT->tabnodes = (JCctreenode *)malloc(N * sizeof(JCctreenode));
  CT->tabsoncells = (JCsoncell *)malloc(2*N * sizeof(JCsoncell));
  CT->flags = NULL;
  // CT->flags = (uint8_t *)calloc(N, sizeof(char));
  // If needed, must be allocated separatly
  //  memset(CT->flags, 0, N);
  if ((CT == NULL) || (CT->tabnodes == NULL) || (CT->tabsoncells == NULL))
  { 
    fprintf(stderr, "%s : malloc failed\n", F_NAME);
    return NULL;
  }
  CT->nbnodes = N;
  CT->nbsoncells = 0;
  return CT;
} // componentTreeAlloc()

/* ==================================== */
void componentTreeFree(JCctree * CT)
/* ==================================== */
{
  free(CT->tabnodes);
  free(CT->tabsoncells);
  if(CT->flags != NULL)free(CT->flags);
  free(CT);
} // ComponentTreeFree()


/* ==================================== */
void componentTreePrint(JCctree *CT, double *Alt, int32_t root, int32_t *Att1)
/* ==================================== */
{
  JCsoncell *s; 
  char c;
  if(root <= CT->nbnodes/2) c = 'S'; else c = 'C';
  printf("---------------------\n%c%d\nAltitude: %lf\n", c, root, Alt[root]);
  if( Att1 != NULL )
    printf("Attribut 1: %d\n",Att1[root]);
  
  printf("Fils: ");

  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next){
    if(s->son <= CT->nbnodes/2) c = 'S'; else c = 'C';
    printf("%c%d (altitude %lf) ",c,s->son, Alt[s->son]);
  }
  printf("\n---------------------\n");
  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next)
    componentTreePrint(CT, Alt, s->son, Att1);
}

/* ==================================== */
void calculReversePointer(JCctree *CT, int32_t root)  
/* ==================================== */
{
  JCsoncell *s; 
  for(s = CT->tabnodes[root].sonlist; s != NULL; s = s->next) 
  {
    calculReversePointer(CT, s->son);
    CT->tabnodes[s->son].father = root;
  }
}

/* The algorithm is an original adaptation/improvment to edge-weighted
   graph of the algorithm of Najman and Couprie published in IEEE IP
   in 2006 */
int32_t saliencyTree(/* INPUTS */
		     graphe *g, 
		     /* OUTPUTS */
		     JCctree ** SaliencyTree, /* Component tree */
		     double *STaltitude  /* weights associated to the
					    nodes of the component
					    tree : Must be allocated
					    elsewhere */
		     )
{
  int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
  JCctree *ST;
  int32_t *clefs; 
  int32_t *STmap;
  JCsoncell * newsoncell1;
  JCsoncell * newsoncell2;
  Tarjan *T;
  int32_t taille = g->nbsom;
  int32_t narcs = g->ind;
  
  if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de malloc\n"); 
  }
   
  if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
    fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
    exit(0);
  }
  
  for(i = 0; i < narcs; i++)
    clefs[i] = i; 
  d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
  
  if( (ST = componentTreeAlloc((2*taille))) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
    exit(0);
  }
  ST->nbnodes = taille;
  ST->nbsoncells = 0;
  
  T = CreeTarjan(taille);
  for(i = 0; i < taille; i++)
  {
    STmap[i] = i;
    ST->tabnodes[i].nbsons = 0;
    ST->tabnodes[i].sonlist = NULL;
    TarjanMakeSet(T,i);
  }

  for(i = 0; i < narcs; i++){ 
    // for each edge of the MST taken in increasing order of altitude
    n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
    if(n1 != n2) {
      /* If the two vertices do not belong to a same connected
	 component at a lowest level */
      // Which component of ST n1 and n2 belongs to ?
      x1 = STmap[n1]; x2 = STmap[n2];    
      // Create a new component
      z = ST->nbnodes; ST->nbnodes++; 
      nbsoncellsloc = ST->nbsoncells;
      ST->tabnodes[z].nbsons = 2;
      // the altitude of the new component is the altitude of the edge
      // under consideration
      STaltitude[z] = getweight(g, clefs[i]);
      // add x1 and x2 to the lists of sons of the new component
      newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);  
      newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
      ST->tabnodes[z].sonlist = newsoncell1;
      newsoncell1->son = x1;
      newsoncell1->next = newsoncell2;
      newsoncell2->son = x2;
      newsoncell2->next = NULL;
      ST->tabnodes[z].lastson = newsoncell2;
      ST->nbsoncells += 2
;
    
      // then link n1 and n2
      k = TarjanLink(T, n1, n2);
      STmap[k] = z;
    }//if(...)
  }//for(...)

  
  for(i = 0; i < taille; i++)
    STaltitude[i] = 0; /* altitudes of the leaves is 0 */

  /* Construct father relationship */
  ST->tabnodes[ST->nbnodes-1].father = -1;
  ST->root = ST->nbnodes - 1;
  calculReversePointer(ST, ST->root); 


  (*SaliencyTree) = ST;


  //leti : quiero ver como esta conformado el arbol 

  //JCctreenode

  
  printf("Imprimiendo arbol en el saliencyTree\n");

  printf(" taille :  %d\n", taille );

  
  /* liberation de la memoire */
  free(STmap); 
  free(clefs);
  TarjanTermine(T);
  return ST->nbnodes - 1;
}

void compactRec(JCctree *CT, double *Val, int32_t root){
  JCsoncell *tmp;
  JCsoncell *f1 = CT->tabnodes[root].sonlist;
  JCsoncell *f2 = CT->tabnodes[root].lastson;
  
  if (f1 == NULL){ /* A leaf of the tree, that is a vertex of the initial graph */
    Val[root] = Val[CT->tabnodes[root].father];
  } else {
    compactRec(CT, Val, f1->son);
    compactRec(CT, Val, f2->son);

    if( (Val[root] == Val[f1->son]) && (f1->son > CT->nbnodes/2) ){
      /* Then f1 must be compacted */
      CT->tabnodes[root].sonlist = CT->tabnodes[f1->son].sonlist;
      tmp = CT->tabnodes[f1->son].lastson;
    } else tmp = f1;

    if( (Val[root] == Val[f2->son]) && (f2->son > CT->nbnodes/2) ){
      /* Then f1 must be compacted */
      tmp->next = CT->tabnodes[f2->son].sonlist;
      CT->tabnodes[root].lastson = CT->tabnodes[f2->son].lastson;
    } else  tmp->next = f2;
  }
}

int32_t componentTree(/* INPUTS */
		     graphe *g, 
		     /* OUTPUTS */
		     JCctree **CompTree, /* Component tree */
		     double *CTaltitude  /* weights associated to the
					    nodes of the component
					    tree : Must be allocated
					    elsewhere */
		      )
{
  saliencyTree(g, CompTree, CTaltitude);
  compactRec((*CompTree), CTaltitude, (*CompTree)->root);
  calculReversePointer((*CompTree), (*CompTree)->root);
  return 1;
}


/* returns the area of the component root and recursively compute in
   the table Surf the area of any descendent of root */
double surfComponents(JCctree *CT, double *Surf, int32_t root)
{
  JCsoncell *cell;
  if(root <= CT->nbnodes/2){
    return Surf[root];
  } else {
    Surf[root] = 0.0;
    for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
      Surf[root] += surfComponents(CT, Surf, cell->son);
    return Surf[root];
  }
}

/* returns the volume of the last slice of a component and recusively
   compute in the table Vol this value for any descendant of root. The
   table Vol is supposed to contain the surface of the
   components.  */
void volLastSlice(JCctree *CT, double *Alt, double *Vol, int32_t root)
{
  JCsoncell *cell;
  Vol[root] = Vol[root] * (Alt[(CT->tabnodes[root]).father] - Alt[root]);
  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
    if(cell->son > CT->nbnodes/2)
      volLastSlice(CT,Alt,Vol, cell->son);

}

double felzenswalbRec(JCctree *CT, double *Alt, double *Vol, int32_t root)
{
  JCsoncell *cell;
  double r,v;
  Vol[root] = Vol[root] * (Alt[(CT->tabnodes[root]).father] - Alt[root]);
  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
    if(cell->son > CT->nbnodes/2){
      v = felzenswalbRec(CT,Alt,Vol, cell->son);
      if (v > Vol[root]){
	Vol[root] = v;
      }
    }
  return Vol[root];
}

int32_t felzenswalbRec2(JCctree *CT, double *Alt, double *Vol, int32_t root, double v)
{
  JCsoncell *cell;
  double fmax = -1.0;
  int32_t smax = -1, nbmin = 0;
  

  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
    if( (cell->son > CT->nbnodes/2) && (Vol[cell->son] > fmax) ){ 
      fmax = Vol[cell->son];
      smax = cell->son;
    }

  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
    if(cell->son > CT->nbnodes/2)
      if(cell->son == smax)
	nbmin += felzenswalbRec2(CT, Alt, Vol, cell->son, smax);
      else nbmin += felzenswalbRec2(CT, Alt, Vol, cell->son, Vol[root]);
  if(fmax != -1){
    Vol[root] = 0.0;
    return nbmin;
  }
  else {
    Vol[root] = v;
    return 1;
  }
}


/* returns the volume of a component and recusively compute in the
   table Vol this value for any descendant of root. The table Vol is
   supposed to initially contain the surface of the components.  */
double volumeRec(JCctree *CT, double *Alt, double *Vol, int32_t root)
{
  JCsoncell *cell;
  double volR = 0;
  
  volR = 0;

  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next)
    if(cell->son > CT->nbnodes/2)
      volR += volumeRec(CT,Alt,Vol, cell->son);
  
  Vol[root] = volR + Vol[root] * (Alt[(CT->tabnodes[root]).father] - Alt[root]);
  return Vol[root];
}

/* Returns the "lowest" descendent of the component root and
   recursively compute in the table minSon the "lowest" descendent of
   any descendent of root. The term lowest refers to the total order
   made by the altitude of the vertices in the tree and their index in
   the graph. This value is required for computing the dynamics */
int32_t minSonComponents(JCctree *CT, int32_t *minSon, double *alt, int32_t root)
{
  JCsoncell *cell;
  int32_t tmp;
  if(root <= CT->nbnodes/2){
    minSon[root] = root;
    return root;
  } else {
    minSon[root] = root;
    for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next){
      tmp = minSonComponents(CT, minSon, alt, cell->son);
      if (alt[tmp] < alt[minSon[root]]) minSon[root] = tmp;
      else if( (alt[tmp] ==  alt[minSon[root]]) && (tmp < minSon[root]))
	minSon[root] = tmp;				   
    }
    return minSon[root];
  }
}

/* Compute recursively the dynamics of all minima that are descendent
   of root in the table Dyna. The value altPass is the value of the
   first pass point above root. Returns the number of minima contained
   in the subtree of root. */ 
int32_t dynaRec(JCctree *CT, int32_t *minSon, double *Alt, double *Dyna, int32_t root, double altPass)
{
  int32_t nbMin = 0, ms = minSon[root];
  JCsoncell *cell;
  
  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next){
    if(cell->son > CT->nbnodes/2){ /* Then son is a component of the tree */
      if(ms == minSon[cell->son])
	nbMin += dynaRec(CT, minSon, Alt, Dyna, cell->son, altPass);
      else
	nbMin += dynaRec(CT, minSon, Alt, Dyna, cell->son, Alt[root]);
    }
  }
  if( nbMin != 0){
    Dyna[root] = 0.0; 
    return nbMin;
  } else {
    Dyna[root] = altPass - Alt[ms];
    return 1;
  }
}

/* Compute the extinction value of a component and recursively the
   extinction values of all its descendent (stored in Att, thus
   override the previous values of Att) - extinction values of non
   leaves components is considered to be 0 - Returns the number of
   minima contained in the subtree of root */
int32_t extinction(JCctree *CT, double *Att, int32_t root)
{
  JCsoncell *cell;
  int32_t maxson = -1;
  double maxV = 0.0;
    
  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next){
    if(cell->son > (CT->nbnodes/2)){ // it is a component and not a vertex
      if(Att[cell->son] > maxV){
	maxV = Att[cell->son];
	maxson = cell->son;
      }
    }
  }
  
  if(maxson != -1) { 
    Att[maxson] = Att[root]; 
    Att[root] = 0.0;
  } else return 1;
  
  //  maxV = 0.0;
  maxson = 0;
  for(cell = CT->tabnodes[root].sonlist; cell != NULL; cell = cell->next){
    if(cell->son > (CT->nbnodes/2)){
     maxson += extinction(CT, Att, cell->son);
    }
  }
  return maxson;
}

/* Attribute to each vertex of the minima the order of this minima
   with respect to the extinctions values stored in Surf */
double *orderExtinctions(JCctree *CT, int32_t nbmin, double *Surf, double *Av)
{
  int32_t *Key;
  int32_t i, j, nbnodes = CT->nbnodes;
  double *Bitmap;

  Key = (int32_t *)calloc(nbmin, sizeof(int32_t));
  if(Key == NULL){
    //fprintf(stderr, "orderExtinctions : cannot allocated (%d bytes) memory", nbmin*sizeof(int32_t));
    exit(1);
  }
  
  j = 0;
  for(i = (nbnodes/2) + 1; i < nbnodes; i++)
    /* for all components of CT (ie, vertices of the initial graph excluded) */
    if(Surf[i] > 0){
      Key[j] = i;
      j++;
    }
  /* Then sort the j keys */
  d_TriRapideStochastique(Key, Surf,0,nbmin-1);
  
  Bitmap = (double *)calloc(nbmin, sizeof(double));
  if(Bitmap == NULL){
    //fprintf(stderr, "orderExtinctions : cannot allocated (%d bytes) memory", nbmin*sizeof(double));
    exit(1);
  }
  
  for(i =0; i < nbmin; i++) {
    Bitmap[i] = Surf[Key[i]];
    Surf[Key[i]] = (double)(i+1); /* the first value is 1, not 0, thus 1 is added to i */
  }
  for(i = 0; i <= nbnodes/2; i++){
    if( Surf[CT->tabnodes[i].father] != 0){
      /* Vertex i belongs to a minimum */
      Av[i] = Surf[CT->tabnodes[i].father];
    }
  }
  free(Key);
  return Bitmap;
}

int32_t lextinctionvalues(graphe *g, int32_t mode, double *Av, double **Bitmap)
{
  JCctree *CT;
  double *Val, *Att;
  int32_t i, nbmin;
  int32_t *MinSon;

  Val = (double *)calloc(2*g->nbsom, sizeof(double));
  if(Val == NULL){ 
    fprintf(stderr, "extinctionvalues: Cannot allocate memory\n"); 
    exit(1);
  }

  componentTree(g,&CT,Val);

  /*  componentTreePrint(CT,Val, CT->root, NULL); */
  

  Att = (double *)calloc(2*g->nbsom, sizeof(double));
  if(Att == NULL){ 
    fprintf(stderr, "extinctionvalues: Cannot allocate memory\n"); 
    exit(1);
  }
  /* Initialize the surface of the vertices with the values stored in Av */
  memcpy(Att, Av, sizeof(double) * g->nbsom);
  switch(mode/2){
  case 0:
    surfComponents(CT, Att, CT->root);  
    /*
      printf("########################\n");
      printf("########################\n");
      printf("########################\n");
      componentTreePrint(CT,Val, CT->root, Att);
    */
    nbmin = extinction(CT, Att, CT->root);
    
    break;
  case 1:
    MinSon = (int32_t *)calloc(2 * g->nbsom, sizeof(int32_t));
    if(MinSon == NULL){
      //fprintf(stderr, "lextinctionvalues: Cannot allocate memory (%d bytes)\n", 2 * g->nbsom * sizeof(int32_t)); 
      exit(1);
    }
    minSonComponents(CT, MinSon, Val, CT->root);
    nbmin = dynaRec(CT, MinSon, Val, Att, CT->root, Val[CT->root]);
    
    free(MinSon);
    break;
  case 2:
    surfComponents(CT, Att, CT->root);
    volumeRec(CT, Val, Att, CT->root);
    nbmin = extinction(CT, Att, CT->root);
    fprintf(stderr," volume done avec nbmin = %d\n",nbmin);
    break;
  case 3:
    surfComponents(CT, Att, CT->root);
    fprintf(stderr,"surfComponents done\n");
    felzenswalbRec(CT, Val, Att, CT->root);
    fprintf(stderr,"felzenswalbRec done\n");
    nbmin = felzenswalbRec2(CT, Val, Att, CT->root, Att[CT->root]); 
    fprintf(stderr,"felzenswalbRec2 done avec nbmin = %d\n",nbmin);
    break;
  default:
    fprintf(stderr,"lextinctionvalues: case mode == %d not yet implemented\n", mode);
    exit(0);
  }
  /*
  printf("########################\n");
  printf("########################\n");
  printf("########################\n");
  componentTreePrint(CT,Val, CT->root, Surf);  
  printf("nbmin = %d et nbsom %d \n", nbmin, 2*g->nbsom);
  */

  for(i = 0; i < g->nbsom; i++) Av[i] = 0.0;
  if( mode%2 == 0){ /* Then outputs are directly extinction values */
    for(i = 0; i < g->nbsom; i++){
      if( Att[CT->tabnodes[i].father] != 0){
	/* Vertex i belongs to a minimum */
	Av[i] = (double)Att[CT->tabnodes[i].father];
      }
    }
  } else { /* Then the outputs are orders of the minima with respect to extinction values */
    (*Bitmap) = orderExtinctions(CT, nbmin, Att, Av);
  }
  
  free(Val); free(Att);
  componentTreeFree(CT);
  return nbmin;
}


/* The algorithm is an original adaptation/improvment to edge-weighted
 graph of the algorithm of Najman and Couprie published in IEEE IP
 in 2006 / Its description is provided in Cousty, Najman, Perret
 ISMM2013 and in Najman, Cousty, Perret ISMM 2013. */
int32_t BPTMST(/* INPUTS */
              graphe *g,
              /* OUTPUTS */
              JCctree ** SaliencyTree, /* Component tree - BPTAO  */
              double *STaltitude,  /* weights associated to the
                                    nodes of the BPTAO: Must
                                    be allocated elsewhere */
              MST* mst
              )
{
    int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
    JCctree *ST;
    int32_t *clefs;
    int32_t *STmap;
    JCsoncell * newsoncell1;
    JCsoncell * newsoncell2;
    Tarjan *T;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    
    if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de malloc\n");
    }
    
    if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
        fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
        exit(0);
    }
    
    setSizeMST(mst, narcs);
    
    
    for(i = 0; i < narcs; i++)
        clefs[i] = i;
    
    
    
    d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
    
    
    
    
    if( (ST = componentTreeAlloc((2*taille))) == NULL){
        fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
        exit(0);
    }
    ST->nbnodes = taille;
    ST->nbsoncells = 0;
    T = CreeTarjan(taille);
    for(i = 0; i < taille; i++)
    {
        STmap[i] = i;
        ST->tabnodes[i].nbsons = 0;
        ST->tabnodes[i].sonlist = NULL;
        TarjanMakeSet(T,i);
    }
    
    for(i = 0; i < narcs; i++){
        // for each edge of the graph or MST taken in increasing order of altitude
        n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
        if(n1 != n2) {
            /* If the two vertices do not belong to a same connected
             component at a lowest level */
            // Which component of ST n1 and n2 belongs to ?
            x1 = STmap[n1]; x2 = STmap[n2];
            // Create a new component
            z = ST->nbnodes; ST->nbnodes++;
            nbsoncellsloc = ST->nbsoncells;
            ST->tabnodes[z].nbsons = 2;
            // the altitude of the new component is the altitude of the edge
            // under consideration
            STaltitude[z] = getweight(g, clefs[i]);
            // add x1 and x2 to the lists of sons of the new component
            newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);
            newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
            ST->tabnodes[z].sonlist = newsoncell1;
            newsoncell1->son = x1;
            newsoncell1->next = newsoncell2;
            newsoncell2->son = x2;
            newsoncell2->next = NULL;
            ST->tabnodes[z].lastson = newsoncell2;
            ST->nbsoncells += 2;
            
            // then link n1 and n2
            k = TarjanLink(T, n1, n2);
            STmap[k] = z;
            
            //printf("MST node z:%i edge x,y %i - %i , compl,compr %i,%i \n",z, g->tete[clefs[i]],g->queue[clefs[i]], n1, n2);
            mst->MSTedges[mst->size].x= g->tete[clefs[i]];
            mst->MSTedges[mst->size].y= g->queue[clefs[i]];
            mst->MSTedges[mst->size].weight = getweight(g, clefs[i]);//already calculated weight of edge
            mst->MSTedges[mst->size].idx = clefs[i];
            //printf("MST edge %i - %i :  %f \n", g->tete[clefs[i]], g->queue[clefs[i]],mst->MSTedges[mst->size].weight);
            
            
            mst->size+=1;
            
        }//if(...)
        else{
            // printf("Discarded edge x,y %i - %i , components %i,%i\n", g->tete[clefs[i]],g->queue[clefs[i]],n1,n2 );
        }
    }//for(...)
    
    

    
    (*SaliencyTree) = ST;
    
    /* liberation de la memoire */
    free(STmap);
    free(clefs);
    TarjanTermine(T);
    return ST->nbnodes - 1;
}





////////////////// Edward


/** Calculates MST, outputs JCctree St, MST
 */
int getMST(graphe *g, MST** mst ){
    
    double *STAltitudes;
    JCctree *ST;
    
    *mst = (MST*)malloc(sizeof(MST));
    
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr,"getQBT cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    
    // To obtain the binary partition tree by altitude ordering
    BPTMST(g, &ST, STAltitudes, *mst);
    free(STAltitudes);
    componentTreeFree(ST);
    return 0;
}


////////----------------------------- With Attributes ---------------


int initializeTreeAttributes(/* INPUTS */
                   graphe *g,
                   /* OUTPUTS */
                   JCctree **CompTree,/* CompTree :  initalized component tree - QFZ */
                   double  *STaltitude,
                   double * STAltitudesSegmentation,
                   int *STArea,
                   double *STMedia,
                   double *STMaxWeight
                   ){
    
    int i ,nbsoncellsloc;
    JCsoncell * newsoncell1, *lastchild;
    JCctree *CT = *CompTree; //component tree
    int32_t taille = g->nbsom;
    
    CT->nbnodes = taille;
    CT->nbsoncells = 0;
    
    for(i = 0; i < taille; i++){
        CT->tabnodes[i].nbsons = 0;
        CT->tabnodes[i].sonlist = NULL;
        CT->tabnodes[i].lastson = NULL;
        CT->tabnodes[i].nodeAsSon = NULL;
    }
    
    for(i = 0; i < taille; i++){
        STaltitude[i] = 0; /* altitudes of the leaves is 0 */
        STArea[i] = 1; /* altitudes of the leaves is 0 */
        STMedia[i] = 1;  /* media of the leaves is 0 by now  leticia */
        STMaxWeight = 0;
    }

    int z;
    z = CT->nbnodes; CT->nbnodes++;
    STaltitude[z] = 9999999999;
    CT->tabnodes[z].nbsons = taille;
    CT->tabnodes[z].nodeAsSon = NULL;
    
    //// Create root node z of the tree and attach all vertices as their initial children.
    CT->root = z;
    //cout<<"z : "<<z;
    nbsoncellsloc = CT->nbsoncells;
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = 0;
    newsoncell1->next = NULL;
    newsoncell1->last = NULL;
    
    CT->tabnodes[z].sonlist = newsoncell1;
    CT->tabnodes[z].lastson = newsoncell1;
    CT->tabnodes[z].father = -1;
    CT->tabnodes[0].father = z;
    CT->tabnodes[0].nodeAsSon = newsoncell1;
    CT->nbsoncells += 1;
    
    STArea[z]=STArea[z]+ STArea[0];
    STMedia[z]=STMedia[z] + STMedia[0];
    
    lastchild =newsoncell1;
    for(i = 1; i < taille; i++){
        nbsoncellsloc = CT->nbsoncells;
        newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
        newsoncell1->son = i;
        newsoncell1->next = NULL;
        STArea[z]=STArea[z]+ STArea[i];
        CT->tabnodes[i].father = CT->root;
        CT->nbsoncells += 1;
        CT->tabnodes[i].nodeAsSon = newsoncell1;
        newsoncell1->last = lastchild;
        lastchild->next =newsoncell1;
        lastchild = newsoncell1;
    }
    CT->tabnodes[z].lastson = lastchild;
    
    
    for(i =g->nbsom ; i< 2*g->nbsom; i++ ){
        STAltitudesSegmentation[i] = STaltitude[i] ;
    }
    
    
    return 0;
}


void findTransition(/* INPUTS */
                    JCctree *QFZ,
                    graphe *g,
                    double  *STaltitude,
                    int leafIdx, /* leaf index for Edge{x,y}*/
                    int edgeIdx,
                    /* OUTPUTS */
                    int *c,
                    int *p
                    ){
    
    int idxSon = leafIdx;
    int idxFather =  QFZ->tabnodes[idxSon].father;
    
   // if(STaltitude[idxSon] <= g->weight[edgeIdx]+1 && !(STaltitude[idxFather]<= g->weight[edgeIdx]+1)){
    double lambda = g->weight[edgeIdx]+1;
    //double lambda = g->weight[edgeIdx];
     if(STaltitude[idxSon] <= lambda && !(STaltitude[idxFather]<= lambda)){
        *c = idxSon; *p = idxFather;
        return;
    }
    else{
        findTransition(QFZ, g,STaltitude, idxFather , edgeIdx , c,p );
    }
    
}

void detach(JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            const int idxP,   // idxS is dettached from idxS
            const int idxS ){
    
    JCctree * CT = *CompTree;
    JCsoncell *node = CT->tabnodes[idxS].nodeAsSon;
    
    if(CT->tabnodes[idxP].sonlist == node && CT->tabnodes[idxP].lastson == node ){
        // We are dettaching the last son of idxP, we have to mark it to delete
        CT->tabnodes[idxP].sonlist = NULL;
        CT->tabnodes[idxP].lastson = NULL;
        CT->tabnodes[idxP].nbsons = 0;
       // printf("idxP %i se quedo sin hijos \n", idxP);
        return; //detach
    }
    if(CT->tabnodes[idxP].sonlist == node && CT->tabnodes[idxP].lastson != node ){
        // There are still nodes left
        CT->tabnodes[idxP].sonlist = node->next;
        node->next->last = NULL;
        CT->tabnodes[idxP].nbsons-=1;
        return;
    }
    else{
        // it is the last;
        if(CT->tabnodes[idxP].lastson == node){
            CT->tabnodes[idxP].lastson = node->last;
            node->last->next = NULL;
            CT->tabnodes[idxP].nbsons-=1;
        }// it is in the middle;
        else{
            node->last->next = node->next;
            node->next->last = node->last;
            CT->tabnodes[idxP].nbsons-=1;
        }
    }
}


void attach(/* INPUTS */
            JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            const int idxP,
            const int idxS //  idxS is attached to idxP
){
    JCctree * CT = *CompTree;
    //int  nbsoncellsloc;
    JCsoncell * newsoncell1 ;
    
    if(CT->tabnodes[idxS].father == idxP)return;//no need to attach as it is already son of idxP
    if(idxS == idxP){return;}//no attaching to its self}

    CT->tabnodes[idxP].nbsons += 1;
    
    detach(CompTree, STaltitudes , CT->tabnodes[idxS].father,  idxS);
    
    newsoncell1 = CT->tabnodes[idxS].nodeAsSon;
    newsoncell1->next = NULL;
    newsoncell1->last = CT->tabnodes[idxP].lastson;
    
    if(CT->tabnodes[idxP].lastson != NULL){
        CT->tabnodes[idxP].lastson->next = newsoncell1 ;
        CT->tabnodes[idxP].lastson = newsoncell1;
    }
    else{
        CT->tabnodes[idxP].lastson = newsoncell1;
        CT->tabnodes[idxP].sonlist = newsoncell1;
    }
    CT->tabnodes[idxS].father = idxP;
    
}



void mergeAttributes(/* INPUTS */
           JCctree **CompTree, /* Component tree - QFZ */
           double  *STaltitudes,
           int *STArea,
           double *STMaxWeight,
           int p1,
           int p2,
           int *np,
           int *area1,
           int *area2
           ){
    
    if(p1==p2){*np = p1; return;}
    
    JCctree * CT = *CompTree;
    int larger, smaller, i ;
    if(CT->tabnodes[p1].nbsons > CT->tabnodes[p2].nbsons){
        larger = p1;
        smaller = p2;
        
        //////// Attributes  //////
        int mp = STArea[larger] ;
        STArea[larger] = STArea[larger] + STArea[smaller] - *area1 - *area2;
        *area1 = mp;
        mp = STArea[smaller] ;
        STArea[smaller] = 0;
        *area2 = mp;
        
        if(STMaxWeight[smaller] > STMaxWeight[larger])
            STMaxWeight[larger] = STMaxWeight[smaller];
        ////////////////////////////
        
    }
    else{
        larger = p2;
        smaller = p1;
        
        //////// Attributes  //////
        int mp = STArea[larger] ;
        STArea[larger] = STArea[larger] + STArea[smaller] - *area1 - *area2;
        *area2 = mp;
        mp = STArea[smaller] ;
        STArea[smaller] = 0;
        *area1 = mp;
        
        if(STMaxWeight[smaller] > STMaxWeight[larger])
            STMaxWeight[larger] = STMaxWeight[smaller];
        
        ////////////////////////////
        
    }
    // Children of the smaller one will go into the larger one.
    JCsoncell * child = CT->tabnodes[smaller].sonlist;
    
    if(CT->tabnodes[larger].lastson){
        CT->tabnodes[larger].lastson->next = child ;
        child->last = CT->tabnodes[larger].lastson;
    }
    else{
        // in case larger does not have a last son
        CT->tabnodes[larger].sonlist = child ;
        printf("larger %i has  %i sons, merging with %i with %i sons \n", larger, CT->tabnodes[larger].nbsons,smaller , CT->tabnodes[smaller].nbsons );
    }
    
    CT->tabnodes[larger].lastson = CT->tabnodes[smaller].lastson;
    
    for(i=0; i< CT->tabnodes[smaller].nbsons; i++ ){
        CT->tabnodes[child->son].father = larger;
        CT->tabnodes[larger].nbsons += 1;
        child = child->next;
    }
    
    detach(CompTree, STaltitudes , CT->tabnodes[smaller].father,  smaller);
    
    

    
    CT->tabnodes[smaller].sonlist = NULL;
    CT->tabnodes[smaller].lastson = NULL;
    CT->tabnodes[smaller].father = smaller; //marked for delete
    CT->tabnodes[smaller].nbsons = 0;
    
    *np = larger;
    
}

void attachAttributes(/* INPUTS */
            JCctree **CompTree, /* Component tree - QFZ */
            double  *STaltitudes,
            double *STAltitudesSegmentation,
            int *STArea,
            double *STMedia,
            double *STMaxWeight,
            const int idxP,  //idxP representa el padre
            const int idxS, //  idxS is attached to idxP // representa el hijo 
            int *area
){
    JCctree * CT = *CompTree;
    //int  nbsoncellsloc;
    JCsoncell * newsoncell1 ;
    
    if(idxS == idxP){return;}//no attaching to its self}

    //////// Attributes  //////
    
    int mp = STArea[idxP] ;
    STArea[idxP] = STArea[idxP] - *area + STArea[idxS];
    *area = mp;

    if(STMaxWeight[idxS]>STMaxWeight[idxP])
        STMaxWeight[idxP] = STMaxWeight[idxS];
    ////////////////////////////
    
    if(CT->tabnodes[idxS].father == idxP){
        return;//no need to attach as it is already son of idxP
    }

    CT->tabnodes[idxP].nbsons += 1;

    
    detach(CompTree, STaltitudes , CT->tabnodes[idxS].father,  idxS);
    
    newsoncell1 = CT->tabnodes[idxS].nodeAsSon;
    newsoncell1->next = NULL;
    newsoncell1->last = CT->tabnodes[idxP].lastson;
    
    if(CT->tabnodes[idxP].lastson != NULL){
        CT->tabnodes[idxP].lastson->next = newsoncell1 ;
        CT->tabnodes[idxP].lastson = newsoncell1;
    }
    else{
        CT->tabnodes[idxP].lastson = newsoncell1;
        CT->tabnodes[idxP].sonlist = newsoncell1;
    }
    CT->tabnodes[idxS].father = idxP;
    

    //cout<<"until updating media : "<<endl;
    // update media based on the new sons of idxP

    STMedia[idxP] = (STMaxWeight[idxP]/(STArea[idxP]));

  }


void createNodeAttributes( /* INPUTS */
                JCctree **CompTree, /* Component tree - QFZ */
                double  *STaltitudes,
                double *STAltitudesSegmentation,
                int *STArea,
                double *STMedia,
                double *STMaxWeight,
                graphe *g,
                int edgeIdx, /* Edge index for {x,y}  */
                int c1,
                int c2,
                /* OUTPUTS */
                int *n,
                int *area1, int *area2){
    JCctree * CT = *CompTree;
    int z, nbsoncellsloc;
    double lambda;
    lambda = g->weight[edgeIdx]+1;
   //  lambda = g->weight[edgeIdx];
    
    if(STaltitudes[c1] ==  STaltitudes[c2] &&STaltitudes[c1] == lambda  && CT->tabnodes[c1].nbsons > 0 && CT->tabnodes[c2].nbsons >0 ){
        mergeAttributes(CompTree , STaltitudes, STArea, STMaxWeight, c1,c2, n, area1, area2);
        if( g->originalWeight[edgeIdx] > STMaxWeight[*n]){
            STMaxWeight[*n] = g->originalWeight[edgeIdx];
        }
        
        return ;
    }

    if(STaltitudes[c1] == lambda && STaltitudes[c2] < lambda ){
        attachAttributes(CompTree , STaltitudes, STAltitudesSegmentation, STArea, STMedia, STMaxWeight, c1 , c2, area1);
        
        *area2 = STArea[c2];
        // Max weight
        if( g->originalWeight[edgeIdx] > STMaxWeight[c1]){
            STMaxWeight[c1] = g->originalWeight[edgeIdx];
        }
        *n = c1;
        return ;
    }

    if(STaltitudes[c2] == lambda && STaltitudes[c1] < lambda ){
        attachAttributes(CompTree , STaltitudes,STAltitudesSegmentation ,  STArea, STMedia, STMaxWeight, c2 , c1, area2);

        *area1 = STArea[c1];
        // Max weight
        if( g->originalWeight[edgeIdx] > STMaxWeight[c2]){
            STMaxWeight[c2] = g->originalWeight[edgeIdx];
        }
        
        *n = c2;
        return ;
    }

    
    JCsoncell *newsoncell1;
    z = CT->nbnodes; CT->nbnodes++;
    
    nbsoncellsloc = CT->nbsoncells;
    
    CT->tabnodes[z].lastson = NULL;
    CT->tabnodes[z].sonlist = NULL;
    CT->tabnodes[z].nbsons = 0;
    
    STaltitudes[z] = lambda;
    STAltitudesSegmentation[z] = lambda-1; // For compatibilty with normal QFZ 
    
    attachAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea, STMedia, STMaxWeight, z, c1, area1);
    attachAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea, STMedia, STMaxWeight, z, c2, area2);
    /////////////// Attributes /////////////
    //Area
    *area1 = STArea[c1];
    *area2 = STArea[c2];
    //// get the maximum weight
    STMaxWeight[z] = g->originalWeight[edgeIdx];
    if(STMaxWeight[c1] > STMaxWeight[z]){
        STMaxWeight[z] = STMaxWeight[c1];
    }
    if(STMaxWeight[c2] > STMaxWeight[z]){
        STMaxWeight[z] = STMaxWeight[c2];
    }
    ///////////////////////////////////////
    

    

    
    newsoncell1 = &(CT->tabsoncells[nbsoncellsloc]);
    newsoncell1->son = z;
    newsoncell1->next =NULL;
    newsoncell1->last =NULL;
    
    CT->tabnodes[z].father = CT->root;
    CT->tabnodes[CT->tabnodes[z].father].nbsons += 1;
    CT->tabnodes[z].nodeAsSon = newsoncell1;
    
    if(CT->tabnodes[CT->root].lastson != NULL){
        CT->tabnodes[CT->tabnodes[z].father].lastson->next = newsoncell1 ;
        newsoncell1->last = CT->tabnodes[CT->tabnodes[z].father].lastson ;
        CT->tabnodes[CT->tabnodes[z].father].lastson = newsoncell1;
    }
    else{ // root was left with no children
        CT->tabnodes[CT->root].lastson = newsoncell1 ;
        CT->tabnodes[CT->root].sonlist = newsoncell1 ;
    }
    *n = z;
    CT->nbsoncells+=1;
    
    //ahora necesito darle valor a la media de n
    STMedia[*n] = (int)(STMaxWeight[*n]/(STArea[c1] + STArea[c2]));
        
}


int treeUpdateAttributes( /* INPUTS */
               JCctree **CompTree, /* Component tree - QFZ */
               double  *STaltitudes,
               double *STAltitudesSegmentation,
               int* STArea,
               double* STMedia,
               double *STMaxWeight,
               graphe *g,
               int edgeIdx /* Edge index for {x,y}  */
/* OUTPUTS */
//       CompTree
){
    
    int c1,p1,c2,p2, n, tmp;
    int p1p, p2p, np;
    int area1, area2;
    area1 = area2 = 0;
    //printf("{%i,%i} and w=%f \n" , g->tete[edgeIdx], g->queue[edgeIdx], g->weight[edgeIdx]);
    findTransition(*CompTree, g,STaltitudes , g->tete[edgeIdx] ,  edgeIdx  , &c1,&p1);
    findTransition(*CompTree, g,STaltitudes , g->queue[edgeIdx] ,  edgeIdx  , &c2,&p2);
    
    
    //printf("c1 %i c2 %i \n" ,  c1, c2);
    
    createNodeAttributes(CompTree, STaltitudes, STAltitudesSegmentation, STArea, STMedia, STMaxWeight, g, edgeIdx , c1,c2 , &n, &area1, &area2);

    while(p1 != p2){

        int tmp_p1 = p1;
        int tmp_p2 = p2;
        if(STaltitudes[p2] < STaltitudes[p1]){
            // SWAP
            tmp = p1;p1= p2; p2 = tmp;
            tmp = area2;
            area2 = area1;
            area1 = tmp;

        }
        if(STaltitudes[p1] < STaltitudes[p2]){
            p1p = (*CompTree)->tabnodes[ p1 ].father ;
            attachAttributes(CompTree ,STaltitudes,STAltitudesSegmentation, STArea, STMedia, STMaxWeight, p1 , n , &area1);
            n = p1;
            p1 = p1p;
        }
        else{
            p1p = (*CompTree)->tabnodes[ p1 ].father;  p2p = (*CompTree)->tabnodes[ p2 ].father;
            //Merging
            mergeAttributes(CompTree , STaltitudes, STArea, STMaxWeight, p1,p2, &np, &area1, &area2);
            int dummy = 0;
            attachAttributes(CompTree , STaltitudes,STAltitudesSegmentation, STArea, STMedia, STMaxWeight, np , n , &dummy);
            n = np;
            p1 = p1p;
            p2 = p2p;
        }

       
        

    }
    
   /*
    for(int i =g->nbsom ; i< 2*g->nbsom; i++ ){
        if(STaltitudes[i]!= INFINITO) STAltitudesSegmentation[i] = STaltitudes[i] - 1;        
    }*/
    
    return 0;
}


/* The algorithm is an original adaptation/improvment to edge-weighted
   graph of the algorithm of Najman and Couprie published in IEEE IP
   in 2006 */
int32_t saliencyTree2(/* INPUTS */
		     graphe *g, 
		     /* OUTPUTS */
		     JCctree ** SaliencyTree, /* Component tree */
		     double *STaltitude,  /* weights associated to the
					    nodes of the component
					    tree : Must be allocated
					    elsewhere */
            int ** pMST
		     )
{

    
  int32_t i,x1,x2,n1,n2,z,k, STx, STy, nbsoncellsloc;
  JCctree *ST;
  int32_t *clefs; 
  int32_t *STmap;
  JCsoncell * newsoncell1;
  JCsoncell * newsoncell2;
  Tarjan *T;
  int32_t taille = g->nbsom;
  int32_t narcs = g->ind;
  int nbarcsMST = taille-1;
  int nbarcsFound = 0;
  int * MST;
    if((MST = (int32_t *)malloc(sizeof(int32_t) * nbarcsMST)) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de malloc\n"); 
    exit(0);
  }
   *pMST = MST;
  if((STmap = (int32_t *)malloc(sizeof(int32_t) * taille)) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de malloc\n"); 
    exit(0);
  }
   
  if((clefs = (int32_t*)malloc(sizeof(int32_t) * narcs)) == NULL){
    fprintf(stderr,"jcSaliencyTree: erreur de malloc\n");
    exit(0);
  }
  
  for(i = 0; i < narcs; i++)
    clefs[i] = i; 
  d_TriRapideStochastique(clefs, g->weight, 0, narcs-1);
  if( (ST = componentTreeAlloc((2*taille))) == NULL){
    fprintf(stderr, "jcSalliancyTree: erreur de ComponentTreeAlloc\n");
    exit(0);
  }
  ST->nbnodes = taille;
  ST->nbsoncells = 0;
  T = CreeTarjan(taille);
  for(i = 0; i < taille; i++)
  {
    STmap[i] = i;
    ST->tabnodes[i].nbsons = 0;
    ST->tabnodes[i].sonlist = NULL;
    TarjanMakeSet(T,i);
  }

  for(i = 0; i < narcs && nbarcsFound<nbarcsMST; i++){ 
    // for each edge of the MST taken in increasing order of altitude
    n1 = TarjanFind(T, g->tete[clefs[i]]);  n2 = TarjanFind(T,g->queue[clefs[i]]);
    if(n1 != n2) {
      /* If the two vertices do not belong to a same connected
	 component at a lowest level */
      // Which component of ST n1 and n2 belongs to ?
      x1 = STmap[n1]; x2 = STmap[n2];    
      // Create a new component
      z = ST->nbnodes; ST->nbnodes++; 
      nbsoncellsloc = ST->nbsoncells;
      ST->tabnodes[z].nbsons = 2;
      // the altitude of the new component is the altitude of the edge
      // under consideration
      STaltitude[z] = getweight(g, clefs[i]);
      // add x1 and x2 to the lists of sons of the new component
      newsoncell1 = &(ST->tabsoncells[nbsoncellsloc]);  
      newsoncell2 = &(ST->tabsoncells[nbsoncellsloc+1]);
      ST->tabnodes[z].sonlist = newsoncell1;
      newsoncell1->son = x1;
      newsoncell1->next = newsoncell2;
      newsoncell2->son = x2;
      newsoncell2->next = NULL;
      ST->tabnodes[z].lastson = newsoncell2;
      ST->nbsoncells += 2
;
    
      // then link n1 and n2
      k = TarjanLink(T, n1, n2);
      STmap[k] = z;
     MST[nbarcsFound]=clefs[i];
     ++nbarcsFound;
    }//if(...)
  }//for(...)

  
  for(i = 0; i < taille; i++)
    STaltitude[i] = 0; /* altitudes of the leaves is 0 */

  /* Construct father relationship */
  ST->tabnodes[ST->nbnodes-1].father = -1;
  ST->root = ST->nbnodes - 1;
  calculReversePointer(ST, ST->root); 


  (*SaliencyTree) = ST;

  /* liberation de la memoire */
  free(STmap); 
  free(clefs);
  TarjanTermine(T);
  return ST->nbnodes - 1;
}


graphe * MSTToGraph(graphe *g, int * MST, double * weights)
{
    int nbsom =g->nbsom;
    int nbarcs = g->nbsom-1;
    graphe * g2 =initgraphe(nbsom, 2*nbarcs);
    int i;
    for(i=0; i<nbarcs; ++i)
    {
        int tete = g->tete[MST[i]];
        int queue  = g->queue[MST[i]];
        double w = weights[i+nbsom];
        addarete(g2, tete, queue, w);
    }
    return g2;
}
