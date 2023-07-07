#ifndef MST_h
#define MST_h

//#include <stdlib.h>
//#include <sys/types.h>
//#include <stdio.h>

#define INFINITO 9999999999
#define MSTConstant 0

typedef struct{
    int x, y;
    int idx;
    double weight;
}edge;


typedef edge* edges;

typedef struct{
    edges MSTedges;
    int size;
    
}MST;



void setSizeMST(MST* mst, int n){
    
    
    if(( mst->MSTedges = (edge*)malloc(sizeof(edge) * n)) == NULL){
        fprintf(stderr,"setSizeMST: allocation error\n");
        exit(0);
    }
    mst->size=0;

}


#endif /* MST_h */
