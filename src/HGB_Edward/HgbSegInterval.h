
#include "include/lcomptree.h"
#include "include/list.h"



double computeD1( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double Dif){
    
    double  sc1_c2, sc2_c1, D;
    sc1_c2 = (Dif -  STMaxWeightCx )* AreaCx;
    sc2_c1 = (Dif -  STMaxWeightCy )* AreaCy;
    D = max(sc1_c2 , sc2_c1); //if(D<0)D=0; // thersholding dissimilarity measure

    return D;
}
double computeD2( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double Dif){
    
    double  sc1_c2, sc2_c1, D; //volumen
    
    sc1_c2 = (Dif -  STMaxWeightCx )* (AreaCx*STMaxWeightCx);
    sc2_c1 = (Dif -  STMaxWeightCy )* (AreaCy*STMaxWeightCy);
    D = max(sc1_c2 , sc2_c1); //if(D<0)D=0; // thersholding dissimilarity measure
    if (D<=0) D=1;
    return D;
}

double computeD3( int AreaCx, int AreaCy, double STMediaCx, double STMediaCy, double Dif){
    
    double  sc1_c2, sc2_c1, D;  
    sc1_c2 = (Dif -  STMediaCx )* AreaCx;
    sc2_c1 = (Dif -  STMediaCy )* AreaCy;
    D = max(sc1_c2 , sc2_c1); //if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}

double computeD4( int AreaCx, int AreaCy, double STMediaCx, double STMediaCy, double Dif ){
    
    double  sc1_c2, sc2_c1, D;
    sc1_c2 = (Dif -  STMediaCx )* AreaCx*STMediaCx;
    sc2_c1 = (Dif -  STMediaCy )* AreaCy*STMediaCy;
    D = max(sc1_c2 , sc2_c1); //if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}

double computeD( int AreaCx, int AreaCy, double STMaxWeightCx, double STMaxWeightCy, double STMediaCx, double STMediaCy, double Dif, int typeD = 1){
    
    double  sc1_c2, sc2_c1, D;
    
    sc1_c2 = (Dif -  STMaxWeightCx )* AreaCx;
    sc2_c1 = (Dif -  STMaxWeightCy )* AreaCy;
    
    D = max(sc1_c2 , sc2_c1);

    switch (typeD)
    {
        case 1:
            D = computeD1(AreaCx, AreaCy, STMaxWeightCx, STMaxWeightCy, Dif);   //inicial             
            break;
        case 2:
            D = computeD2(AreaCx, AreaCy, STMaxWeightCx, STMaxWeightCy, Dif); //volumen
            break;
        case 3:
            D = computeD3(AreaCx, AreaCy, STMediaCx, STMediaCy, Dif);   //media
            break;    
        case 4:
            D = computeD4(AreaCx, AreaCy, STMediaCx, STMediaCy, Dif); //volumen usando media
            break;    
        default:
            D = computeD1(AreaCx, AreaCy, STMaxWeightCx, STMaxWeightCy, Dif);  //default : inicial
            break;
    }
    
    if(D<0)D=0; // thersholding dissimilarity measure
    
    return D;
}

void filterInterval( NodeList** list, int threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list;
    ite = tmp;
    
    while(ite!=NULL){
        ite = tmp->next->next;
        if( fabs(tmp->next->value - tmp->value) < threshold){
            if( tmp->last!=NULL ){
                tmp->last->next = tmp->next->next;
                if(ite!=NULL){
                    ite->last = tmp->last;
                }
                
                free(tmp->next);
                free(tmp);
            }
            else{
                ite->last = NULL;
                *list = ite;
                free(tmp->next);
                free(tmp);
            }
        }
        tmp = ite;
    }
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
     
}


void filterInterval2( NodeList** list, int threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
//    if( (*lastNode)->value - (*list)->value <= threshold){
//        printf("size of interval less than threshold \n");
//        return;
//    }
    
    double sumLenght = 0 ;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumLenght = 0;
        while(1){
            
            sumLenght+= fabs(tmp->next->value - tmp->value);
            if(tmp->next->bound!=1)
                tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
            }
        }
        
        ite = lastp->next;
        
        
        // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumLenght < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
}

void filterIntervalVolume( NodeList** list, int threshold,  NodeList** lastNode){
    
    NodeList* ite;
    NodeList* tmp = *list, *tmp2;
    NodeList* firstp = *list;
    NodeList* lastp = *list;
    ite = tmp;
    
    double a, b, c;
    double sumVol = 0 ;
    
    while(ite!=NULL){
        
        firstp = tmp;
        sumVol = 0;
        while(1){
            c = tmp->next->value - tmp->value;
            a = fabs(tmp->dis - tmp->value);
            b = fabs(tmp->dis - tmp->next->value);
            sumVol+= c*(a+b)/2;
            if(tmp->next->bound!=1)
               tmp = tmp->next->next;
            else{
                lastp = tmp->next;
                break;
                }
        }

        ite = lastp->next;
       // printf("Volume for %.6f  %.6f is \t %.6f \n", firstp->value, lastp->value, sumVol );
        
        if( sumVol < threshold){
            if( firstp->last!=NULL ){
                firstp->last->next = lastp->next;
                if(ite!=NULL){
                    ite->last = firstp->last;
                }
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
            else{
                ite->last = NULL;
                *list = ite;
                while(firstp!=lastp){
                    tmp2 = firstp->next;
                    free(firstp);
                    firstp=tmp2;
                }
            }
        }
        tmp = ite;
    }
    
    ite = *list;
    while(ite!=NULL){
        *lastNode = ite;
        ite = ite->next;
    }
   // printf("scale %.8f\n", (*lastNode)->value  );
}

int  getNegativeIntervals(NodeList* positiveInterval, NodeList** negativeInterval, NodeList** lastNinterval, double *total_length){
    
    int counter = 0;
    //double lbound = -999999999999, ubound;
    double lbound = 0, ubound;
    double tlength = 0;
    
    NodeList* tmp = positiveInterval;
    NodeList* lastNode = NULL;
    while(tmp != NULL){
        ubound = tmp->value-0.001;
        addValueNormal( negativeInterval, &lastNode, lbound );
        addValueNormal( negativeInterval, &lastNode, ubound );
        
        tlength+=fabs(ubound-lbound);
        
        if((tmp->next)->next!= NULL)
            lbound = (tmp->next)->value + 0.001;
        tmp=(tmp->next)->next;
        counter++;
    }
    *lastNinterval = lastNode;
    *total_length = tlength;
    
    
    
    return counter;
}

double getRankNegative(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            rank_val = ite->value - obs;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
}
double getRankPositive(NodeList* list, double p){
    
    double rank_val = -1111;
    NodeList* ite = list;
    double total_l = 0, p_obs;
    
    while(ite!=NULL){
        total_l+=(  ite->next->value - ite->value );
        
        ite = ite->next->next;
    }
    //p_obs = floor( total_l*p );
    p_obs = ( total_l*p );
    
    ite = list;
    while(ite!=NULL){
        if( (  ite->next->value - ite->value ) - p_obs  >=  0 ){
            rank_val = ite->value+p_obs;
            break;
        }
        else{
            p_obs-=(  ite->next->value - ite->value );
            ite = ite->next->next;
        }
        
    }
    return rank_val;
}

double getRankPositiveBound(NodeList* list, double p){
    
    double rank_val = -1111;
    NodeList* ite = list;
    double total_l = 0, p_obs;
    
    while(ite!=NULL){
        total_l+=(  ite->next->value - ite->value );
        
        ite = ite->next->next;
    }
    //p_obs = floor( total_l*p );
    p_obs = ( total_l*p );
    
    ite = list;
    while(ite!=NULL){
        if( (  ite->next->value - ite->value ) - p_obs  >=  0 ){
            rank_val = ite->value;
            break;
        }
        else{
            p_obs-=(  ite->next->value - ite->value );
            ite = ite->next->next;
        }
        
    }
    return rank_val;
}

void printInterval(int x,int y, NodeList* list){
    NodeList* tmp = list;
    printf("Edge %i %i\n", x, y);
    while(tmp!=NULL){
      //  if(  (tmp->next->value - tmp->value) < 9000000000.000)
      //  printf("%.3f\t%.3f\t%.3f\n", tmp->value, tmp->next->value , (tmp->next->value - tmp->value));
        printf("%.3f\t%.3f \n", tmp->value, tmp->next->value  );
        
        tmp = tmp->next->next;
    }
}

double getRankNegativeBound(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            
            if(ite->last->value!=0)
                //rank_val = ite->last->value ;
                rank_val = ite->value ;
            else
                rank_val = ite->value;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
}

double getRankNegativeVolume(NodeList* list, NodeList* lastNode, double obs){
    
    double rank_val = -1111;
    
    NodeList* ite;
    NodeList* tmp;
    
    ite = list;
    ite->value = 0 ;
    double length = 0;
    
    while(ite!=NULL){
        length+= fabs( ite->next->value - ite->value );
        ite=ite->next->next;
    }
    
    obs = floor( length*obs );
    
    
    tmp = lastNode;
    ite = tmp;
    
    while(ite!=NULL){
        if( ( (ite->value - ite->last->value) - obs ) >  0 ){
            rank_val = ite->value - obs;
            break;
        }
        else{
            obs-=(ite->value - ite->last->value);
            ite = ite->last->last;
        }
    }
    return rank_val;
}

double minimizationLambdasInterval(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double *STMedia, double* STMaxWeight, int edgeIdx , int op, int* numberInterval, int *max_numberInterval, int typeComputeD = 1)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    
    double lambdaUpper = -9999;
    double lambdaLower = -9999;
    double total_l;
    NodeList* positiveInterval = NULL;
    NodeList* negativeInterval = NULL;    
    NodeList* LastNode = positiveInterval;
    NodeList* lastNinterval = NULL;
    int cn_intervals = -1;
    
    cx = x;
    cy = y;
    int T = -1;
    double p_obs=0;
    lambda = STAltitudes[cx];
    double prevLambda=-1, lambda_star =-9999, l_1,l_2;
    int cinterval=0;
    
      //  printf("Edge %i,%i \n",x,y );
    do{
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
        
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
          //  printf( "%.2f,%i,%i,%.2f,%.2f,%.2f\n " ,  lambda ,STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy]  ,  D );
            if(STAltitudes[px ] == lambda )
                cx =px;
            if(STAltitudes[py ] == lambda )
                cy =py;
        }
        if(D <= prevLambda) lambdaLower = prevLambda + 0.0001;
        else
            lambdaLower = D;
        
        double prevD = D;
        
        while( D<=lambda  &&( px!=-1&&py!=-1)){
            prevLambda = lambda;
            
            if(px!=-1&&py!=-1)
                lambda = min(STAltitudes[px ] , STAltitudes[py] );
            else{
               if(py!=-1&& px == -1 )
                   lambda = STAltitudes[py];
                else
                    if(px!=-1&& py == -1 )
                        lambda = STAltitudes[px];
                    else
                        lambda =STAltitudes[cx]; // arrived to root
            }
            
            prevD = D;
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
            
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            
            if(px!=-1 && STAltitudes[px ] == lambda )
                cx =px;
            if(py!=-1 && STAltitudes[py ] == lambda )
                cy =py;
            
            
        }
        if(prevD < prevLambda )
            lambdaUpper = prevLambda + 0.0001;
        else lambdaUpper = prevLambda;
        
        addValueNormal(&positiveInterval, &LastNode, lambdaLower);
        addValueNormal(&positiveInterval, &LastNode, lambdaUpper);
       // if(cinterval>=1)
       //    printf("  bounds <%.3f , %.3f > \n", lambdaLower, lambdaUpper );
        *numberInterval+=1; // counting number of intervals
        cinterval+=1;
    }
    while(px!=-1 && py!=-1);
    
    //if(cinterval>=1){
        //printInterval(x,y,positiveInterval);
        //printf(" intervals %i\n", cinterval);
        //}
    switch(op){
        case 1: //Min: min value on the first observation interval
            lambda_star  = positiveInterval->value;
            break;
        case 2: // Max, first value on the last observation interval
            lambda_star  = LastNode->last->value;
            
                            //printInterval(x,y,positiveInterval);
                            //printf("lambda %.3f\n", lambda_star);
            
            break;
        case 3:  // min lowet bound of filtered positive intervals
            T=10;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 4:
            T=100;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            break;
        case 5:
            T=500;
            filterInterval(&positiveInterval,T,&LastNode  );
            lambda_star  = positiveInterval->value;
            
            break;
        case 6:
            T=10;
            
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;            
            break;
        case 7:
            T=100;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 8:
            T=500;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 9:
            T=700;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
        case 10:
            T=900;
            getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            filterInterval(&negativeInterval,T,&lastNinterval  );
            lambda_star  = lastNinterval->value;
            break;
            
        case 11:
            //P = 0.3
            //P = 0.1
            //P = 0.01
            //P = 0.05
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.01 );
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
//            if( cinterval >3){
//            printInterval(x,y,negativeInterval);
//            printf("lambda %.3f\n",lambda_star );
//            }
            
            break;
        case 12:
            //P = 0.3;
            //P = 0.1
            //P=0.01
            //P = 0.005
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.005 );
           
            //for boundaries
            //lambda_star = getRankNegativeBound(negativeInterval, lastNinterval,p_obs);
            //For rank value only
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
           // printf("lambda %.3f \n", lambda_star);
            
            break;
            
        case 13:
            // 0.01 on positive
            //p_obs = floor( total_l*0.01 );
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.9);
            break;
        case 14:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.99);
            break;
            
        case 15:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.01);
            break;
        
        case 16:
            LastNode->value = LastNode->last->value + 1;
            lambda_star = getRankPositive(positiveInterval, 0.001);
//            if(cinterval>1){
//                printInterval(x,y,positiveInterval);
//                printf("lambda %.3f\n", lambda_star);
//            }
            break;
        
        case 17:
            
            
            LastNode->value = LastNode->last->value + 10;
            lambda_star = getRankPositiveBound(positiveInterval, 0.99999);
            
            
                            printInterval(x,y,positiveInterval);
                            printf("lambda %.3f\n", lambda_star);
            
            
            break;
            
            
        case 20:
            
        //Silvio's :
            //identify negative intervals
           // getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            //lambda_star = lastNinterval->value;
           //lambda_star  = positiveInterval->value;
            cn_intervals = getNegativeIntervals(positiveInterval, &negativeInterval, &lastNinterval, &total_l);
            p_obs = floor( total_l*0.01 );
            lambda_star = getRankNegative(negativeInterval, lastNinterval,p_obs);
            
            //if(cinterval>1)
               // printInterval(x,y,negativeInterval);
         //  printf("lambda %.3f \n", lambda_star);
            break;
        case 21:
            
            lambda_star  = LastNode->last->value;
            break;
    }
    
    if(cinterval>  *max_numberInterval)
    {
        *max_numberInterval = cinterval;
        
    }
   // printf("  ->  edge %i,%i   lambda %.3f \n",x,y, lambda_star);
    
    while(positiveInterval!= NULL){
        LastNode = positiveInterval;
        positiveInterval = positiveInterval->next;
        free(LastNode);
    }
    
    return lambda_star;
}


double scaleSelection(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes, int *STArea, double *STMedia, double* STMaxWeight, int edgeIdx , int op, int th, double prank, int typeComputeD = 1)
{
    JCctree* CT = *CompTree;
    int cx, cy, px, py;
    px = py = -9999;
    double  D=0, lambda;
    
    int x = mst->MSTedges[edgeIdx].x;
    int y = mst->MSTedges[edgeIdx].y;
    
    double Dif=  mst->MSTedges[edgeIdx].weight;
    
    
    double lambdaUpper = -9999;
    double lambdaLower = -9999;
   // double total_l;
   // NodeList* positiveInterval = NULL;
    NodeList* positiveIntervalVol = NULL;
    NodeList* negativeIntervalVol = NULL;
    
    //NodeList* negativeInterval = NULL;
    //NodeList* LastNode = positiveInterval;
    NodeList* ln_neg = negativeIntervalVol;
    NodeList* ln_pos = positiveIntervalVol;
  //  NodeList* lastNinterval = NULL;
  //  int cn_intervals = -1;
    
    cx = x;
    cy = y;
    int T = -1;
    double p_obs=0;
    
    
    lambda = STAltitudes[cx];
    double prevLambda=-99999999, lambda_star =-9999;
    int cinterval=0;
    double prevD = -1;
    
    
   // printf("Edge %i,%i   intervals \n",x,y );
    do{
        //D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], Dif );
        D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
        //cout<<endl<<"despues del primer computeD "<<endl;
        //cout<<"D : "<<D<<", lambda : "<<lambda <<endl;

        prevD = -1;
        prevLambda = lambda;
        while( D >lambda ){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            prevLambda = lambda;
            lambda = min(STAltitudes[px ] , STAltitudes[py] );
            //D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
            if(D>lambda  ){
                addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
                addValueDis(&negativeIntervalVol, &ln_neg, lambda, D);
                prevD = D;
                //prevLambda = lambda;
                if(STAltitudes[px ] == lambda )
                    cx =px;
                if(STAltitudes[py ] == lambda )
                    cy =py;
            }
        }
        if(D <= prevLambda) {
            lambdaLower = prevLambda + 0.009;
            //lambdaLower = prevLambda ;
        }
        else{
            lambdaLower = D ;
            addValueDis(&negativeIntervalVol, &ln_neg, prevLambda, D);
            addValueDis(&negativeIntervalVol, &ln_neg, lambdaLower, D);
        }
        
        ln_neg->bound = 1;
        prevLambda = lambdaLower;
        prevD = -99999999;
        
//        if(lambda == prevLambda ){
//            if(px!=-1 && STAltitudes[px ] == lambda )
//                cx =px;
//            if(py!=-1 && STAltitudes[py ] == lambda )
//                cy =py;
//        }
        
        
        while( D<=lambda  &&( px!=-1&&py!=-1)){
            px = CT->tabnodes[cx].father;
            py = CT->tabnodes[cy].father;
            
            if(px!=-1&&py!=-1)
                lambda = min(STAltitudes[px ] , STAltitudes[py] );
            else{
                if(py!=-1&& px == -1 )
                    lambda = STAltitudes[py];
                else
                    if(px!=-1&& py == -1 )
                        lambda = STAltitudes[px];
                    else{
                         break;
                    }
            }
            
            //D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy] , Dif);
            D = computeD(STArea[cx], STArea[cy] ,STMaxWeight[cx], STMaxWeight[cy], STMedia[cx], STMedia[cy], Dif, typeComputeD);
            if(D<=lambda  ){
                addValueDis(&positiveIntervalVol, &ln_pos, prevLambda, D);
                addValueDis(&positiveIntervalVol, &ln_pos, lambda, D);
                prevD = D;
                prevLambda = lambda;
                if(px!=-1 && STAltitudes[px ] == lambda )
                    cx =px;
                if(py!=-1 && STAltitudes[py ] == lambda )
                     cy =py;
            }
//            else if(prevLambda == lambda ){
//                    if(px!=-1 && STAltitudes[px ] == lambda )
//                        cx =px;
//                    if(py!=-1 && STAltitudes[py ] == lambda )
//                        cy =py;
//                 }
            

        }
        if(prevD < prevLambda )
            lambdaUpper = prevLambda + 0.009;
            //lambdaUpper = prevLambda ;
        else lambdaUpper = prevLambda;
        
        ln_pos->bound = 1;
        
        //addValueNormal(&positiveInterval, &LastNode, lambdaLower);
        //addValueNormal(&positiveInterval, &LastNode, lambdaUpper);
        
        cinterval+=1;
    }
    while(px!=-1 && py!=-1);

    switch(op){
        case 1: //Min: min value on positive observation intervals
            
            lambda_star  = positiveIntervalVol->value;
            
            
            break;
        case 2: // Max,  last upper bound on negative intervals
            
//            printf("\t Neg \n");
//            printIntervalDis(x,y,negativeIntervalVol);
//            printf("\t Pos \n");
//            printIntervalDis(x,y,positiveIntervalVol);
//            printf("\t All positive intervals \n" );
//            printInterval(x,y,positiveInterval );
//
//            printf("lambda last note, %.3f \t %.3f  \n ", LastNode->last->value , ln_neg->value);
            lambda_star = ln_neg->value;
            
            break;
        case 3:  // On positive intervals, apply length thresholding and  min-rule
            //T=10;
            T = th;
            filterInterval2(&positiveIntervalVol,T,&ln_pos  );
            lambda_star  = positiveIntervalVol->value;
            
           // printf("lambda star note, %.3f   \n ", lambda_star);
            break;
            
        case 4:
            // On negative intervals, apply length thresholding and  max-rule
            T=th;            
           // printf("\t Neg \n");
            //printIntervalDis(x,y,negativeIntervalVol);
            
            negativeIntervalVol->value = -9999999999;  // force negative to have a large interval by the right (to deal with case of only one negative interval)
            
            filterInterval2(&negativeIntervalVol,T,&ln_neg  );
            
            //printf("\t Neg Filtered\n");
            //printIntervalDis(x,y,negativeIntervalVol);
            
            lambda_star  = ln_neg->value;
            break;
        case 5:
            // On positive intervals, apply rank and  min-rule
            ln_pos->value =ln_pos->last->value+1; // limit the upper bound 
            lambda_star = getRankPositive(positiveIntervalVol, prank);
            break;
        case 6:
            // On negative intervals, apply rank and  max-rule
            lambda_star = getRankNegativeVolume(negativeIntervalVol, ln_neg,prank);
            break;
            
        case 7:
            
            // On positive intervals, apply area and  min-rule
            T = th;
            filterIntervalVolume(&positiveIntervalVol, T, &ln_pos);
            lambda_star =positiveIntervalVol->value;
            
            break;
        case 8:
            // On negative intervals, apply area and  max-rule
            T = th;
            negativeIntervalVol->value = -99999999;
            filterIntervalVolume(&negativeIntervalVol, T, &ln_neg);
            lambda_star = ln_neg->value;
            break;
        case 9:
             // On negative intervals, apply area, then ranking and  max-rule
            T = th;
            negativeIntervalVol->value = -99999999;
            filterIntervalVolume(&negativeIntervalVol, T, &ln_neg);
            p_obs = prank;
            lambda_star = getRankNegativeVolume(negativeIntervalVol, ln_neg,p_obs);
            break;
        case 10:
            // On positive intervals, get lowerbound on last interval
            
       //     printIntervalDis(x,y,positiveIntervalVol);
            lambda_star = ln_pos->last->value ;
          //  printf("lambda %.3f\n",lambda_star );
            
            break;
    }
    
    
//    while(positiveInterval!= NULL){
//        LastNode = positiveInterval;
//        positiveInterval = positiveInterval->next;
//        free(LastNode);
//    }

    while(positiveIntervalVol!= NULL){
        ln_pos = positiveIntervalVol;
        positiveIntervalVol = positiveIntervalVol->next;
        free(ln_pos);
    }
    
    while(negativeIntervalVol!= NULL){
        ln_neg = negativeIntervalVol;
        negativeIntervalVol = negativeIntervalVol->next;
        free(ln_neg);
    }
    
    return lambda_star;
}


void reweightToMax(graphe *g, MST* mst, JCctree **CompTree,  double  *STAltitudes)
{
    JCctree* CT = *CompTree;
    
    
    double   lambda;
    int i;
    
     for(i=0;i<mst->size;i++){
        
         
        int cx, cy, px, py;
        int x = mst->MSTedges[i].x;
        int y = mst->MSTedges[i].y;
         
         cx = x;
         cy = y;
        lambda = STAltitudes[cx];
         double newf = -9999;
            while(  cx!=cy ){ // search common parent
                px = CT->tabnodes[cx].father;
                py = CT->tabnodes[cy].father;
                
                if(px!=-1&&py!=-1)
                    lambda = min(STAltitudes[px ] , STAltitudes[py] );
                else{
                    if(py!=-1&& px == -1 )
                        lambda = STAltitudes[py];
                    else  if(px!=-1&& py == -1 )
                            lambda = STAltitudes[px];
                          else
                            lambda =STAltitudes[cx]; // arrived to root
                }
                px = CT->tabnodes[cx].father;
                py = CT->tabnodes[cy].father;
                if(px!=-1 && STAltitudes[px ] == lambda )
                    cx =px;
                if(py!=-1 && STAltitudes[py ] == lambda )
                    cy =py;
            }
         
         px = CT->tabnodes[cx].father;
         
         if(px == CT->root || px==-1 ){
             newf= STAltitudes[cx]+1;
         }
         else{
             
             newf= STAltitudes[px];
         }
         
         setweight(g,i,newf);
         
           //  printf("edge %i,%i   lambda %.3f \n",x,y, newf);
         
         
     }
    
    
}



void incrementalSegmentationInterval(graphe *g, int op){
    
    
  //  clock_t startAll, endALL; double cpu_time_usedALL;

    int i; double  newf;
    MST *mst; graphe *newG;
    
   
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double* STMedia = NULL; // STMedia
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    STMedia     = (double *)calloc(2*g->nbsom, sizeof(double)); // STMedia
    
    
    //////////////////  Incremental Segmentation  //////////////
    
    
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMedia , STMaxWeight );
    
    
    /// Set all edges to large values for later SM computation ////
    
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    //startAll = clock();
    
    int numberInterval=0;
    int max_numberInterval=0;
    
    for(i=0;i<mst->size;i++){
        NodeList* lambdaStarInterval;
        ListInit(&lambdaStarInterval);
        
       
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 1 ;
        newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMedia , STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;

        setweight(newG,i,newf);

        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMedia, STMaxWeight, newG, i);

    }
    

    //cout << "finish for insid mst->size : "<< endl;
    //endALL = clock();
    //cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
    
    //printf("\t mean number of intervals %i, max number of intervals  %i\n", numberInterval/mst->size, max_numberInterval );

    //cout<< "mean number of intervals : "<< numberInterval/mst->size<<  ", max number of intervals :" << max_numberInterval<<endl;
    
    
    if(op>=20){
        //Set values to Max on last interval
        reweightToMax(newG, mst,  &CompTree ,STAltitudesSegmentation);
    }
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
   
    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);
}

void intervalSegmentation(graphe *g, int op, int th, double prank, int typeComputeD){
    
    
    //  clock_t startAll, endALL; double cpu_time_usedALL;
    int i; double  newf;
    MST *mst; graphe *newG;
    
    getMST(g, &mst); // The MST is computed and stored in mst
    graphFromMST(&newG, mst,g->nbsom ,mst->size, INFINITO); // Edges on MST are already sorted.
    
    /////// Incremental QFZ structures////////////
    JCctree* CompTree = NULL; // Component tree - QFZ
    double* STAltitudes;
    double * STAltitudesSegmentation;
    int* STArea = NULL;
    double* STMedia = NULL; // STMedia
    double * STMaxWeight = NULL;
    
    if( (CompTree = componentTreeAlloc((g->nbsom*2))) == NULL){
        fprintf(stderr, "CT: allocation error\n");
        exit(0);
    }
    
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudes == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudes\n");
        exit(0);
    }
    STAltitudesSegmentation = (double *)calloc(2*g->nbsom, sizeof(double));
    if (STAltitudesSegmentation == NULL){
        fprintf(stderr," cannot allocate memory for STAltitudesSegmentation\n");
        exit(0);
    }
    
    STMaxWeight = (double *)calloc(2*g->nbsom, sizeof(double));
    STArea      = (int *)calloc(2*g->nbsom, sizeof(int));
    STMedia     = (double *)calloc(2*g->nbsom, sizeof(double)); // STMedia
    
    
    //////////////////  Incremental Segmentation  //////////////
    
    //cout<<"antes de initializeTreeAttributes "<< endl;
    initializeTreeAttributes(newG,&CompTree,STAltitudes,STAltitudesSegmentation, STArea, STMedia, STMaxWeight );
    
    
    /// Set all edges to large values for later SM computation ////
    
    
    for(i=0;i < g->nbmaxar; i++){
        setweight(g, i, INFINITO);
    }
    //startAll = clock();
    
    //cout<< endl << "despues de setear pesos --setweight "<< endl;
    int numberInterval=0;
    int max_numberInterval=0;
    
    for(i=0;i<mst->size;i++){
        //NodeList* lambdaStarInterval;
        //ListInit(&lambdaStarInterval);
        
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 1 ;
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;
        //newf =  scaleSelection(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMaxWeight, i ,  op, th, prank) - 0.000001 ;
        //newf =  minimizationLambdasInterval(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMedia , STMaxWeight, i ,  op, &numberInterval, &max_numberInterval) - 0.0001 ;
        
        //newf =  scaleSelection(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMedia, STMaxWeight, i ,  op, th, prank, typeComputeD) - 1  ;
        newf =  scaleSelection(newG, mst,  &CompTree ,STAltitudesSegmentation , STArea, STMedia, STMaxWeight, i ,  op, th, prank, typeComputeD) - 1; //0.0001;
        
        setweight(newG,i,newf);
        treeUpdateAttributes(&CompTree, STAltitudes,STAltitudesSegmentation,  STArea, STMedia, STMaxWeight, newG, i);
    }
    //cout << "despues del setweight"<<endl;
    //endALL = clock();
    //cpu_time_usedALL = ((double) (endALL - startAll)) / CLOCKS_PER_SEC;
    
    //printf("\t mean number of intervals %i, max number of intervals  %i\n", numberInterval/mst->size, max_numberInterval );
    
    //    if(op>=20){
    //        //Set values to Max on last interval
    //        reweightToMax(newG, mst,  &CompTree ,STAltitudesSegmentation);
    //    }
    
    //////////////////////////////////////////////////////////////
    
    for(i=0;i<mst->size;i++)
        setweight(g, mst->MSTedges[i].idx, getweight(newG, i) + 1);
    

    terminegraphe(newG);
    free(STAltitudes);
    free(STMaxWeight);
    free(STAltitudesSegmentation);
    free(STArea);
    free(mst->MSTedges);
    free(mst);
    
    componentTreeFree(CompTree);

}
