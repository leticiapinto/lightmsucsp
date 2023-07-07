
#include "functions.h"

using namespace cv;
using namespace std;



int getEdge(int node)
{   
    return node; 
}

int weight_node(vector <double>  &altitude, int node)
{
    return altitude[getEdge(node)];
}

void canonizeQBT(vector< Element* > &qbt_parents, vector< Element* > &qct_parents, vector <double>  &altitude, int vertex_number, int QCT_size)
{
    //cout<<endl<<vertex_number<<endl;
    //for all nodes n of Q BT do Q CT .parent[n]:=Q BT .parent[n]; Q CT .size+=1;
    QCT_size = 0;
    //cout <<"qbt_parents : "<<qbt_parents.size()<<endl;
    for (int i = 0; i < qbt_parents.size(); i++) 
    {
        qct_parents.push_back(qbt_parents[i]);
        QCT_size++;
    }
    //cout <<"qct_parents : "<<qct_parents.size()<<endl<< "altitude : " <<altitude.size()<<endl;
    //print_vector(altitude);
    //for each non-leaf and non-root node n of Q BT by decreasing order do
    for(int i = qbt_parents.size()-2; i >= vertex_number; i--)
    {
        ////cout<<"node i : " <<i;
        if(qbt_parents[i]->parent!=NULL && qbt_parents[i]->children.size()>0)
        {
            int p = qct_parents[i]->parent->index;
            ////cout<<"; weight_node("<<p<<") : "<<weight_node(p)<<", weight_node("<<i<<") : "<<weight_node(i)<<endl;
            if(weight_node(altitude, p) == weight_node(altitude, i))
            {
                for(int j =0; j < qct_parents[i]->children.size(); j++)
                {
                    //el padre de estos hijos es p
                    int tmp_child = qct_parents[i]->children[j]->index;
                    ////cout<<"tmp_child : "<<tmp_child<<endl;
                    //actualizar hijos de p
                    qct_parents[p]->addChild(qct_parents[tmp_child]);    //le agrego este hijo al nuevo padre
                    qct_parents[tmp_child]->parent = qct_parents[p]; //cambio el padre de estos hijos
                    //qct_parents[i]->removeChild(qct_parents[tmp_child]); //le quito los hijos a este nodo que se elimina
                }
                qct_parents[p]->removeChild(qct_parents[i]); //quito a este nodo de su padre 
                qct_parents[i]->removeAllChildren();
                qct_parents[i]->parent= NULL;
                qct_parents[i]->in_tree = 0;
                ////cout<<" nodo i"<<i<<endl;  
            }
        }
    }

    //cout << "each non-leaf and non-root node n " <<endl;
    //cout<< "qct_parents : "<<qct_parents.size()<<endl;
    for(int i = 0; i < qct_parents.size(); i++) // If needed, build the list of children
    {   
        ////cout<<"i :"<< i<<endl;
        ////cout<<"qct_parents[i]->parent : "<<qct_parents[i]->parent<<endl;
        if(qct_parents[i]->parent!= NULL)
        {
            int p = qct_parents[i]->parent->index; ////cout<< qct_parents[i]->parent->index<<endl;
            qct_parents[p]->addChild(qct_parents[i]);
        }
    }

    //cout<<endl<<"---- size "<<qct_parents.size()<< endl ;
}


int level_node(vector< Element* > &parents, int q )
{
    int i = 0;
    while(parents[q]->parent!=NULL)
    {
        q= parents[q]->parent->index;
        i++;
    }
    return i;

}


void iterate_pre(vector< Element* > &parents, vector<int> &euler, vector<int> &level, int i){

    if( parents[i]->parent==NULL && parents[i]->children.size()==0 )
    {
            ////cout<<"parents[i] "<<parents[i]->index;
            return;
    }

    ////cout<<"i" <<i <<": "<<parents[i]->in_tree<<endl;
    ////cout << parents[i]->index<<  " (" << level_node(parents, parents[i]->index)<< "), "; // do something with<fre the value 
    euler.push_back(parents[i]->index); 
    level.push_back(level_node(parents, parents[i]->index));
    
    for(int j = 0 ; j < parents[i]->children.size(); j++)
    {
        iterate_pre(parents, euler, level, parents[i]->children[j]->index);
        if(parents[i]->children[j]->parent!= NULL)
        {
            ////cout << parents[i]->children[j]->parent->index<<  " (" << level_node(parents, parents[i]->children[j]->parent->index)<< "), ";
            euler.push_back(parents[i]->children[j]->parent->index); 
            level.push_back(level_node(parents, parents[i]->children[j]->parent->index));
        }
    }
}


int LCA(int u, int v, vector<int > &R, vector<int > euler, vector<int > level )
{
    int index_u = R[u]; int index_v = R[v]; int min_,max_;
    if(index_u<= index_v)
    {
        min_ = index_u; max_ = index_v;
    }else{
        max_ = index_u; min_ = index_v;
    }

    int min_l = level[min_]; 
    int lca = euler[min_];
    for(int i=min_; i<=max_;i++)
    {
        if(level[i]< min_l)
        {
            min_l= level[i];
            lca = euler[i];
        }
    }

    return lca;
}

bool isleaf(const vector< Element* > &parents, vector<int> &level_down, int node)
{

    if(parents[node]->children.size()==0  && level_down[node]==0) return true;
    return false;
}


void fathertobottom(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node)
{
    //father to bottom
    leaves[node] = 0;
    int max_level_children = 0;     //int size_children = parents[i]->children.size();
    for(int j = 0; j < parents[node]->children.size(); j++)
    {
        int tmp_child = parents[node]->children[j]->index;     
        ////cout<<"tmp_child : "<<tmp_child<< ", leaves : "<< leaves[tmp_child]<<endl;
        if(level_down[tmp_child]!=-1)
        {
            if(isleaf(parents, level_down, tmp_child)){
                leaves[node] = leaves[node] + 1;
            }else{
                leaves[node] += leaves[tmp_child];
            }
            //level
            if(level_down[tmp_child] >= max_level_children)
            {
                max_level_children = level_down[tmp_child];
            }
            level_down[node] = max_level_children+1;
        }
    }
    
}

void fathertotop(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int father, int root)
{
    //father to top
    Element* p = parents[father];
    while (p->index != root)
    {
        Element* actual = p;
        //leaves[actual->index] -= 1;
        fathertobottom(parents, leaves, level_down, p->index);
        p = p->parent;
        /*
        if(actual->children.size()==0) // si no tiene hijos  
        {
            actual->parent = NULL;
            level_down[actual->index] = -1;
            int max_level = -2;
            for(int i = 0; i < p->children.size(); i++)
            {
                if(level_down[p->children[i]->index] > max_level ) {max_level = level_down[p->children[i]->index];}
            }
            level_down[p->index] = max_level + 1; 
        }
        */
       //fathertobottom();
       
    }
}

bool reducenode(vector< Element* > &parents, int node)
{
    if(parents[node]->children.size() == 1) return true;
    return false;
}

void removefromfather(vector< Element* > &parents, int node, int father )
{
    vector<Element* > aux_children;
    for (int j = parents[father]->children.size()- 1; j >= 0; j--)
    {
        Element* actualnode = parents[father]->children[j];
        parents[father]->children.pop_back();
        if(actualnode->index!= node ) {aux_children.push_back(actualnode);}
        else {break;}
    }

    for(int j = 0; j < aux_children.size(); j++)
    {
        parents[father]->children.push_back(aux_children[j]);
    }
}

void updatenodeinformation(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node, int root)
{
    //debemos preguntar si es necesario reducir el arbol, es decir quitar este nodo y unir el padre(node) - hijo(node)
    //reducenode(parents, node)? //cout << "True": //cout << "False";
    ////cout << endl;

    if(reducenode(parents, node) && !isleaf(parents, level_down, parents[node]->children[0]->index))
    {
        ////cout <<"reduce this node "<< endl;
        int father_node = parents[node]->parent->index;
        int son_node = parents[node]->children[0]->index;

        //likeando padre-hijo
        parents[son_node]->parent = parents[father_node];
        removefromfather(parents, node, father_node);
        parents[father_node]->addChild(parents[son_node]);

        //clean node
        parents[node]->parent = NULL;
        parents[node]->removeAllChildren();
        leaves[node] = -1;
        level_down[node] = -1;

        //el nuevo node es el padre de node
        node = father_node;
    }
    ////cout << "fathertobottom:  "<< endl;
    fathertobottom(parents, leaves, level_down, node);
    ////cout << "fathertotop: " << endl;
    fathertotop(parents, leaves, level_down, node, root);
}

void count_leaves(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int leaves_pixels, int root)
{
    //cout << "Iniciando conteo de hojas por nodo en el arbol"<<endl;
    ////cout<<"parents.size()"<< parents.size() <<endl;
    //cout <<"leaves_pixels: "<< leaves_pixels << endl;
    for(int i = leaves_pixels; i < parents.size(); i++)
    {
        if(level_down[i]!=-1)
        {
            //updatenodeinformation(parents, leaves, level_down, i, root);
            fathertobottom(parents, leaves, level_down, i);
            /*
            leaves[i] = 0;
            int max_level_children = 0;     //int size_children = parents[i]->children.size();
            for(int j = 0; j < parents[i]->children.size(); j++)
            {
                int tmp_child = parents[i]->children[j]->index;     
                ////cout<<"tmp_child : "<<tmp_child<< ", leaves : "<< leaves[tmp_child]<<endl;
                if(level_down[tmp_child]!=-1)
                {
                    if(isleaf(parents, level_down, tmp_child)){
                        leaves[i] = leaves[i] + 1;
                    }else{
                        leaves[i] += leaves[tmp_child];
                    }
                    //level
                    if(level_down[tmp_child] >= max_level_children)
                    {
                        max_level_children = level_down[tmp_child];
                    }
                    level_down[i] = max_level_children+1;
                }
            }
            */
        }
    }
    ////cout<<"leaves "<<endl;    //print_vector(leaves);    //print_vector(level_down);
}

void childrenFilter(vector< Element* > &parents, vector<int> &leaves, int leaves_pixels, int parameter_number)
{
    //print_vector(leaves); //ahora se exactamente cuantos nodo hoja tiene cada nodo
    //cout <<"Elimino con parametro : "<<parameter_number << " hasta parametro 1" << endl;
    for(int i = parameter_number; i > 0; i--)
    {
        ////cout<<"elimino con parametro "<<i<<endl;
        for(int j = leaves_pixels; j < leaves.size(); j++)
        {
            if(leaves[j] == i && parents[j]->children.size()>0)  //busco nodos con la cantidad de hojas igual al parametro i
            {
                ////cout <<"j : "<<j<<", parent : "<< parents[j]->parent->index <<endl;
                //eliminar node j
                int p = parents[j]->parent->index; //p padre del node j 
                int size_children = parents[j]->children.size(); //cantidad de hijos
                ////cout<<"size_children : "<<size_children<<endl;
                for(int k =0; k < size_children; k++)
                {
                    //el padre de estos hijos va a ser p 
                    int tmp_child = parents[j]->children[k]->index;
                    ////cout<<"tmp_child del que voy agredar : "<<tmp_child<<endl;
                    parents[p]->addChild(parents[tmp_child]);    //le agrego este hijo al nuevo padre
                    parents[tmp_child]->parent = parents[p]; //cambio el padre de estos hijos
                    //parents[j]->removeChild(parents[tmp_child]); //le quito los hijos a este nodo que se elimina
                }
                
                parents[p]->removeChild(parents[j]); //quito a este nodo de su padre
                parents[j]->removeAllChildren(); 
                parents[j]->parent= NULL;
                parents[i]->in_tree = 0;
                leaves[j] = -1;
            }
        }
        //print_vector(leaves);
    }

    //despues eliminar nodos hoja que su padre sea 

    int root = parents[parents.size()-1]->index;
    //cout << "root : "<< root << endl;
    
    for(int i = 0; i < leaves_pixels; i++)
    {
        if(parents[i]->parent->index == root) //nodo hoja que necesito eliminar
        {
            //cout<<"nodo hoja a eliminar : "<< parents[i]->index << endl;
            parents[root]->removeChild(parents[i]); //quitarle los hijos al root            
            parents[i]->parent = NULL; //quitarle el padre
        }
    }
    
}

void iterate_cut(vector< Element* > &parents, int i, int evaluate_node, int color, vector<int> &myflag, vector<int> &vector_color, int leaves_pixels){
    
    if( parents[i]->parent==NULL && parents[i]->children.size()==0 )
    {
            ////cout<<"parents[i] "<<parents[i]->index;
            return;
    }
    ////cout << parents[i]->index<<", "; // do something with the value 
    if(parents[i]->index < leaves_pixels && !myflag[parents[i]->index])
    {
        vector_color[parents[i]->index] = color; 
    }
    myflag[parents[i]->index] = 1;
    
    for(int j = 0 ; j < parents[i]->children.size(); j++)
    {
        
        if(parents[i]->children[j]->parent!= NULL)
        {
            iterate_cut(parents, parents[i]->children[j]->index, evaluate_node, color, myflag, vector_color, leaves_pixels);
            ////cout << parents[i]->children[j]->parent->index<<", ";
            if(parents[i]->children[j]->parent->index < leaves_pixels && parents[i]->children[j]->parent->index)
            {
                vector_color[parents[i]->children[j]->parent->index] = color;
            }
            myflag[parents[i]->children[j]->parent->index] = 1;
        }
        
    }
}





void iterate_coloring(vector< Element* > &parents, vector<int> level_down, vector<int> vector_color, int node, int color){

    if(isleaf(parents, level_down,  node))
    {
        vector_color[node] = color;
        
    }else{

        for(int i = 0; i < parents[node]->children.size(); i++)
        {   
            iterate_coloring(parents, level_down, vector_color, parents[node]->children[i]->index, color);
        }
    }
    
}


void ColorRegions(vector< Element* > &parents, vector<int> &R, vector<int> &vector_color, vector<int> &level, int max, int cut_tree, int leaves_pixels, bool option)
{
    //cout<<"#######  coloring regions #####"<<endl;
    //Empiezo a colorear
    //usar flag si el nodo ha sido visitado o no 
    //primero coloreo todas las hojas que tengan como hierarchy <= parametro
    //luego coloreo las hojas que no han sido coloreadas agregando cada hoja aislada a un color diferente

    //int param = 1;    //int max = *max_element(level.begin(), level.end());
    //cout<<"Iniciando coloreado de regiones "<<endl;
    int color =0;
    vector<int> myflag(parents.size(),0);

    //tomar regiones con nivel igual o menor del parametro establecido
    //debo empezar desde los nodos padres creados .. no considerar las hojas 
    //partir del numero de hojas 

    //cout << "cut_tree: "<< cut_tree<< endl;
        
        for(int j = cut_tree; j > 0; j--)
        {
            for(int i = leaves_pixels; i < parents.size(); i++)
            {
                ////cout<<"j : "<<j<<", i : "<<i<<", level["<<i<<"] = "<<level[i]<<"\tmyflag[i] "<<myflag[i]<<"\t";
                if((level[i] == j && !myflag[i] && (parents[i]->children.size()>0) ))
                {
                    ////cout<<"\tcolor = " <<color<<endl;
                    iterate_cut(parents, i, i, color, myflag, vector_color, leaves_pixels);
                    color++;
                }
            }
        }

        ////cout<<"color : "<<color<< endl << "fin" <<endl;
        //primero pregunto si este nodo hoja no tiene padre si no tiene padre lo que sigue es darle color background
        int background = vector_color.size();

        for(int i = 0; i < leaves_pixels; i++)
        {
            if(parents[i]->parent==NULL){vector_color[i] = background;}
        }

       //cout<<"color : "<<color<<endl;
        int color_ = color;
        
        for (int i = 0; i < leaves_pixels; i++) //generar regiones de las hojas sueltas 
        {
            if(vector_color[i] == -1)
            {
                if(!option) //todas las hojas sueltas pertenecen a una misma región
                {
                    vector_color[i] = color_;
                }else{
                    vector_color[i] = color;    //cada hoja suelta genera una región diferente
                    color++;
                }   
            }
        }
        //cout<<"color : "<<color<<endl;


}



void visititeratead( const vector< Element* > &parents, vector<int> &newparents, vector<int> &vector_color, vector<int> &myflag,  vector<int> &level_down, int color, int node)
{

    myflag[node] = 1;
    vector_color[node] = color;
    if(isleaf(parents, level_down, node))
    {
        newparents.push_back(node);
    }else{
        for(int i = 0; i < parents[node]->children.size(); i++)
        {   
            ////cout <<"node "<<parents[node]->children[i]->index << " is child from "<< node << endl;
            visititeratead(parents, newparents, vector_color, myflag, level_down, color, parents[node]->children[i]->index);
        }
    }
    
}

int cutleveltree(const vector< Element* > &parents, vector< vector<int> > &newregions,  vector<int> &level_down,  vector<int> &vector_color, int max, int cut_tree, int leaves_pixels)
{
    vector<int> myflag(parents.size(),0);
    int color =0;
    
    for(int j = cut_tree; j > 0; j--)
    {
        for(int i = parents.size() -1; i >leaves_pixels; i--)
        {
            if(level_down[i] == j && !myflag[i] && (parents[i]->children.size()>0) && level_down[i]!=-1 )
            {
                ////cout<<"j : "<<j<<", i : "<<i<<", level_down["<<i<<"] = "<<level_down[i]<<"\tmyflag[i] "<<myflag[i]<<endl;
                vector<int> newparents;
                visititeratead(parents, newparents, vector_color, myflag, level_down, color, i);
                sort(newparents.begin(), newparents.end());
                newregions.push_back(newparents);
                color++;
            }
        }
    }

    vector<int> aux;
    for(int i = 0; i < leaves_pixels; i++)
    {
        if(!myflag[i])
        {
            vector_color[i] = color;
            aux.push_back(i);
        }
    }
    if(aux.size()>0) newregions.push_back(aux);
    //cout << "color: "<< color << endl;

    return color;        
}


int ColorRegionsHGB(const vector< Element* > &parents, vector<int> &vector_color, vector<int> &level_down, int max, int cut_tree, int leaves_pixels, bool option, int width, int height)
{
    //cout<<"#######  coloring regions #####"<<endl;
    //Empiezo a colorear
    //usar flag si el nodo ha sido visitado o no 
    //primero coloreo todas las hojas que tengan como hierarchy <= parametro
    //luego coloreo las hojas que no han sido coloreadas agregando cada hoja aislada a un color diferente

    //std::cout << "size-enter-jojojo: " << parents.size() << std::endl;

    //int param = 1;
    //int 
    //max = *max_element(level.begin(), level.end());
    //cout<<"Iniciando coloreado de regiones "<<endl;
    int color =0;
    vector<int> myflag(parents.size(),0);

    //print(parents);
    
    if(cut_tree==0)// corte a nivel de hojas
    {
        for(int i = 0; i < leaves_pixels; i++)
        {
            vector_color[i] = color++;
        }
    }
    else{
        //tomar regiones con nivel igual o menor del parametro establecido
        //debo empezar desde los nodos padres creados .. no considerar las hojas 
        //partir del numero de hojas 
        //cout << "cut_tree: "<< cut_tree<< endl;

        vector<vector<int>> newregions;
        color = cutleveltree(parents, newregions, level_down, vector_color, max, cut_tree, leaves_pixels);

        //cout << "color : "<< color << endl;

        //cout << "vector_color.size() "<< vector_color.size() << endl;
        vector <int > regionorder;    


        vector< pair <int,int> > regionlist;
        // index, middle_item


        int max_sizeregion = 0;
        int max_region = -1;
        for(int i = 0; i < newregions.size(); i++)
        {
            vector<int> aux = newregions[i];
            int actual_size = aux.size();
            if(actual_size > max_sizeregion) 
            {
                max_sizeregion = actual_size;
                max_region = i;
            }
            
        }

        int color =  vector_color[newregions[max_region][0]];
        int central_pixel = width/2*(height-1);
        int cx = central_pixel%width;
        int cy = central_pixel/width;
        //cout << "central_pixel:  "<< central_pixel << endl;
        //cout << "cx: " << cx <<", cy: "<< cy << endl;
        int max_distance = sqrt(pow(width/6, 2) + pow(height/6, 2));
        //cout << "max_distance: "<< max_distance << endl;

        vector < vector<int> > afterclean;
        int min_distance = width*height;
        int region  = -1;
        for(int i = 0; i < newregions.size(); i++)
        {
            vector<int> aux = newregions[i];

            if(aux.size() != max_sizeregion)
            {
                int actual_size = aux.size();
                int middle_item = aux[actual_size/2];
                ////cout << "size : "<< aux.size() << ", the middle_item : " << middle_item << endl;
                int actual_cx = middle_item%width;
                int actual_cy = middle_item/width;
                int actual_distance = sqrt(pow(abs(cx - actual_cx), 2) + pow(abs(cy - actual_cy), 2));

                if(actual_distance < min_distance)
                {   
                    min_distance = actual_distance;
                    region = i;
                }
            }
        }

        //si necesito para mostrar el proceso mantengo el vector color
        //por ahora lo voy a borrar

        vector_color.clear();
        vector_color.resize(leaves_pixels,0);

        vector<int> aux = newregions[region];
        
        //cout << "min pixel " << aux[0] << endl;
        
        int temp_i = int(aux[0]/width);
        int temp_j = int(aux[0]%width);

        //cout << "i: " << temp_i << ", j: " << temp_j << endl;

        //cout << "max pixel " << aux[aux.size()-1] << endl;
        temp_i = int(aux[aux.size()-1]/width);
        temp_j = int(aux[aux.size()-1]%width);

        cout << "i: " << temp_i << ", j: " << temp_j << endl;

        //save if it make sense ... dividing the image in 8 

        cout << "width "<< width << ", height " << height << endl;
        bool discard = false;
        int max_4 = int(width*(height/4));
        cout << "max_4 " << max_4 << endl;
        if(aux[aux.size()-1] <  max_4){discard = true;}   //if max in 1 2 3 or 4 i don't count on it 
        int min_9 = 3*max_4 + 1;
        cout << "min_9 " << min_9 << endl;
        if(aux[0] > min_9){discard = true;}            //if min in 8 9 10 or 11 

        if(discard == false)
        {
            //If it wasn't discarded then proceed to save this information
            cout << "saving this one ######################### " << endl;
            //open file and save this cut and the points 
            string text;
            ofstream file;
            file.open("test.txt");
            cout << "Type some text" << endl;
            getline(cin, text);
            file << text;
        }

        for(int j = 0; j < aux.size(); j++)
        {
            //cout << aux[j] << endl;
            vector_color[aux[j]] = 1;
        }

        //ubicar color de la region background para luego pasar las regiones que no son el corazon a el color del background
        /*En este caso funcionan las lineas porque las regiones se pueden apreciar manualmente pero pueden ordenar la region 2 y 3 para ser borradas en el filtrado */


    }



    return color;
}


int ColorRegionsHGB2(string filename, int op, int mycut, int typeComputeD, const vector< Element* > &parents, vector<int> &vector_color, vector<int> &level_down, int max, int cut_tree, int leaves_pixels, bool option, int width, int height)
{
    //cout<<"#######  coloring regions #####"<<endl;
    //Empiezo a colorear
    //usar flag si el nodo ha sido visitado o no 
    //primero coloreo todas las hojas que tengan como hierarchy <= parametro
    //luego coloreo las hojas que no han sido coloreadas agregando cada hoja aislada a un color diferente

    //std::cout << "size-enter-jojojo: " << parents.size() << std::endl;

    //int param = 1;
    //int 
    //max = *max_element(level.begin(), level.end());
    //cout<<"Iniciando coloreado de regiones "<<endl;
    int color =0;
    vector<int> myflag(parents.size(),0);

    //print(parents);
    
    if(cut_tree==0)// corte a nivel de hojas
    {
        for(int i = 0; i < leaves_pixels; i++)
        {
            vector_color[i] = color++;
        }
    }
    else{
        //tomar regiones con nivel igual o menor del parametro establecido
        //debo empezar desde los nodos padres creados .. no considerar las hojas 
        //partir del numero de hojas 
        //cout << "cut_tree: "<< cut_tree<< endl;

        vector<vector<int>> newregions;
        color = cutleveltree(parents, newregions, level_down, vector_color, max, cut_tree, leaves_pixels);

        //cout << "color : "<< color << endl;

        //cout << "vector_color.size() "<< vector_color.size() << endl;
        vector <int > regionorder;    


        vector< pair <int,int> > regionlist;
        // index, middle_item


        int max_sizeregion = 0;
        int max_region = -1;
        for(int i = 0; i < newregions.size(); i++)
        {
            vector<int> aux = newregions[i];
            int actual_size = aux.size();
            if(actual_size > max_sizeregion) 
            {
                max_sizeregion = actual_size;
                max_region = i;
            }
            
        }

        int color =  vector_color[newregions[max_region][0]];
        int central_pixel = width/2*(height-1);
        int cx = central_pixel%width;
        int cy = central_pixel/width;
        //cout << "central_pixel:  "<< central_pixel << endl;
        //cout << "cx: " << cx <<", cy: "<< cy << endl;
        int max_distance = sqrt(pow(width/6, 2) + pow(height/6, 2));
        //cout << "max_distance: "<< max_distance << endl;

        vector < vector<int> > afterclean;
        int min_distance = width*height;
        int region  = -1;
        for(int i = 0; i < newregions.size(); i++)
        {
            vector<int> aux = newregions[i];

            if(aux.size() != max_sizeregion)
            {
                int actual_size = aux.size();
                int middle_item = aux[actual_size/2];
                ////cout << "size : "<< aux.size() << ", the middle_item : " << middle_item << endl;
                int actual_cx = middle_item%width;
                int actual_cy = middle_item/width;
                int actual_distance = sqrt(pow(abs(cx - actual_cx), 2) + pow(abs(cy - actual_cy), 2));

                if(actual_distance < min_distance)
                {   
                    min_distance = actual_distance;
                    region = i;
                }
            }
        }

        //si necesito para mostrar el proceso mantengo el vector color
        //por ahora lo voy a borrar

        vector_color.clear();
        vector_color.resize(leaves_pixels,0);

        vector<int> aux = newregions[region];
        

        //save if it make sense ... dividing the image in 8 

        //cout << "width "<< width << ", height " << height << endl;
        bool discard = false;
        int max_4 = int(width*(height/4));
        int min_9 = 3*max_4 + 1;

        int minimum_i = width*height*100;
        int minimum_j = width*height*100;
        int maximum_i = -width*height*100;
        int maximum_j = -width*height*100;

        for(int j = 0; j < aux.size(); j++)
        {
            //cout << aux[j] << endl;
            vector_color[aux[j]] = 1;

            int current_i = int(aux[j]/width);
            int current_j = int(aux[j]%width);

            //cout << "i: " << current_i << ", j: " << current_j << endl;

            if(current_i < minimum_i){minimum_i = current_i;}
            if(current_j < minimum_j){minimum_j = current_j;}
            if(current_i > maximum_i){maximum_i = current_i;}
            if(current_j > maximum_j){maximum_j = current_j;}
        }

        //cout << "max_4 " << max_4 << ", min_9 " << min_9 <<  endl;
        if((maximum_i*height + maximum_j) < max_4){discard = true;}   //if max in 1 2 3 or 4 i don't count on it 
        if((minimum_i*height + minimum_j) > min_9){discard = true;}            //if min in 8 9 10 or 11 
        //cout << "min_i " << minimum_i << ", min_j " << minimum_j << ", max_i " << maximum_i << ", max_j " << maximum_j <<  endl;
        if(discard == false)
        {
            //If it wasn't discarded then proceed to save this information
            //cout << "saving this one ######################### "<< filename << endl;
            //cout << "min pixel " << aux[0] << endl;
            
            string saving_file = "/data/leticia/training_da/";
            //cout << "hgb2: filename: " <<filename << endl;

            //filename = filename.substr(28);
            filename = filename.substr(33);
            string filename2 = filename.erase(filename.size() - 4);
            //cout << "hgb2: filename2: " <<filename2 << endl;

            //open file and save this cut and the points 
            //string text;
            fstream file;
            file.open(saving_file + filename2 + ".txt", std::ios_base::app);
            //cout << "Type some text" << endl;
            //getline(cin, text);
            //op cut_level D xmin,ymin xmax,ymax
            string text = to_string(width) + "\t" +to_string(height) +"\t"+ to_string(mycut) + "\t"  + to_string(op)  + "\t" + to_string(typeComputeD) + "\t" + to_string(minimum_i) + "\t" + to_string(minimum_j) + "\t" + to_string(maximum_i) + "\t" + to_string(maximum_j) + "\n";
            //file.write(text.data(), text.size());
            
            file << text;
            file.close();
        }

        
        //ubicar color de la region background para luego pasar las regiones que no son el corazon a el color del background
        /*En este caso funcionan las lineas porque las regiones se pueden apreciar manualmente pero pueden ordenar la region 2 y 3 para ser borradas en el filtrado */


    }



    return color;
}

void updateTree(vector< Element* > &parents, vector <int> &euler, vector <int> &level, vector <int> &R, vector <int> &leaves, vector <int> &level_down, int leaves_pixels, int root )
{
    //Para el saliency map uso LCA euler, level, R
    //update tree information
    //cout<<"Inicia actualización de información del arbol"<<endl;
    euler.clear();
    level.clear();
    R.clear();
    R.resize(parents.size(),-1);
    
    iterate_pre(parents, euler, level, root);

    for(int i =0; i < euler.size(); i++)
    {
        int index = euler[i];
        if(R[index]==-1){R[index] = i;}
    }
    //Para corte y coloreado uso level_down y leaves
    leaves.clear();
    level_down.clear();
    leaves.resize(parents.size(),0);
    level_down.resize(parents.size(),0);
    count_leaves(parents, leaves, level_down , leaves_pixels, root);

}

Graph saliencymapQCT(Graph g, vector<int> &vector_color, vector< Element* > &parents, int leaves_pixels, int level_tree ) // leaves_pixels : cantidad de hojas 
{

    int parameter_number = 500; //la regla de cantidad de elementos que voy a considerar .. si tengo menor o igual a parameter_number los elimino 
    int root= parents[parents.size()-1]->index;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler, level;
    iterate_pre( parents, euler, level, root);
    //cout <<"despues de iterate_pre "<< endl;
    double max = *max_element(level.begin(), level.end());
    //cout<<"max : "<<max <<endl;
    vector<int> R(parents.size(), -1);
    //cout <<"size euler : "<<euler.size() <<", size level : "<<level.size() <<endl;

    for(int i =0; i < euler.size(); i++)
    {
        int index = euler[i];
        if(R[index]==-1){R[index] = i;}
    }
    //for(int i = 0; i < R.size();i++){//cout<<"R["<<i<<"]: " <<R[i]<<", ";}
    //Procesamiento:

    //debo procesar y luego crear un grafo con los mismos edges pero diferentes pesos
    // saliency map basado en algoritmo 1 del paper
    
    //cout<<"##############################   Saliency map ##############################"<<endl;
    Graph saliency(g.V,g.E); 
    //saliency.V = g.V;
    //saliency.E = saliency.V;
    vector< pair<double, iPair> >::iterator it;    
    for (it=g.edges.begin(); it!=g.edges.end(); it++)
    //for (it=MST.begin(); it!=MST.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        //int w = it->first;  //S[{x, y}] := level[LCA(T, {x}, {y})] -1; //level means altitude
        int z = LCA(u,v,R, euler, level);
        //Saliency.push_back({g.altitude[z], {u, v}});  //de la version inicial del saliency map basado en regiones
        saliency.addEdge(u, v, g.altitude[z]);
        ////cout<< "u : "<<u <<", v : " << v << ", alitud :  "<<g.altitude[z]<<endl;
    }
    //pasar del grafo saliency al binary tree usando kruskal para reconstruir parents
    
    //guardar el saliency 
    //int mst_wt = saliency.kruskalMST();
    int QCT_size = 0;
    //parents.clear(); //borrar parents
    //canonizeQBT(saliency.parents_graph, parents, saliency.altitude, leaves_pixels, QCT_size);    

    //procesamiento despues del saliency map
    vector<int> leaves(parents.size(),0); //uso vector auxiliar de hojas 
    vector<int> level_down(parents.size(),0);
    //cout<<"parents size : "<<parents.size()<<endl;

    count_leaves(parents, leaves, level_down , leaves_pixels, root);
    //print_vector(leaves);

    //cout<<"####   children filter ####"<<endl;
    childrenFilter(parents, leaves, leaves_pixels, parameter_number); //children filtery
    updateTree(parents, euler, level, R, leaves, level_down, leaves_pixels, root); //update information
  
    max = *max_element(level.begin(), level.end());
    //cout<<"max : "<<max<<endl;
    vector_color.resize(leaves_pixels,-1);
    //cout<<"vector_color size "<<vector_color.size()<<endl;
    
    //para el grafo creado el nivel maximo es 4  el root
    //el coloreado es segun el corte que se hace en el arbol el mayor nivel es el nivel del root que genera una sola region 
    //si el corte es debajo de la raiz 
    double max_level = *max_element(level_down.begin(), level_down.end());
    //cout<<"level down max : "<<max_level<<endl;
    //cout<<"leaves_pixels "<<leaves_pixels<<endl;
    //int level_tree = 7; //representa el corte en el arbol para el pintado //9116 para 1_graph
    bool separate_regions = true; //true: genera regiones separadas por cada hoja suelta //false: el conj de hojas sueltas genera una sola region
    //ColorRegions(parents, R, vector_color, level_down, max, cut_tree, leaves_pixels,separate_regions); 
    //cambiando level_tree = max_level -1;
    level_tree = max_level -1;
    //int width = 500;    int height = 500;
    ColorRegions(parents, R, vector_color, level_down, max, level_tree, leaves_pixels, false); 
    //cout<<endl<<"end saliency map"<<endl;    

    return saliency;
}

Mat ProcessGraph(int type, Mat image, string mygraphfile, string finalgraphfile, int level_tree)
{
    //cout << endl <<"Iniciando ProcessGraph ------------------------------------------------------------------- "<<endl;
    int width = 0, height = 0;
    int V = 0, E = 0;

    Graph g(V,E);
    Graph saliency(V,E);
    vector<int> vector_color;

    switch(type) {
        case 0  :
            //cout<<"Reading graph " << endl; //leer desde imagen 
            g = readGraph(V, E, width, height, mygraphfile);
            break; //optional
        case 1 :
            //cout<<"Saving graph " << endl; //leer desde .graph
            g = saveGraph(image, V, E, height, width);
            break; //optional
        
        default : //Optional
            cout<<"This type doesn't exist."<<endl;
    }

    // processing a graph
    //cout<<"height : "<< height<< ", width : " << width << "V : " << V << ", E : " << E << endl;
    //cout << "Edges of MST are \n";
    Graph pr_g = g;
    int mst_wt  = pr_g.kruskalMST();
    //cout <<endl<<"Weight of MST is " << mst_wt<<"\n";
    //cout <<"count : "<< pr_g.cnt<<"\n";
    //cout<<"size of MST "<<pr_g.MST.size()<<endl; 
    int vertex_number = (width)*(height);
    //cout<<"vertex_number "<<vertex_number<<endl;

    int QCT_size = 0;
    vector< Element* >  qct_parents;

    //cout<<"qbt"<<endl;
    //print(pr_g.parents_graph);
    //cout << endl << "altitud"<<endl;
    //print_vector(pr_g.altitude);
    //cout<<endl;
    canonizeQBT(pr_g.parents_graph, qct_parents, pr_g.altitude, vertex_number, QCT_size);  
    //cout<<" canonizeQBT "<< endl;
    saliency = saliencymapQCT(pr_g, vector_color, qct_parents, vertex_number, level_tree); //Saliency
    //cout << " saliency "<<endl;
    int sizecolors = 100;
    Mat final_image = show_coloring(vector_color, sizecolors, width, height);
    saveGraphFile(finalgraphfile, saliency, height, width);
    return final_image;

}

// Depth-first preprocessing
int32_t LCApreprocessDepthFirst(JCctree *CT, int32_t node, int32_t depth, int32_t *nbr, int32_t *rep, int32_t *Euler, int32_t *Represent, int32_t *Depth, int32_t *Number)
{
  int32_t son;
  JCsoncell *sc;
  if (CT->tabnodes[node].nbsons > -1) {
    (*nbr)++;
    Euler[*nbr] = node;
    Number[node] = *nbr;
    Depth[node] = depth;
    Represent[*nbr] = node;
    (*rep)++;
    for (sc = CT->tabnodes[node].sonlist; sc != NULL; sc = sc->next)    {
      son = sc->son;
      LCApreprocessDepthFirst(CT, son, depth+1, nbr, rep, Euler, Represent, Depth, Number);
      Euler[++(*nbr)] = node;
    }
  }
  free(sc);


  return *nbr;
}



int32_t ** LCApreprocess(JCctree *CT, int32_t *Euler, int32_t *Depth, int32_t *Represent, int32_t *Number, int32_t *nbR, int32_t *lognR)
{
  //O(n.log(n)) preprocessing
  int32_t nbr, rep, nbNodes;
  nbr = -1; // Initialization number of euler nodes
  rep = 0;
  nbr = LCApreprocessDepthFirst(CT, CT->root, 0, &nbr, &rep, Euler, Represent, Depth, Number);
  nbNodes = rep;

  // Check that the number of nodes in the tree was correct
  assert((nbr+1) == (2*nbNodes-1));

  int32_t nbRepresent = 2*nbNodes-1;
  int32_t logn = (int32_t)(ceil(log((double)(nbRepresent))/log(2.0)));
  *nbR = nbRepresent;
  *lognR = logn; 

  int32_t i,j,k1,k2;
  int32_t *minim = (int32_t *)calloc(logn*nbRepresent, sizeof(int32_t));
  int32_t **Minim = (int32_t **)calloc(logn, sizeof(int32_t*));
  Minim[0] = minim;

  for (i=0; i<nbRepresent-1; i++) {
    if (Depth[Euler[i]] < Depth[Euler[i+1]]) {
      Minim[0][i] = i;
    } else {
      Minim[0][i] = i+1;
    }
  }
  Minim[0][nbRepresent-1] = nbRepresent-1;

  for (j=1; j<logn; j++) {
    k1 = 1<<(j-1);
    k2 = k1<<1;
    Minim[j] = &minim[j*nbRepresent];
    for (i=0; i<nbRepresent; i++) {
      if ((i+ k2) >= nbRepresent) {
  Minim[j][i] = nbRepresent-1;
      } else {
          if (Depth[Euler[Minim[j-1][i]]] <= Depth[Euler[Minim[j-1][i+k1]]]) {
            Minim[j][i] = Minim[j-1][i];
          } else {
            Minim[j][i] = Minim[j-1][i+k1];
          }
      }
    }
  }
  return Minim;
}


int32_t LowComAncFast(int32_t n1, int32_t n2, int32_t *Euler, int32_t *Number, int32_t *Depth, int32_t **Minim)
{
  int32_t ii, jj, kk, k;

  ii = Number[n1];
  jj = Number[n2];
  if (ii == jj)
    return ii;

  if (ii > jj) {
    kk = jj;
    jj = ii;
    ii = kk;
  }

  k = (int32_t)(log((double)(jj - ii))/log(2.));

  if (Depth[Euler[Minim[k][ii]]] < Depth[Euler[Minim[k][jj-(1<<(k))]]]) {
    return Number[Euler[Minim[k][ii]]];
  } else {
    return Number[Euler[Minim[k][jj-(1<<k)]]];
  }
}

void autoremovenodes(vector< Element* > &parents, vector<int> &level_down, Element* q )
{
   
    while(parents[q->index]->children.size() == 0){
        
        level_down[q->index] = -1;
        int remove_node = q->index; 
        parents[remove_node]->children.clear();

        int father = parents[remove_node]->parent->index;
        vector<int> aux;
        for(int j = 0; j < parents[father]->children.size(); j++)
        {
            int actual_value = parents[father]->children[j]->index;
            if(actual_value!=remove_node) 
            {
                aux.push_back(actual_value);
            }
        }
        parents[father]->children.clear();
        for(int j = 0; j < aux.size(); j++){parents[father]->addChild(parents[aux[j]]);}

        q = parents[remove_node]->parent;
        parents[remove_node]->parent = NULL;
        
    }
}

void sendtoRoot(vector< Element* > &parents, int node, int root)
{
    parents[root]->addChild(parents[node]);
    parents[node]->parent = parents[root];
}

void iterate_remove(vector< Element* > &parents, vector<int> &leaves, vector<int> &level_down, int node, int root){

    ////cout << "node: "<< node << endl;
    if(isleaf(parents, level_down, node))
    {
        ////cout<<"node "<< node << " is leaf"<<endl;
        sendtoRoot(parents, node, root);
    }else{
        ////cout <<"node "<<node << " is not leaf"<< endl;
        // guardar al padre y hacer un pop_front al padre
        level_down[node] = -1;
        leaves[node] = -1;
        //parents[node]->parent = NULL;
        
        for(int i = 0; i < parents[node]->children.size(); i++)
        {   
            ////cout <<"node "<<parents[node]->children[i]->index << " is child from "<< node << endl;
            iterate_remove(parents, leaves, level_down, parents[node]->children[i]->index, root);
        }
        //padre de node sin hijo de 
        //parents[node]->removeAllChildren();
        
        
    }
    
}


void removePixels(vector< Element* > &parents, vector<int> &level_down, vector<int> &leaves, int node, int root)
{
    int father = parents[node]->parent->index; 
    if(father == root) return;
    sendtoRoot(parents, node, root);
    removefromfather(parents, node, father);
    updatenodeinformation(parents, leaves, level_down, father, root);

    /*
    Element* p = parents[father];
    while (p->index != root)
    {
        Element* actual = p;
        leaves[actual->index] -= 1;
        p = p->parent;
        if(actual->children.size()==0) // si no tiene hijos  
        {
            actual->parent = NULL;
            level_down[actual->index] = -1;
            int max_level = -2;
            for(int i = 0; i < p->children.size(); i++)
            {
                if(level_down[p->children[i]->index] > max_level ) {max_level = level_down[p->children[i]->index];}
            }
            level_down[p->index] = max_level + 1; 
        }
    }
    */

}


void regionFilter(vector< Element* > &parents, vector<int> &leaves, int leaves_pixels, vector<int> &level_down, int width, int height)
{

    //cout<< "#################  parents.size(): "<<parents.size() << endl;
    //int width = 256;  //int height = 216;
    int root = parents.size()-1; 
    //cout << "root : "<< root << endl;

    //cout << "parents[root]->children.size(): " <<  parents[root]->children.size() << endl;

    // Quitar regiones que esten fuera un area especifica de la region ####################################################################
    
    //int axmin = width/(13.7);       int aymin = height/(11.4);      int areamin = axmin*aymin/5;
    //int axmax = width/(13.7/2.5);   int aymax = (height/(11.4/2));  int areamax = axmax*aymax*12;
    
    for(int i = 0; i < parents[root]->children.size(); i++)
    {
        int actual_node = parents[root]->children[i]->index;
        //cout << "nodo hijo ["<< i << "] : " << actual_node << endl; 
    }

    //cout << "---------------------------------------"<<endl;
    
    /* 
    int areamin = 200;   int areamax = 55000; //19600 : 1ra vesion // 20900 : 2da version;  //pixeles : 55296
    //cout << "areamin: " << areamin << ", areamax: " << areamax << endl;

    
    for(int i=leaves.size()-3; i >leaves_pixels; i--)
    {
        //quitar nodo retirar subarbol y pasar nodos hoja a la raiz // funcion para saber que nodos hoja tiene un nodo en especifico nodo hijo 
        if( (leaves [i] > areamax  ) && level_down[i] !=-1  && leaves[i] !=-1 ) 
        { 
            int node = i;
            int father = parents[node]->parent->index;
            //cout << "leaves["<< i << "] = "<< leaves[i]<< ", level_down[" << i <<"] = " << level_down[i] <<endl;

            iterate_remove(parents, leaves, level_down, node, root);
            removefromfather(parents, node, father);
            parents[node]->removeAllChildren();
            parents[node]->parent = NULL;
            updatenodeinformation(parents, leaves, level_down, father, root);
        }
    }
    */


    //cout << "---------------------------------------"<<endl;    
    ////cout<< "despues de quitar regiones  mayores a "<< areamax << endl;
    

    // Solo se debe considerar una region 

    // tomo los 10 primeros nodos y elimino los nodos que son hijos

    //cout << "cantidad de hijos del root : "<< parents[root]->children.size() << endl;

    //if(parents[root]->children[i]->index > leaves_pixels)
    //{
        ////cout << "-- "<< parents[root]->children[i]->index << endl;
        //
        ////tomo a su hijo y pregunto a que cuadrante pertenece y evaluo
    //}
    

    // Quitar regiones fuera del cuadrado  ####################################################################
    
    int x = width/5;    int y = height/5;  
    for(int node = 0, lastnode = leaves_pixels-1;  node < y*width; node++, lastnode--)  //factor arriba y factor abajo
    {
        removePixels(parents, level_down, leaves, node, root);
        removePixels(parents, level_down, leaves, lastnode, root);
    }
    for(int nodey = y; nodey <= 4*y; nodey++)   //bordes laterales
    {
        for(int nodex = 4*x; nodex < width; nodex++)
        {
            int node     = nodey*width + nodex;
            int secondnode = node - 4*x;
            removePixels(parents, level_down, leaves, node, root);
            removePixels(parents, level_down, leaves, secondnode, root);
        }
    }
    

}




Mat sm_code(string mygraphfile, int level_tree)
{
    //cout << endl << "Iniciando sm code ------------------------------------------------------------------- "<<endl;
    int32_t rs, cs, rsIm, csIm, ds, N, mode,i,j,x,y;
    pcell p;
    double *Fv;
    float *F, tmp;
    
    graphe *g;
    g = ReadGraphe(mygraphfile, &Fv);

    rs = g->rs; cs = g->cs; N = rs * cs;
    if( (rs == -1) || (cs == -1) ) {
        //fprintf(stderr, "%s: %s is not the graph of an image\n", argv[0], argv[1]);
        exit(0);
    }

    //cout<<endl<<"cs : "<<cs <<", rs : "<<rs <<endl;
    //convertir grafo a saliencyTree
    int32_t logn, nbRepresent, u,/* x, y,*/ c1;
    int32_t nbarcs = g->ind;
    double *STAltitudes;
    JCctree *ST;
    /* Data structure for fast computation of LCA (least common ancestor) */
    int32_t *Euler, *Depth, *Represent, *Number, **Minim;
    STAltitudes = (double *)calloc(2*g->nbsom, sizeof(double));
    //if (STAltitudes == NULL){
    //  fprintf(stderr,"ultramopen cannot allocate memory for STAltitudes\n");
    //  exit(0);
    //}
    
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;

    //cout<<" :  taille : "<<taille<<", narcs : "<<narcs<<endl;
    saliencyTree(g, &ST, STAltitudes);  //aqui leti aquii //STAltitude vector que guarda los pesos de los edges

    //imprimiendo arbol
    ////cout<<" :  taille : "<<taille<<", narcs : "<<narcs << " g->ind : "<< g->ind << ", g->nbmaxar : "<< g->nbmaxar<<endl;  ////cout << " ST->nbnodes : " << ST->nbnodes<< ",  nbsoncells : " << ST->nbsoncells << "ST->root : "<< ST->root<<endl;

    /*for(int k = 0; k < ST->nbnodes; k++)
    {
        if (ST->tabnodes[k].nbsons == -1 ) //cout << "ST->tabnodes["<< k << "].nbsons is a deleted node"<<endl;
        if (ST->tabnodes[k].nbsons < 2 &&  ST->tabnodes[k].nbsons >0) //cout << k<< ST->tabnodes[k].nbsons<<endl;
    }*/

    vector< Element* >  parents_sm;  //guardar como hijos propios usando mi definicion de vector parents
    vector<int> vector_color; //vector color para usar en el coloreado
    for(int k = 0; k < taille; k++) //guardar nodos hojas 
    {
        Element* parent = new Element(k);
        parents_sm.push_back(parent);
    }
    
    for(int k = taille; k < ST->nbnodes; k++) //crear nodos padres y sus hijos a la vez
    {
        JCsoncell *tmp_sm;
        tmp_sm = ST->tabnodes[k].sonlist;
        
        Element* parent = new Element(k);
        parent->addChild(parents_sm[tmp_sm[0].son]);
        parent->addChild(parents_sm[tmp_sm[1].son]);
        parents_sm.push_back(parent);
        
    }
    
    //print(parents_sm);
    int parameter_number = 2; //la regla de cantidad de elementos que voy a considerar .. si tengo menor o igual a parameter_number los elimino 
    int root= ST->root;
    int leaves_pixels = taille;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler, level;
    vector<int> R(parents_sm.size(), -1);
    //Procesamiento:
    vector<int> leaves(parents_sm.size(),0); //uso vector auxiliar de hojas 
    vector<int> level_down(parents_sm.size(),0);
    //cout<<"parents size : "<<parents_sm.size()<<endl;

    count_leaves(parents_sm, leaves, level_down , leaves_pixels, root);
    //print_vector(leaves);
    ////cout<<"##############################   children filter ##############################"<<endl;
    //childrenFilter(parents, leaves, leaves_pixels, parameter_number); //children filtery
    updateTree(parents_sm, euler, level, R, leaves, level_down, leaves_pixels, root); //update information
    //cout <<"size euler : "<<euler.size() <<", size level : "<<level.size() <<endl;
    double max = *max_element(level.begin(), level.end());
    //cout<<"max : "<<max<<endl;

    
    vector_color.resize(leaves_pixels,-1);
    //cout<<"vector_color size "<<vector_color.size()<<endl;
    
    //para el grafo creado el nivel maximo es 4  el root
    //el coloreado es segun el corte que se hace en el arbol el mayor nivel es el nivel del root que genera una sola region 
    //si el corte es debajo de la raiz 
    double max_level = *max_element(level_down.begin(), level_down.end());
    //cout<<"level down max : "<<max_level<<endl;
    //cout<<"leaves_pixels "<<leaves_pixels<<endl;

    //haciendo cambio de level_tree = max_level -1;
    level_tree = max_level -30;
    ColorRegions(parents_sm, R, vector_color, level_down, max, level_tree, leaves_pixels, false); // el ultimo parametro representa el corte en el arbol para el pintado

    /*Euler = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
    Represent = (int32_t *)calloc(2*ST->nbnodes-1, sizeof(int32_t));
    Depth = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
    Number = (int32_t *)calloc(ST->nbnodes, sizeof(int32_t));
  
    if ((Euler == NULL) || (Represent == NULL)  || (Depth == NULL) || (Number == NULL)) {
    //fprintf(stderr, "%s : malloc failed\n", F_NAME);
        return;
    }
    //cout<< endl << "---------------------------------------------------------" <<endl;
    Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);
    //cout<< " Depth[ST->root] :" << Depth[ST->root] << ", Number[ST->root] : " << Number[ST->root] << endl;
    //cout<< " Depth[taille-1] :" << Depth[taille-1] << ", Number[taille-1] : " << Number[taille-1] << endl;
    //cout<< " Depth[0] :" << Depth[0] << ", Number[0] : " << Number[0] << endl;
    //cout<< " STAltitudes[110590] : "<< STAltitudes[110590]<<endl;

    // inicia algoritmo de saliency map  #############################################################  GENERAR EXTERNO
    // For any edge of g
    for(u = 0; u < nbarcs; u++){
        x = g->tete[u]; y = g->queue[u];
        c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
        //cout<<"c1 "<<c1<<endl;
        setweight(g,u,STAltitudes[c1]);
    }
    */
    int sizecolors = 100;
    Mat final_image = show_coloring(vector_color, sizecolors, rs, cs);

    return final_image;
}
void usage(char *arg){
    printf("#################################################################\n\n");
    printf("USAGE: %s IN.graph  OUT.graph TYPE \n",arg);
    printf("Types \n");
    printf("1  Min  \n");
    printf("2  Max  \n");
    printf("3  Mean  \n");
    printf("4  Median \n");
    printf("#################################################################\n\n");
}

inline int mini(int a, int b){
	return (a<=b)?a:b;
}

inline int maxi(int a, int b){
	return (a>=b)?a:b;
}

inline float minf(float a, float b){
	return (a<=b)?a:b;
}

inline float maxf(float a, float b){
	return (a>=b)?a:b;
}

inline double mind(double a, double b){
	return (a<=b)?a:b;
}

inline double maxd(double a, double b){
	return (a>=b)?a:b;
}

inline void * myCalloc(size_t num, size_t size)
{
	void * ref = calloc(num, size);
	if (ref==NULL)
	{ fprintf(stderr,"Cannot allocate enough memory\n"); exit(1); }
	return ref;
}

int * computeNodeArea(JCctree *CT)
{
	int * area = (int *)myCalloc(CT->nbnodes, sizeof(int));
	int i;
	JCctreenode *tabnodes = CT->tabnodes;
	for(i=0;i<CT->nbnodes;++i)
	{
		if(i<=CT->nbnodes/2)
		{
			area[i]=1;
		}else{
			JCsoncell * soncell;
			for( soncell = tabnodes[i].sonlist; soncell!=NULL; soncell=soncell->next)
			{
				int child = soncell->son;
				area[i]+=area[child];	
			}
            
		}
	}
    
	return area;
}



void saliency(JCctree *ST, graphe *g, double * STAltitudes)
{  
  int32_t logn, nbRepresent, u, x, y, c1;
  int32_t nbarcs = g->ind;
  int32_t *Euler, *Depth, *Represent, *Number, **Minim; /* Data structure for fast computation of LCA (least common ancestor) */
  Euler = (int32_t *)myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
  Represent = (int32_t *)myCalloc(2*ST->nbnodes-1, sizeof(int32_t));
  Depth = (int32_t *)myCalloc(ST->nbnodes, sizeof(int32_t));
  Number = (int32_t *)myCalloc(ST->nbnodes, sizeof(int32_t));
  Minim = LCApreprocess(ST, Euler, Depth, Represent, Number, &nbRepresent, &logn);

  for(u = 0; u < nbarcs; u++){ // For any edge of g
    x = g->tete[u]; y = g->queue[u];
    c1 = Represent[LowComAncFast(x, y, Euler, Number, Depth, Minim)];
    setweight(g,u,STAltitudes[c1]);
  }  
  free(Euler);
  free(Represent);
  free(Depth);
  free(Number);
  free(Minim[0]);
  free(Minim);


}


JCctree * UpperBoundBenjamin(graphe *g, double ratio)
{

    int nbPix = g->nbsom;   //char * stop;
    int minArea = (int)(ratio* nbPix);  //int minArea = (int)strtol(area, &stop, 10);
    //cout <<"Area : "<< minArea << endl;
	double *level = (double *)myCalloc(2*g->nbsom, sizeof(double));
	JCctree *CT;
    int * MST;
    
  	saliencyTree2(g,&CT,level, &MST);
	int * nodeArea = computeNodeArea(CT);
    int nbLeaves = CT->nbnodes/2+1;
    JCctreenode *tabnodes = CT->tabnodes;
   
    for(int i=nbLeaves ; i<CT->nbnodes; ++i)
    {
        JCsoncell * soncell = tabnodes[i].sonlist;
		int child1 = soncell->son;
		int child2 = soncell->next->son;
        if(nodeArea[i] < minArea || nodeArea[child1] < minArea || nodeArea[child2] < minArea)
            level[i]=0;
        
    }
    
    graphe * g2 = MSTToGraph(g, MST, level);
    free(MST);
    componentTreeFree(CT);
    saliencyTree2(g2,&CT,level, &MST);
    saliency(CT, g, level);
    //SaveGraphe(g, argv[argc-1],Av);
    free(MST);
	free(nodeArea);
	//terminegraphe(g);
    terminegraphe(g2);
	//free(Av);
	free(level); 

    return CT;
}


template<char delimiter>
class WordDelimitedBy : public std::string
{};

Mat HGB_Edward_old(string myfile, string imagefile, int op, string finalgraphfile, int level_tree, int typeComputeD)
{

    
    //cout << endl << "Iniciando HGB_Edward code, op = "<< op <<" -------------------------------------------------------------- "<<endl;
    //cout << endl << myfile << endl;
    graphe *g;
    double *Av;
    JCctree *CT;

    if(!(op>=1&&op<=21)) {
        //cout<<"Wrong type!, valid types"<<endl<<"1 Min "<<endl<< "2 Max " << endl << " 3 "<<endl<<" 4 "<< " 5 "<< " 6 "<< " 7 "<< " 8 "<<" ... "<<endl;
        exit(1);
    }
    
    g = ReadGraphe(myfile, &Av);

    //read and save my graph
    int width = 0, height = 0;
    int V = 0, E = 0;

    Graph g2(V,E);
    g2 = readGraph(V, E, width, height, myfile);
    int mst_wt  = g2.kruskalMST();

    //cout << "iniciando image to graph ------------------------------------------------"<<endl;

    //ImagetoGraph(imagefile, &Av);

    int th = 10; // 10 -> 500
    double param = 0.01; // 0.01 -> 0.05
    //cout << "typeComputeD : " << typeComputeD<< endl;
    intervalSegmentation(g, op, th, param, typeComputeD  ); // the new HGB with interval //new HGB

    //incrementalSegmentationInterval(g, op);  //linea 1 llamando al HGB //Old HGB
    double ratio = 0.005;
    CT = UpperBoundBenjamin(g, ratio); //linea 2   llamando a benjamin

    //cout << "nbsons:  "<< CT->tabnodes[CT->root].nbsons << "**********************************************************************" << endl;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    int rs = g->rs; int cs = g->cs; int N = rs * cs;
    //cout<<"rs : "<<rs<<", cs: "<<cs<<", N: "<<N<<endl;

    //cout << "guardando archivo del grafo antes de aplicar coloreado y filtrado de regiones : " << finalgraphfile << endl;
    //SaveGraphe(g, finalgraphfile, Av);
    //mysavegraphe(g, finalgraphfile, Av);

    //cout<<"aqui despues de mysavegraphe"<<endl;

    //aqui empiezo a guardar el grafo para poder modificarlo  ###############################

    vector< Element* >  parents_hgb;  //guardar como hijos propios usando mi definicion de vector parents 
    vector<int> vector_color; //vector color para usar en el coloreado
    //cout<< "taille: "<<taille<<endl;

    for(int k = 0; k < taille; k++) //guardar nodos hojas 
    {
        Element* parent = new Element(k);
        parents_hgb.push_back(parent);
    }
    
    //cout<<"despues del primer for"<< endl << "CT->nbnodes: "<<CT->nbnodes<< endl <<"CT->root: "<<CT->root<<endl;

    for(int k = taille; k < CT->nbnodes; k++) //crear nodos padres y sus hijos a la vez
    {
        JCsoncell *tmp_sm;
        tmp_sm = CT->tabnodes[k].sonlist;
        Element* parent = new Element(k);
        parent->addChild(parents_hgb[tmp_sm[0].son]);
        parent->addChild(parents_hgb[tmp_sm[1].son]);
        parents_hgb.push_back(parent);
    }

    int parameter_number = 2; //la regla de cantidad de elementos que voy a considerar .. si tengo menor o igual a parameter_number los elimino 
    int root= CT->root;
    int leaves_pixels = taille;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler, level, R(parents_hgb.size(), -1);
    terminegraphe(g);
    componentTreeFree(CT);
    
    //Despues que ya termine de guardar el grafo 
    //aplicar saliency map 
    //aqui empieza todo el ultimo cambio que me pidio edward master thesis thesis!!

    //cout<<"##############################   Saliency map ##############################"<<endl;
    
    /*
    //cout<< "V: "<< g2.V << ", E: " << g2.E<< endl;

    //Graph saliency(g2.V,g2.E);
    vector< pair<double, iPair> >::iterator it;

    //hacer preprocesamiento de LCA

    vector< pair<int, iPair> > saliency_;

    int root_saliency= parents_hgb[parents_hgb.size()-1]->index;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler_saliency, level_saliency;
    iterate_pre( parents_hgb, euler_saliency, level_saliency, root);
    //cout <<"despues de iterate_pre "<< endl;
    double max_saliency = *max_element(level_saliency.begin(), level_saliency.end());
    //cout<<"max_saliency : "<<max_saliency <<endl;
    vector<int> R_saliency(parents_hgb.size(), -1);
    //cout <<"size euler : "<<euler_saliency.size() <<", size level : "<<level_saliency.size() <<endl;

    for(int i =0; i < euler_saliency.size(); i++)
    {
        int index = euler_saliency[i];
        if(R_saliency[index]==-1){R_saliency[index] = i;}
    }

    for (it=g2.edges.begin(); it!=g2.edges.end(); it++)
    {
        int u = it->second.first;
        int v = it->second.second;
        ////cout<< "u : "<<u <<", v : " << v << endl;
        int z = LCA(u,v,R_saliency, euler_saliency, level_saliency);
        ////cout << ", altitud :  "<<g2.altitude[z]<<endl;
        //saliency.addEdge(u, v, g2.altitude[z]);
        
        iPair p = make_pair (u,v);
        saliency_.push_back(make_pair(g2.altitude[z], p));
        
    }

    Mat saliency_image = generate_image(saliency_, width, height);
    */



    //std::string delimiter = "";
    //std::string token = myfile.substr(0, myfile.find(delimiter)); // token is "scott"
    
    //string _patient = token;
    ////cout<< "patient: " << _patient << endl;
    //imwrite("saliency_HGB_"+ to_string(op) + "_typeD_" + to_string(typeComputeD) + ".png", saliency_image);

    //cout<< "despues del saliency" << endl;

    vector<int> leaves(parents_hgb.size(),0), level_down(parents_hgb.size(),0); //uso vector auxiliar de hojas 
    
    //cout<<"parents size : "<<parents_hgb.size()<<endl;
    count_leaves(parents_hgb, leaves, level_down , leaves_pixels, root);  //print_vector(leaves);

    double max = *max_element(level_down.begin(), level_down.end());
    ////cout<<"level down max : "<< max_level <<  ", leaves_pixels : "  <<  leaves_pixels  << endl;
    level_tree = max - level_tree;
    //cout << "max: "<< max << ", level_tree: "<< level_tree <<  endl;
    //double max = *max_element(level.begin(), level.end());  //cout<<"max : "<<max<<endl;
    vector_color.resize(parents_hgb.size(),-1);  //cout<<"vector_color size "<<vector_color.size()<<endl;

    //cortar en el nivel level_tree 

    int colors = ColorRegionsHGB(parents_hgb, vector_color, level_down, max, level_tree, leaves_pixels, false, rs, cs); // el ultimo parametro representa el corte en el arbol para el pintado
    
    //cout << "rs: "<< rs << ", cs: "<< cs << endl;
    
    
    //Mat final_image = show_coloring(vector_color, colors, rs, cs);
    Mat final_image = show_grayscale_coloring(vector_color, colors, rs, cs);
    return final_image;
}


Mat HGB_Edward(string myfile, string imagefile, int op, string finalgraphfile, int level_tree, int typeComputeD)
{

    Mat final_image = Mat(3, 2, CV_8UC3);
    return final_image;
}

Mat generateboundbox2(string filename, int op, int level_tree, int typeComputeD )
{

    
    //Read file with information 
    graphe *g;
    double *Av;
    JCctree *CT;
    cout << filename << "\t" << op << "\t" << level_tree << "\t" << typeComputeD << endl;
    int mycut = level_tree;
    //read and save my graph
    g = TexttoGraphe(filename, &Av); //Create graph and return graph
    string line;

    //apply HGB 

    
    int th = 10; // 10 -> 500
    double param = 0.01; // 0.01 -> 0.05
    //cout << "typeComputeD : " << typeComputeD<< endl;
    intervalSegmentation(g, op, th, param, typeComputeD  ); // the new HGB with interval //new HGB
    
    //incrementalSegmentationInterval(g, op);  //linea 1 llamando al HGB //Old HGB
    double ratio = 0.005;
    CT = UpperBoundBenjamin(g, ratio); //linea 2   llamando a benjamin
    int root= CT->root;
    //cout << "nbsons:  "<< CT->tabnodes[CT->root].nbsons << "**********************************************************************" << endl;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    int rs = g->rs; int cs = g->cs; int N = rs * cs;
    terminegraphe(g);
    vector< Element* >  parents_hgb;  //guardar como hijos propios usando mi definicion de vector parents 
    vector<int> vector_color; //vector color para usar en el coloreado
    //cout<< "taille: "<<taille<<endl;

    for(int k = 0; k < taille; k++) //guardar nodos hojas 
    {
        Element* parent = new Element(k);
        parents_hgb.push_back(parent);
    }
    
    //cout<<"despues del primer for"<< endl << "CT->nbnodes: "<<CT->nbnodes<< endl <<"CT->root: "<<CT->root<<endl;
    int endvalue =  CT->nbnodes;
    for(int k = taille; k < endvalue; k++) //crear nodos padres y sus hijos a la vez
    {
        JCsoncell *tmp_sm;
        tmp_sm = CT->tabnodes[k].sonlist;
        Element* parent = new Element(k);
        parent->addChild(parents_hgb[tmp_sm[0].son]);
        parent->addChild(parents_hgb[tmp_sm[1].son]);
        parents_hgb.push_back(parent);
    }
    componentTreeFree(CT);

    int parameter_number = 2; //la regla de cantidad de elementos que voy a considerar .. si tengo menor o igual a parameter_number los elimino 
    
    int leaves_pixels = taille;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler, level, R(parents_hgb.size(), -1);
    
    vector<int> leaves(parents_hgb.size(),0), level_down(parents_hgb.size(),0); //uso vector auxiliar de hojas 
    count_leaves(parents_hgb, leaves, level_down , leaves_pixels, root);  //print_vector(leaves);
    double max = *max_element(level_down.begin(), level_down.end());
    level_tree = max - level_tree; //bajo level tree desde el root 
    vector_color.resize(parents_hgb.size(),-1);  //cout<<"vector_color size "<<vector_color.size()<<endl;
    
    //cortar en el nivel level_tree 
    
    int colors = ColorRegionsHGB2(filename, op, mycut, typeComputeD, parents_hgb, vector_color, level_down, max, level_tree, leaves_pixels, false, rs, cs); // el ultimo parametro representa el corte en el arbol para el pintado
    
    Mat final_image = show_grayscale_coloring(vector_color, colors, rs, cs);
    //imwrite( '/home/leticia/Documents/Windows_UCSP/results/hgb/image.jpg', final_image );
    //Mat final_image

    
    return final_image;
    /*  ************************************************************************************
    */  //************************************************************************************
    //found the regions (xmin, ymin) (xmax, ymax) that get the max region possible or the whole image is returned
}



void generateboundbox3(string filename, int op, int level_tree, int typeComputeD )
{

    
    //Read file with information 
    graphe *g;
    double *Av;
    JCctree *CT;
    cout << filename << "\t" << op << "\t" << level_tree << "\t" << typeComputeD << endl;
    int mycut = level_tree;
    //read and save my graph
    g = TexttoGraphe(filename, &Av); //Create graph and return graph
    string line;

    //apply HGB 

    
    int th = 10; // 10 -> 500
    double param = 0.01; // 0.01 -> 0.05
    //cout << "typeComputeD : " << typeComputeD<< endl;
    intervalSegmentation(g, op, th, param, typeComputeD  ); // the new HGB with interval //new HGB
    
    //incrementalSegmentationInterval(g, op);  //linea 1 llamando al HGB //Old HGB
    double ratio = 0.005;
    CT = UpperBoundBenjamin(g, ratio); //linea 2   llamando a benjamin
    int root= CT->root;
    //cout << "nbsons:  "<< CT->tabnodes[CT->root].nbsons << "**********************************************************************" << endl;
    int32_t taille = g->nbsom;
    int32_t narcs = g->ind;
    int rs = g->rs; int cs = g->cs; int N = rs * cs;
    terminegraphe(g);
    vector< Element* >  parents_hgb;  //guardar como hijos propios usando mi definicion de vector parents 
    vector<int> vector_color; //vector color para usar en el coloreado
    //cout<< "taille: "<<taille<<endl;

    for(int k = 0; k < taille; k++) //guardar nodos hojas 
    {
        Element* parent = new Element(k);
        parents_hgb.push_back(parent);
    }
    
    //cout<<"despues del primer for"<< endl << "CT->nbnodes: "<<CT->nbnodes<< endl <<"CT->root: "<<CT->root<<endl;
    int endvalue =  CT->nbnodes;
    for(int k = taille; k < endvalue; k++) //crear nodos padres y sus hijos a la vez
    {
        JCsoncell *tmp_sm;
        tmp_sm = CT->tabnodes[k].sonlist;
        Element* parent = new Element(k);
        parent->addChild(parents_hgb[tmp_sm[0].son]);
        parent->addChild(parents_hgb[tmp_sm[1].son]);
        parents_hgb.push_back(parent);
        
    }
    componentTreeFree(CT);

    int parameter_number = 2; //la regla de cantidad de elementos que voy a considerar .. si tengo menor o igual a parameter_number los elimino 
    
    int leaves_pixels = taille;
    //cout<<endl<<"root : "<<root<<endl;
    vector <int> euler, level, R(parents_hgb.size(), -1);
    
    vector<int> leaves(parents_hgb.size(),0), level_down(parents_hgb.size(),0); //uso vector auxiliar de hojas 
    count_leaves(parents_hgb, leaves, level_down , leaves_pixels, root);  //print_vector(leaves);
    double max = *max_element(level_down.begin(), level_down.end());
    level_tree = max - level_tree; //bajo level tree desde el root 
    vector_color.resize(parents_hgb.size(),-1);  //cout<<"vector_color size "<<vector_color.size()<<endl;
    
    //cortar en el nivel level_tree 
    
    int colors = ColorRegionsHGB2(filename, op, mycut, typeComputeD, parents_hgb, vector_color, level_down, max, level_tree, leaves_pixels, false, rs, cs); // el ultimo parametro representa el corte en el arbol para el pintado

}