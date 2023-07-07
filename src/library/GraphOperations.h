int getEdge(int node)
{   
    return node; 
}

int weight_node(int node)
{
    return altitude[getEdge(node)];
}

void canonizeQBT(vector< Element* > &qbt_parents,  vector< Element* > &qct_parents, int vertex_number, int QCT_size)
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
        //cout<<"node i : " <<i;
        if(qbt_parents[i]->parent!=NULL && qbt_parents[i]->children.size()>0)
        {
            int p = qct_parents[i]->parent->index;
            //cout<<"; weight_node("<<p<<") : "<<weight_node(p)<<", weight_node("<<i<<") : "<<weight_node(i)<<endl;
            if(weight_node(p) == weight_node(i))
            {
                for(int j =0; j < qct_parents[i]->children.size(); j++)
                {
                    //el padre de estos hijos es p
                    int tmp_child = qct_parents[i]->children[j]->index;
                    //cout<<"tmp_child : "<<tmp_child<<endl;
                    //actualizar hijos de p
                    qct_parents[p]->addChild(qct_parents[tmp_child]);    //le agrego este hijo al nuevo padre
                    qct_parents[tmp_child]->parent = qct_parents[p]; //cambio el padre de estos hijos
                    //qct_parents[i]->removeChild(qct_parents[tmp_child]); //le quito los hijos a este nodo que se elimina
                }
                qct_parents[p]->removeChild(qct_parents[i]); //quito a este nodo de su padre 
                qct_parents[i]->removeAllChildren();
                qct_parents[i]->parent= NULL;
                qct_parents[i]->in_tree = 0;
                //cout<<" nodo i"<<i<<endl;  
            }
            
        }

    }

    //exit(0);
    //print(qct_parents);
    //cout << "each non-leaf and non-root node n " <<endl;
    // If needed, build the list of children
    //cout<< "qct_parents : "<<qct_parents.size()<<endl;
    for(int i = 0; i < qct_parents.size(); i++)
    {   
        //cout<<"i :"<< i<<endl;
        //cout<<"qct_parents[i]->parent : "<<qct_parents[i]->parent<<endl;
        if(qct_parents[i]->parent!= NULL)
        {
            int p = qct_parents[i]->parent->index; //cout<< qct_parents[i]->parent->index<<endl;
            qct_parents[p]->addChild(qct_parents[i]);
        }
    }

    //cout<<endl<<"altitude"<<endl;
    //print_vector(altitude);

    //cout<<"altitude cantidad : "<<cnt_altitude<<", altitude last_altitude : "<<last_altitude<<endl;
    //cout<<endl<<"-----------------------------size "<<qct_parents.size();
}
