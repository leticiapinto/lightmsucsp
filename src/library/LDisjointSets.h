// To represent Disjoint Sets
struct DisjointSets
{

    int n;
    int QBT_size;

    // Constructor.

    DisjointSets(int n)
    {
        // Allocate memory
        this->n = n;
    }
    void qbt_makeSet(int q)
    {
        Element* parent1 = new Element(q);
        qbt_parents.push_back(parent1);
        altitude.push_back(-1);
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
    
    int qbt_union(int cx, int cy, int weight)
    {
        Element* parent1 = new Element(QBT_size);
        parent1->addChild(qbt_parents[cx]);
        parent1->addChild(qbt_parents[cy]);
        qbt_parents.push_back(parent1);
        if(weight!=last_altitude)
        {
            cnt_altitude++;
            last_altitude= weight;
        }
        altitude.push_back(cnt_altitude);

        QBT_size++; 

        return QBT_size -1;
    }
};