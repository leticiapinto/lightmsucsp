dir="/home/leticia/datasets/emidec-dataset-1.0.1"
graph_dir="/data/leticia/EMIDEC/graph_files"

#dir="../ACDCtraining"

for f in "$graph_dir"/*; do
        y=${f%.*}
        remainder="${f##*/}"
        #echo $f
        #echo $remainder
        var2=${remainder%????}
        #echo $y
        #echo $var2
        ./main    "$y".txt "$var2"

done