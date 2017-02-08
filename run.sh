make
./strassens > data

awk -f graph.awk data

#xgraph recurGraph.tr strassenGraph.tr -geometry 600x800
