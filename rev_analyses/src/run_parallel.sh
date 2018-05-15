#!/bin/bash

rm -rf temp_output
mkdir temp_output

for rep in {1..20}
do
    echo "rep = ${rep}; seed(${rep}); source(\"src/chromoploid.Rev\");" | rb > temp_output/${rep}.out &
done
