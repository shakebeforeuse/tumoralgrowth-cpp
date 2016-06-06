#!/bin/bash

BRANCHES="unico multiple unoPorFrontera"
BRANCHES_DYNAMICDOMAIN="unoPorFronteraDominioCreciente"

SIZES="1000 2000 4000 8000 12000"
NGENS="100 1000 2000 4000"
TASKS=16
STEP=2

#For each branch
for branch in $BRANCHES; do

	git checkout $branch
	
	make clean
	make
	
	mkdir -p tmp/$branch
	
	for size in $SIZES; do
		echo "Running for size $size..."
		./speedup $size $TASKS $STEP 1000 > tmp/$branch/s$size.txt
	done
done

#Dynamic domain versions
for branch in $BRANCHES_DYNAMICDOMAIN; do

	git checkout $branch
	
	make clean
	make
	
	mkdir -p tmp/$branch
	
	for it in $NGENS; do
		echo "Running $it generations..."
		./speedup 4000 $TASKS $STEP $it > tmp/$branch/i$it.txt
	done
done

git checkout master
