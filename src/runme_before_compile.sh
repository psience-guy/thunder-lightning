#!/bin/bash
MASTER=thunder-master

for i in include Makefile MACHINES a.GLOBAL b.FUNCTIONS c.SYSTEM d.FUNCTIONS_EXTRA e.FDATA g.XC_FUNCTIONALS h.SOLVESH j.ASSEMBLERS l.SCF o.OUTPUT u.UTIL
do
    if [ -e $i ]
    then
	rm $i
    fi
    ln -s ../../${MASTER}/src/$i
done
