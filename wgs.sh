#!/bin/bash

module load R
chmod 777 main_prog.R
chmod 777 post_edit.R
chmod 777 gmatrix
chmod 777 wombat


if [ $# -eq 0 ]
  then
    echo "No parameter file supplied"
    exit 1
fi

cp $1 parameter.R 
rm -r $1.sol
mkdir $1.sol

name=$(grep  "PHENO" parameter.R | awk '{split($0,a,"\""); print a[4]}')

if [ -z "$2" ] || [ $2 -eq 1 ]
then
	echo 'Analysing a single family'
	./main_prog.R 0
	./gmatrix par.dat > gmat.log
	mv animal.gin $name.gin
	./wombat --blup --dense2 --batch wf.par > wom.log
	mv RnSoln* gmat.log wom.log $1.sol
	mv SumModel.out FixSolutions.out $1.sol

	if [ -f WOMBAT.log ]; then rm WOMBAT.log; fi
	if [ -f recoded1.bin ]; then rm recoded1.bin; fi
	if [ -f recoded2.bin ]; then rm recoded2.bin; fi
	if [ -f recoded3.bin ]; then rm recoded3.bin; fi
	if [ -f SumEstimates.out ]; then rm SumEstimates.out; fi
	if [ -f SumPedigree.ou ]; then rm SumPedigree.ou; fi
	if [ -f Residuals.dat ]; then rm Residuals.dat; fi
	if [ -f $name.gin ]; then rm $name.gin; fi
	if [ -f pheno.dat ]; then rm pheno.dat; fi
	if [ -f map.dat ]; then rm map.dat; fi
	if [ -f pop_100.d ]; then rm pop_100.d; fi
	if [ -f rep.mrk ]; then rm rep.mrk; fi


else 
	echo 'Analysing multiple families'
 	x=$(grep  "PHENO" parameter.R | awk '{split($0,a,"\""); print a[2]}')
	awk '{ a[$2]++ } END { for (b in a) { print b } }' $x > famID

	for i in $(cat famID)
	do
		./main_prog.R $i
		./gmatrix par.dat > gmat.log
		mv animal.gin $name.gin
		./wombat --blup --dense2 --batch wf.par > wom.log

		mv RnSoln* gmat.log wom.log $1.sol
		mv SumModel.out FixSolutions.out $1.sol

		if [ -f WOMBAT.log ]; then rm WOMBAT.log; fi
		if [ -f recoded1.bin ]; then rm recoded1.bin; fi
		if [ -f recoded2.bin ]; then rm recoded2.bin; fi
		if [ -f recoded3.bin ]; then rm recoded3.bin; fi
		if [ -f SumEstimates.out ]; then rm SumEstimates.out; fi
		if [ -f SumPedigree.ou ]; then rm SumPedigree.ou; fi
		if [ -f Residuals.dat ]; then rm Residuals.dat; fi
		if [ -f $name.gin ]; then rm $name.gin; fi
		if [ -f pheno.dat ]; then rm pheno.dat; fi
		if [ -f map.dat ]; then rm map.dat; fi
		if [ -f pop_100.d ]; then rm pop_100.d; fi
		if [ -f rep.mrk ]; then rm rep.mrk; fi

		echo " ... Family $i completed ..."
	done
fi
