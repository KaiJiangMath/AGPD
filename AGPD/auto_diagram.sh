#!/bin/bash
######################################################################
##########                                                  ##########
##########        first version: June.25 2021               ##########
##########        Author: Kai Jiang, Wei Si                 ##########
##########        Email: kaijiang@xtu.edu.cn				##########
##########                                                  ##########
######################################################################
#
# plot phase diagram automatically (based on nohup.sh);
clear
## the running history
if [ ! -d "log/" ];
then
	mkdir log
fi
## the running error
if [ ! -d "err/" ];
then
	mkdir err
fi
## check whether the progress is finished
if [ ! -d "finish/" ];
then
	mkdir finish
fi

## the computational model
# lb, lp, ok, leibler
model=lb
doc=${model}_results

## the total number for data_check; 0: zero time
data_check_tot=(2 1)
## the total number for boundary; 0: zero time
boundary_tot=2

## the split of computational region in tau direction
## complete according to the last element
tau_num_tot="1 1 3 2 3 3"
#tau_num_tot="1"
## the split of computational region in gamma direction
## complete according to the last element
gamma_num_tot="1 1 3 2 3 3"
#gamma_num_tot="1"

## all candidate structures
#pattern_str="lam hex gyroid bcc fcc A15 sigma fdddCube hcp"
#pattern_str="lam hex gyroid bcc fcc A15 sigma"
pattern_str="lam hex gyroid bcc fcc A15"
#pattern_str="lam hex gyroid"
#pattern_str="lam hex 12fold LQ6 LQS6 C3 12i6o 8i10o sq squ sqv sqw Ls"
#pattern_str="lam hex"

## obtain array
pattern_tot=(${pattern_str})
len=${#pattern_tot[@]}


## refine_num: the times of refining space grid; the original step divides 2
if [ $# -eq 0 ]
then
	## delete all files about the boundary coordinates
	pattern_tot=(${pattern_str})
	len=${#pattern_tot[@]}
	for ((i=0;i<${len};i++))
	do
		pattern=${pattern_tot[$i]}
		for ((i1=0;i1<${boundary_tot};i1++))
		do
			file=${doc}/${pattern}/boundary${i1}.txt
			if [ -f ${file} ]; then rm ${file}; fi
		done
	done


	echo "start"
	echo `date`
	## compute a coarse phase diagram
	## nohup can check whether the progresses are finished
	calculate=compute
	./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} 0 "${tau_num_tot}" "${gamma_num_tot}" yes data
	## plot phase diagram
	./nohup.sh ${model} ${doc} "${pattern_str}" tidy 0 1 1 hamilton data
	./nohup.sh ${model} ${doc} "${pattern_str}" reflag 0 1 1 no data
	./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} 0 1 1 no diagram

	## using the suitable convergence value as the initial value
	calculate=data_check
	echo "${model}: data_check (${data_check_tot[0]})"
	for ((j1=0;j1<${data_check_tot[0]};j1++))
	do
		echo "---> ${j1}"
		./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} 0 "${tau_num_tot}" "${gamma_num_tot}" yes data
		## plot phase diagram
		./nohup.sh ${model} ${doc} "${pattern_str}" tidy 0 1 1 hamilton data
		./nohup.sh ${model} ${doc} "${pattern_str}" reflag 0 1 1 no data
		num=$((${j1}+1))
		./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} ${num} 1 1 no diagram
	done

	## refine the phase boundary
	calculate=boundary
	echo "${model}: boundary (${boundary_tot})"
	for ((j3=0;j3<${boundary_tot};j3++))
	do
		echo "boundary ${j3}"
		## estimate all waiting computing parameters
		./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} ${j3} 1 1 no data
		## clear up parameters
		./nohup.sh ${model} ${doc} "${pattern_str}" tidy ${j3} 1 1 parameters data
		## compute a refining phase diagram
		## nohup can check whether the progresses are finished
		./nohup.sh ${model} ${doc} "${pattern_str}" 'compute' ${j3} "${tau_num_tot}" "${gamma_num_tot}" yes data
		## plot phase diagram
		./nohup.sh ${model} ${doc} "${pattern_str}" tidy ${j3} 1 1 hamilton data
		num=$((${data_check_tot[0]}+${j3}+1+${j3}*${data_check_tot[1]}))
		./nohup.sh ${model} ${doc} "${pattern_str}" ${calculate} ${num} 1 1 no diagram
	
		## check data again
		echo "---> data_check (${data_check_tot[1]})"
		for ((j4=0;j4<${data_check_tot[1]};j4++))
		do
			echo "---> ${j4}"
			./nohup.sh ${model} ${doc} "${pattern_str}" 'data_check' ${j3} "${tau_num_tot}" "${gamma_num_tot}" yes data
			## plot phase diagram
			./nohup.sh ${model} ${doc} "${pattern_str}" tidy ${j3} 1 1 hamilton data
			./nohup.sh ${model} ${doc} "${pattern_str}" reflag ${j3} 1 1 no data
			num=$((${data_check_tot[0]}+${j3}+1+${j3}*${data_check_tot[1]}+${j4}+1))
			./nohup.sh ${model} ${doc} "${pattern_str}" 'data_check' ${num} 1 1 no diagram
		done
	done
	echo `date`
elif [ $# -eq 1 ]
then
	## the split of computational region in tau direction
	tau_num_tot=(${tau_num_tot})
	tau_len=${#tau_num_tot[@]}
	if [ ${tau_len} -lt ${len} ]
	then
		ref=${tau_num_tot[$((${tau_len}-1))]}
		for ((i=${tau_len};i<${len};i++))
		do
			tau_num_tot[$i]=${ref}
		done
	fi
	## the split of computational region in gamma direction
	gamma_num_tot=(${gamma_num_tot})
	gamma_len=${#gamma_num_tot[@]}
	if [ ${gamma_len} -lt ${len} ]
	then
		ref=${gamma_num_tot[$((${gamma_len}-1))]}
		for ((i=${gamma_len};i<${len};i++))
		do
			gamma_num_tot[$i]=${ref}
		done
	fi
	
	if [ $1 == 'clean' ]
	then
		for ((i=0;i<${len};i++))
		do
			pattern=${pattern_tot[$i]}
			rm -r ${model}_results/${pattern}/check/
			rm -r ${model}_results/${pattern}/ini_check/
		done
	elif [ $1 == 'log' ]
	then
		for ((i=0;i<${len};i++))
		do
			pattern=${pattern_tot[$i]}
			tau_num=${tau_num_tot[$i]}
			gamma_num=${gamma_num_tot[$i]}
			for ((j1=0;j1<${tau_num};j1++))
			do
				for ((j2=0;j2<${gamma_num};j2++))
				do
					log=log/${pattern}_${j1}_${j2}.log
					tail -20 ${log} | less
				done
			done
		done
	elif [ $1 == 'err' ]
	then
		for ((i=0;i<${len};i++))
		do
			pattern=${pattern_tot[$i]}
			tau_num=${tau_num_tot[$i]}
			gamma_num=${gamma_num_tot[$i]}
			for ((j1=0;j1<${tau_num};j1++))
			do
				for ((j2=0;j2<${gamma_num};j2++))
				do
					err=err/${pattern}_${j1}_${j2}.err
					tail -20 ${err} | less
				done
			done
		done
	elif [ $1 == 'txt' ]
	then
		for ((i=0;i<${len};i++))
		do
			pattern=${pattern_tot[$i]}
			txt=${model}_results/${pattern}/${pattern}_hamilton.txt
			less ${txt}
		done
	elif [ $1 == 'cleanout' ]
	then
		rm -r log err finish *.log
#		rm -r lb_results
	elif [ $1 == 'help' ]
	then
		evince README.pdf
	else
		echo "the parameter: clean, log, err, txt, cleanout or help"
	fi
else
	echo "zero or one parameter"
fi
