#!/bin/bash
# the collection of MATLAB code;


if [ $# -eq 9 ]
then
	model=$1
	doc=$2

	## the candidate structures
	pattern_str=$3		# a string
	pattern_tot=(${pattern_str})
	len=${#pattern_tot[@]}

	## the different running mode
	## including compute, reflag, data_check, diagram_check, boundary
	calculate=$4
	if [ ${calculate} == 'diagram_check' ] || [ ${calculate} == 'boundary' ]
	then
		sflag='diagram_'
	else
		sflag=''
	fi

	## the times of refining space grid; the original step divides 2
	## only for the mode: boundary
	refine_num=$5

	## the split of computational region in tau direction
	tau_num_tot=($6)
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
	gamma_num_tot=($7)
	gamma_len=${#gamma_num_tot[@]}
	if [ ${gamma_len} -lt ${len} ]
	then
		ref=${gamma_num_tot[$((${gamma_len}-1))]}
		for ((i=${gamma_len};i<${len};i++))
		do
			gamma_num_tot[$i]=${ref}
		done
	fi
	
	## calculate: yes or not (parameters, parameters_tot, hamilton for tidy)
	cal_flag=$8


	if [ $9 == 'data' ]
	then
		for ((i=0;i<${len};i++))
		do
			pattern=${pattern_tot[$i]}
			tau_num=${tau_num_tot[$i]}
			gamma_num=${gamma_num_tot[$i]}
			fname0=${doc}/${pattern}/${pattern}_${sflag}hamilton.txt
			fname1=${doc}/${pattern}/${pattern}_${sflag}hamilton0.txt
			fname2=${doc}/${pattern}/${pattern}_${sflag}hamilton1.txt
			if [ ${refine_num} -le 1 ]
			then
				for ((j1=0;j1<${tau_num};j1++))
				do
					for ((j2=0;j2<${gamma_num};j2++))
					do
						log=log/${pattern}_${j1}_${j2}.log
						err=err/${pattern}_${j1}_${j2}.err
						nohup matlab -nosplash -nodesktop -nodisplay -r \
						"pfc('${model}', '${pattern}', '${refine_num}', \
						'${j1}', '${tau_num}', '${j2}', '${gamma_num}', \
						'${calculate}', '${cal_flag}')" 1>${log} 2>${err} &
					done
				done
				if [ -e ${fname0} ]
				then
					if [ ${refine_num} -eq 0 ]
					then
						cp ${fname0} ${fname1}	# 0
					else
						cp ${fname0} ${fname2}	# 1
					fi
				fi
			else
				diffham=finish/diff.txt
				diff ${fname1} ${fname2} > ${diffham}
				if [ ! -z "${diffham}" ]	## check whether these files are different
				then
					for ((j1=0;j1<${tau_num};j1++))
					do
						for ((j2=0;j2<${gamma_num};j2++))
						do
							log=log/${pattern}_${j1}_${j2}.log
							err=err/${pattern}_${j1}_${j2}.err
							nohup matlab -nosplash -nodesktop -nodisplay -r \
							"pfc('${model}', '${pattern}', '${refine_num}', \
							'${j1}', '${tau_num}', '${j2}', '${gamma_num}', \
							'${calculate}', '${cal_flag}')" 1>${log} 2>${err} &
						done
					done
					cp ${fname2} ${fname1}
					cp ${fname0} ${fname2}
				else
					echo "computing ${calculate}(${refine_num}) is unnecessary"
				fi
			fi
		done
		## check whether these progresses are finished
		finish_tot=0
		for ((i=0;i<${len};i++))
		do
			let finish_tot=${finish_tot}+${tau_num_tot[$i]}*${gamma_num_tot[$i]}
		done
		while true
		do
			finish=0
			for ((i=0;i<${len};i++))
			do
				pattern=${pattern_tot[$i]}
				tau_num=${tau_num_tot[$i]}
				gamma_num=${gamma_num_tot[$i]}
				for ((j1=0;j1<${tau_num};j1++))
				do
					for ((j2=0;j2<${gamma_num};j2++))
					do
						finish_file=finish/${pattern}-[${j1}-${tau_num}]-[${j2}-${gamma_num}].txt
						if [ -e ${finish_file} ]
						then
							let finish+=1
						fi
					done
				done
			done
			if [ ${finish} -eq ${finish_tot} ]
			then
				for ((i=0;i<${len};i++))
				do
					pattern=${pattern_tot[$i]}
					rm finish/${pattern}*.txt
				done
				break
			else
				sleep 10s
			fi
		done
		echo "---> ${calculate} ${finish}/${finish_tot}"
	elif [ $9 == 'diagram' ]
	then
		log=log/diagram.log
		err=err/diagram.err
		# here refine_num just for a name saving figure
		nohup matlab -nosplash -nodesktop -nodisplay -r \
		"diagram('${model}', '${pattern_str}', '${refine_num}')" \
		1>${log} 2>${err} &
		## check whether this progress is finished
		while true
		do
			finish_file=finish/diagram.txt
			if [ -e ${finish_file} ]
			then
				rm ${finish_file}
				break
			else
				sleep 10s
			fi
		done
		echo "---> The code (diagram) finish"
	else
		echo "the ninth parameter: data or diagram"
	fi
else
	echo "require nine parameters"
fi
