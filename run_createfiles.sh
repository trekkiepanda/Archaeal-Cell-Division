#!/bin/bash
#call on bashrc to use aliases
shopt -s expand_aliases
#source ~/.bashrc

NAME=$1

if [ "$#" -ne 1 ]; then
	echo "!!name of directroy must be given as input variable!!"
	exit
fi

GAP=1
DISTIMESTEP=6000
ANGLE=27.5


for i in `seq 5 1 8`; do		# sigma
for j in `seq 10 5 40`; do		# packing
for k in `seq 1 1 10`; do		# seed num
	mkdir ./SIGMA=${i}_PACK=${j}_SEED=${k}/
	cp ./in.local ./SIGMA=${i}_PACK=${j}_SEED=${k}/
	cp ./helix.py ./SIGMA=${i}_PACK=${j}_SEED=${k}/
	cp ./spherecoords_48002.xyz ./SIGMA=${i}_PACK=${j}_SEED=${k}/
	cd ./SIGMA=${i}_PACK=${j}_SEED=${k}/
	python3 helix.py ${ANGLE} ${GAP} ${DISTIMESTEP} ${i} ${j} ${k}
	#rm submit_myriad.pbs
	echo "#!/bin/bash -l" >> submit_myriad.pbs
	echo "#$ -S /bin/bash" >> submit_myriad.pbs
	echo "#$ -l h_rt=48:00:00" >> submit_myriad.pbs
	echo "#$ -l s_rt=46:00:00" >> submit_myriad.pbs
	echo "#$ -l mem=2G" >> submit_myriad.pbs
	echo "#$ -l tmpfs=15G" >> submit_myriad.pbs
	echo "#$ -N "${NAME}"_SIGMA="${i}"_PACK="${j}"_SEED="${k} >> submit_myriad.pbs
	echo "#$ -pe mpi 4" >> submit_myriad.pbs
	echo "#$ -wd /home/zcaptly/Scratch/output/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/" >> submit_myriad.pbs
	echo "" >> submit_myriad.pbs
	echo "function Cleanup ()" >> submit_myriad.pbs
	echo "{" >> submit_myriad.pbs
	echo "  trap \"\" SIGUSR1 EXIT # Disable trap now we're in it" >> submit_myriad.pbs
	echo "  # Clean up task" >> submit_myriad.pbs
	echo "  rsync -rltv \$TMPDIR/ \$HOME/Scratch/output/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/" >> submit_myriad.pbs
	echo "  exit 0" >> submit_myriad.pbs
	echo "}" >> submit_myriad.pbs
	echo "trap Cleanup SIGUSR1 EXIT # Enable trap" >> submit_myriad.pbs
	echo "cd \$TMPDIR" >> submit_myriad.pbs
	echo "rsync -rltv \$HOME/input/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/ \$TMPDIR/" >> submit_myriad.pbs
	echo " " >> submit_myriad.pbs
	echo "#module load lammps/16Mar18/basic/intel-2018" >> submit_myriad.pbs
	echo " " >> submit_myriad.pbs
	echo "gerun \$HOME/Scratch/np-toolkit/lammps/src/lmp_mpi -in in.local" >> submit_myriad.pbs
	#rm submit_grace.pbs
	echo "#!/bin/bash -l" >> submit_grace.pbs
	echo "#$ -S /bin/bash" >> submit_grace.pbs
	echo "#$ -l h_rt=48:00:00" >> submit_grace.pbs
#	echo "#$ -l s_rt=46:00:00" >> submit_grace.pbs
	echo "#$ -l mem=2G" >> submit_grace.pbs
	echo "#$ -l tmpfs=15G" >> submit_grace.pbs
	echo "#$ -N "${NAME}"_SIGMA="${i}"_PACK="${j}"_SEED="${k} >> submit_grace.pbs
	echo "#$ -pe mpi 48" >> submit_grace.pbs
	echo "#$ -wd /home/zcaptly/Scratch/output/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/" >> submit_grace.pbs
	echo "" >> submit_grace.pbs
#	echo "function Cleanup ()" >> submit_grace.pbs
#	echo "{" >> submit_grace.pbs
#	echo "  trap \"\" SIGUSR1 EXIT # Disable trap now we're in it" >> submit_grace.pbs
#	echo "  # Clean up task" >> submit_grace.pbs
#	echo "  rsync -rltv \$TMPDIR/ \$HOME/Scratch/output/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/" >> submit_grace.pbs
#	echo "  exit 0" >> submit_grace.pbs
#	echo "}" >> submit_grace.pbs
#	echo "trap Cleanup SIGUSR1 EXIT # Enable trap" >> submit_grace.pbs
#	echo "cd \$TMPDIR" >> submit_grace.pbs
#	echo "rsync -rltv \$HOME/input/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/ \$TMPDIR/" >> submit_grace.pbs
#	echo " " >> submit_grace.pbs
#	echo "#module load lammps/16Mar18/basic/intel-2018" >> submit_grace.pbs
	echo "cd \$HOME/Scratch/output/"${NAME}"/SIGMA="${i}"_PACK="${j}"_SEED="${k}"/" >> submit_grace.pbs
	echo " " >> submit_grace.pbs
	echo "gerun \$HOME/Scratch/np-toolkit/lammps/src/lmp_mpi -in in.local" >> submit_grace.pbs
  cd ..
done
done
done
