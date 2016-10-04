#!/bin/bash

usage(){
      echo -e "\nUsage:  AnnotateR.sh -i <input_directory> -o <output_directory> [options]"
      echo -e "\n\tGiven a directory of fasta format protein or DNA files this script will determine genome"
      echo -e "\tcompleteneess and search sequences against both a CAZy HMM database (dbCAN) and a custom"
      echo -e "\tmetabolic HMM set (Karthik). Parsed output is ready for input to AnnoVisR."
      echo -e "\n\tOptions:"
      echo -e "\t-i: Directory of fasta format files to use as input (required)"
      echo -e "\t-o: Directory where output files will be written (required, will be created if absent)"
      echo -e "\t-p: Input directory contains DNA contigs; Run prodigal in single genome mode to predict proteins"
      echo -e "\t-c: Do not run CheckM for genome completeness (Runs by Default)"
      echo -e "\t-C: Do not run CAZy HMM set against proteins (Runs by Default)"
      echo -e "\t-M: Do not run Karthik's Metabolic HMM set against proteins (Runs by Default)"
      echo -e "\t-t: Number of threads to use (Default: 6)"
      echo -e "\t-e: E-value cutoff for HMM scoring (Default: 10)"
      echo -e "\t-h: Display this message and exit"
      echo -e "\nby: Spencer Diamond September 9, 2016"
}

## Default Arguments and Database Paths ##
threads=6
evalue=10
checkm_run=T
cazy_run=T
metabolic_run=T

cazy_hmm_database="/home/sdiamond/database/CAZy_HMM/dbCAN-fam-HMMs.txt"
metabolic_hmm_database="/home/sdiamond/database/Custom_HMM/Metabolic_HMMs/"
metabolic_hmm_cutoffs="/home/sdiamond/database/Custom_HMM/Metabolic_HMM_Cutoff_Scores.txt"
## Default Arguments and Database Paths ##

while getopts ":hi:o:pcCMt:e:" opt; do
  case $opt in
    h)
      usage
      exit 1
      ;;
    i)
      input_dir=${OPTARG}
      ;;
    o)
      output_dir=${OPTARG}
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
    p)
      prod_run="T"
      ;;
    c)
      checkm_run="F"
      ;;
    C)
      cazy_run="F"
      ;;
    M)
      metabolic_run="F"
      ;;
    t)
      threads=${OPTARG}
      ;;
    :) if [[ $OPTARG == "t" ]]; then
	:
       fi
      ;;
    e)
      evalue=${OPTARG}
      ;;
    :) if [[ $OPTARG == "e" ]]; then
	:
       fi
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [ "$#" -lt 1 ]
	then
	usage
	exit 1
fi

#### Begin Define Functions ####

prodigal_fn(){
	prodigal -i ${contigs} -a ${contigs}.faa -m -p single > /dev/null 2>&1
  echo -e "Genes predicted for ${contigs}..."
}


checkm_fn(){
  checkm lineage_wf -f ${output_dir}/checkm_output/checkm_summary.txt --genes -x faa -t 10 --pplacer_threads 10 ${input_dir} ${output_dir}/checkm_output
}


hmm_CAZy_fn(){ ###Can modify this function to accept modular HMM database inputs possibly use if statements with set up databases
	hmmscan --domtblout ${proteins}_out_domtbl.txt --cpu 6 --domE ${evalue} ${cazy_hmm_database} ${proteins} > /dev/null 2>&1
	wait
	sh /home/sdiamond/database/CAZy_HMM/hmmscan-parser.sh ${proteins}_out_domtbl.txt > ${proteins}_out_domtbl.parse
	wait
	echo "Completed CAZy_HMM results for ${proteins}..."
}


hmm_metabolic_fn(){
	proteome_name=$(echo "$proteins" | sed 's/.faa//')
	mkdir ../${output_dir}/tmp/${proteome_name}
	for i in $(ls -1 $metabolic_hmm_database); do
		cutoff=$(grep "$i" $metabolic_hmm_cutoffs | awk '{print $3}')
		hmmsearch --cpu 6 --tblout ../${output_dir}/tmp/${proteome_name}/${i}_out.txt -T $cutoff ${HMM_Database}/${i} ${proteins} > /dev/null 2>&1
	done
  echo "Completed Metabolic_HMM results for ${proteins}"
}


#### End Define Functions ###

mkdir -p ${output_dir}
cd ${input_dir}


## Code chunk for running prodigal
if [[ $prod_run == "T" ]]; then
  echo -e "Contig DNA as input...\nRunning prodigal on contig files with ${threads} threads"
	ls -1 *.fasta > all_samples.txt
	split all_samples.txt -l $threads batch.
	for batch in $(ls -1 batch.*); do
		for contigs in $(cat $batch); do
			prodigal_fn &
		done
		wait
	done
	rm all_samples.txt batch.*
	rename ".fasta.faa" ".faa" *.faa ###need to find a better way to rename
fi
cd ..


## Code chunk for running CheckM
if [[ $checkm_run == "T" ]]; then
  mkdir -p ${output_dir}/checkm_output
  echo -e "Running CheckM on genomes with 10 threads" ###Can modify later
  # Need pplacer in path
  PATH=/home/sdiamond/bin/pplacer-Linux-v1.1.alpha17/:$PATH ###This will need to be generalized for other users
  checkm_fn
  wait
  echo -e "CheckM Run complete..."
fi


## Code chunk for performing CAZy HMM annotation
if [[ $cazy_run == "T" ]]; then
  echo -e "Beginning CAZy HMM annotation of protein files..."
  mkdir -p ${output_dir}/CAZy_HMM_output
  cd ${input_dir}
  ls -1 *.faa > all_samples.txt
  split all_samples.txt -l $threads batch.
  for batch in $(ls -1 batch.*); do
	   for proteins in $(cat $batch); do
		    hmm_CAZy_fn &
	   done
	   wait
  done
  rm all_samples.txt batch.* *_out_domtbl.txt

  for i in $(ls -1 *.parse); do ## Here we omit CBM models and clean data
    newname=$(echo "$i" | sed 's/.faa_out_domtbl.parse//')
	   grep -e GH -e GT -e AA -e PL -e CE ${i} | sed 's/.hmm//g' > ../${output_dir}/CAZy_HMM_output/${newname}
  done
  rm *.parse
  cd ..
  echo -e "CAZy HMM annotation complete"
fi


## Code chunk for performing CAZy HMM annotation
if [[ $metabolic_run == "T" ]]; then
  echo -e "Beginning Metabolic HMM annotation of protein files..."
  mkdir -p ${output_dir}/tmp
  cd ${input_dir}
  ls -1 *.faa > all_samples.txt
  split all_samples.txt -l $threads batch.

  for batch in $(ls -1 batch.*); do
	   for proteins in $(cat $batch); do
       hmm_metabolic_fn &
     done
     wait
  done
  rm all_samples.txt batch.*

  cd ../${output_dir}/
  mkdir -p Metabolic_HMM_output
  cd tmp
  for i in $(ls -1); do
   cd $i
   touch ../../Metabolic_HMM_output/${i}
   	for tblout in $(ls -1); do
   		cat $tblout | head -n -10 | tail -n +4 | awk '{print $1, $3, $6}' >> ../../Metabolic_HMM_output/${i}
   	done
   cd ..
  done
  cd ..
  rm -r tmp
fi

echo "HMM Annotation of protein files complete..."
echo "Cleaning Up..."


## Code Chunks for cleaning up and producing metadata
## Write Metadata
cd ${input_dir}
for i in $(ls -1 *.faa); do
  genome_name=$(echo "$i" | sed 's/.faa//')
  echo "${genome_name}" >> 1.tmp
  grep -c ">" $i >> 2.tmp
  if [[ $checkm_run == "T" ]]; then
    grep "${genome_name}" ../${output_dir}/checkm_output/checkm_summary.txt | awk '{print $13, $14}' >> 3.tmp
  fi
done
paste *.tmp > ../${output_dir}/Genome_metadata.txt
rm *.tmp

# Move protein fasta files to a folder if prodigal was run
if [[ $prod_run == "T" ]]; then
  mkdir -p ../${output_dir}/Prodigal_out
  mv *.faa ../${output_dir}/Prodigal_out/
fi

#

echo -e "\n...Annotation complete..."
echo -e "HMM hit tables and alignments can be found in ${output_dir}/HMM_output/"
if [[ $checkm_run == "T" ]]; then
  echo -e "CheckM output can be found in ${output_dir}/checkm_output/"
fi
if [[ $prod_run == "T" ]]; then
  echo -e "Proteins predicted by prodigal can be found in ${output_dir}/Prodigal_out/"
fi
echo -e "A table of sample metadata has been written to ${output_dir}/Genome_metadata.txt"
