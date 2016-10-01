#!/bin/bash

usage(){
      echo -e "\nUsage:  CAZY_hmmscan_files.sh -i <input_directory> -o <output_directory> [options]"
      echo -e "\n\tGiven a directory of fasta format protein or DNA files this script will run hmmscan against"
      echo -e "\ta database of CAZy HMMs derived from dbCAN. Output is a domain table parsed for score"
      echo -e "\n\tOptions:"
      echo -e "\t-i: Directory of fasta format files to use as input (required)"
      echo -e "\t-o: Directory where output files will be written (required, will be created if absent)"
      echo -e "\t-p: Input directory contains DNA contigs; Run prodigal in single genome mode to predict proteins"
      echo -e "\t-c: Do not run CheckM for genome completeness (Runs by Default)"
      echo -e "\t-t: Number of threads to use (Default: 6)"
      echo -e "\t-e: E-value cutoff for HMM scoring (Default: 10)"
      echo -e "\t-h: Display this message and exit"
      echo -e "\nby: Spencer Diamond September 9, 2016"
}

## Default Arguments ##
threads=6
evalue=10
checkm_run=T
## Default Arguments ##

while getopts ":hi:o:pct:e:" opt; do
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


hmm_scan_fn(){ ###Can modify this function to accept modular HMM database inputs possibly use if statements with set up databases
	hmmscan --domtblout ${proteins}_out_domtbl.txt --cpu 6 --domE ${evalue} /home/sdiamond/database/CAZy_HMM/dbCAN-fam-HMMs.txt ${proteins} > ${proteins}_out_full.aln
	wait
	echo "Completed hmmscan for ${proteins}..."
	sh /home/sdiamond/database/CAZy_HMM/hmmscan-parser.sh ${proteins}_out_domtbl.txt > ${proteins}_out_domtbl.parse
	wait
	echo "Completed parsing results for ${proteins}..."
}


#### End Define Functions ###

mkdir -p ${output_dir}
cd ${input_dir}

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


if [[ $checkm_run == "T" ]]; then
  mkdir -p ${output_dir}/checkm_output
  echo -e "Running CheckM on genomes with 10 threads" ###Can modify later
  # Need pplacer in path
  PATH=/home/sdiamond/bin/pplacer-Linux-v1.1.alpha17/:$PATH ###This will need to be generalized for other users
  checkm_fn
  wait
fi

mkdir -p ${output_dir}/HMM_output
cd ${input_dir}
ls -1 *.faa > all_samples.txt
split all_samples.txt -l $threads batch.
for batch in $(ls -1 batch.*); do
	for proteins in $(cat $batch); do
		hmm_scan_fn &
	done
	wait
done
rm all_samples.txt batch.* *_out_domtbl.txt
mv *_out_full.aln ../${output_dir}/HMM_output/
mv *_out_domtbl.parse ../${output_dir}/HMM_output/
cd ..
echo "...HMM Scan of protein files complete..."
