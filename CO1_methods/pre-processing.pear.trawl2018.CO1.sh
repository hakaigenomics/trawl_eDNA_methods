#notes for pre-processing data from the trawl 2018 data

#copy data from NAS
rsync -aP /mnt/MiSeq/CO1_Run1_Training_03062020/Run03062020/Alignment_1/20200306_010319/Fastq/* ~/projects/CO1_runs/raw_lanes/CO1_Run1_Training_03062020/raw_data/
#and for other datasets
#CO1_ReRun2_04032020  CO1_Run1_GenomeQuebec_Oct21,2019  CO1_Run1_GenomeQuebec_Sept23,2019  CO1_Run1_Training_03062020

#now organize a folder for this experiment "trawl-2018" with subfolders for raw data, metadata, etc.
#symlink raw data in its raw_lane folder to the raw data folder for trawl-2018
ln -s /home/evan.morien/projects/CO1_runs/raw_lanes/CO1_Run1_GenomeQuebec_Sept23,2019/* ./ #and for other datasets (so now we have four sub folders in raw_data corresponding to the four lanes listed above)

#now demultiplex each one according to which samples were run in which lane
#the files are all "double" multiplexed. they have been demultiplexed into "pools" prior to receiving from illumina. in each pool there are a set of barcodes (all shared between the pools, but unique w/in pools) that distinguish samples.

#to use fastx to demultiplex (barcodes are in-read, no index file) you need 1) a barcode sheet 2) input fastq files. with a list of fastq files and the barcode sheet for the entire run, and a list of pools to iterate through, we can automatically demultiplex the run
#make barcode sheet: take from .xls format is 4 columns, tab separated, this order: sample_name, BarcodeSequence, Tag, Pool
#metadata/barcode_sheet.EXPT.RUNX.txt

#reads need to be stitched (merged) with PEAR, then properly demultiplexed, then properly oriented, then primer trimmed with cutadapt
#working dir: ~/projects/CO1_runs/trawl-2018/raw_data #
~/programs/microbiome_helper/run_pear.pl -g -p 36 -o raw_data/CO1_ReRun2_04032020/stitched_reads/ raw_data/CO1_ReRun2_04032020/*.gz
~/programs/microbiome_helper/run_pear.pl -g -p 36 -o raw_data/CO1_Run1_GenomeQuebec_Oct21,2019/stitched_reads/ raw_data/CO1_Run1_GenomeQuebec_Oct21,2019/*.gz
~/programs/microbiome_helper/run_pear.pl -g -p 36 -o raw_data/CO1_Run1_GenomeQuebec_Sept23,2019/stitched_reads/ raw_data/CO1_Run1_GenomeQuebec_Sept23,2019/*.gz
~/programs/microbiome_helper/run_pear.pl -g -p 36 -o raw_data/CO1_Run1_Training_03062020/stitched_reads/ raw_data/CO1_Run1_Training_03062020/*.gz

#one pool didn't assemble well (we retained only 14% of reads after merging): Run2-Pool9 #i checked and this is the negative control, so that's expected

#make pool list
awk '{print $4}' metadata/barcode_sheet.trawl-2018.run1.txt | sort | uniq > metadata/pools.run1.txt #important: remove the header label from this unique list of pools, or remove the header column from the input file. either way, you only want pool names in this file, nothing else.
awk '{print $4}' metadata/barcode_sheet.trawl-2018.run2.txt | sort | uniq > metadata/pools.run2.txt

#make file list #IMPORTANT: exact matches between pool names in files and metadata MUST be the only kind of matches. i.e. no match for pool03 in metadata and Pool3 in filelist, and Pool3 in filelist shouldnt be allowed to match Pool3 AND Pool31, but capitalization isnt (shouldnt be!) important
find raw_data/CO1_ReRun2_04032020/stitched_reads/*assembled.fastq.gz -maxdepth 1 | grep -i pool | sort > filelist.CO1_ReRun2_04032020.txt
find raw_data/CO1_Run1_GenomeQuebec_Oct21,2019/stitched_reads/*assembled.fastq.gz -max depth 1 | grep -i pool | sort > filelist.CO1_Run1_GenomeQuebec_Oct21,2019.txt
find raw_data/CO1_Run1_GenomeQuebec_Sept23,2019/stitched_reads/*assembled.fastq.gz -maxdepth 1 | grep -i pool | sort > filelist.CO1_Run1_GenomeQuebec_Sept23,2019.txt
find raw_data/CO1_Run1_Training_03062020/stitched_reads/*assembled.fastq.gz -maxdepth 1 | grep -i pool | sort > filelist.CO1_Run1_Training_03062020.txt

#cleanup during testing
#rm demultiplexed/Pool*
#rm filelist.Pool*
#rm ../metadata/barcode_sheet.Pool*
#mkdir demultiplexed
#some lanes need to be combined before analysis, based on sampleID, so we need separate subfolders initially in the demultiplexed folder
#mkdir demultiplexed/CO1_ReRun2_04032020
#mkdir demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019
#mkdir demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019
#mkdir demultiplexed/CO1_Run1_Training_03062020

#here are the barcode file / lane combos #the loop below needs to be run for each combo. consider automating in future.
barcode_sheet.trawl-2018.run1.txt - CO1_Run1_GenomeQuebec_Oct21,2019 & CO1_Run1_GenomeQuebec_Sept23,2019 & CO1_Run1_Training_03062020
barcode_sheet.trawl-2018.run2.txt - CO1_ReRun2_04032020

#currently needs changing per loop run: filelist.X.txt, barcode_sheet.X.txt, pools.X.txt, demultiplexing output folders, fastq-multx.X.log
#iterate over list
#order of operations here is important. you can't just do cutadapt first, it makes it so you can't demultiplex after because we lose the barcodes in the process. we need to demultiplex, then run cutadapt in a second loop.
#the nomenclature from fastx toolkit makes it more difficult to run cutadapt in the same loop as fastx-mult. the % referent is the issue, and i don't feel like fussing around to fix that.
#will do a second loop for cutadapt below.

######ReRun2 (run2)######
while read line; do
        echo "######################################################### demultiplexing "$line;
        grep -w $line metadata/barcode_sheet.trawl-2018.run2.txt | awk '{print $1"\t"$2"\t"$4}' > metadata/sheets_tmp/barcode_sheet.${line}.txt;
        grep -F ${line}_ filelist.CO1_ReRun2_04032020.txt > filelist_tmp/filelist.${line}.txt;
        num_files=`grep assembled filelist_tmp/filelist.${line}.txt | wc -l`;
        #num_samp=$((num_files / 2))
        #echo $num_files;
        START=1;
        END=$num_files;
        echo "processing samples for pool ${line}";
        if [ "$END" -eq "$START" ]; then
                infile=`grep assembled filelist_tmp/filelist.${line}.txt | head -n 1`; #doing this accounts for cases where we're missing an R1 or R2 file
                /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_ReRun2_04032020/${line}.%.fastq.gz;
        else
                for i in $(seq $START $END); do
                        echo "demultiplexing" $line "run" $i "of" $END
                        infile=`grep assembled filelist_tmp/filelist.${line}.txt| head -n ${i} | tail -1`;
                        /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_ReRun2_04032020/${line}.%.run${i}.fastq.gz;
                done
        fi
        echo "#########################################################";
done < metadata/pools.run2.txt &>demultiplexed/fastq-multx.CO1_ReRun2_04032020.pear_merged.log

mkdir demultiplexed/CO1_ReRun2_04032020/reorient_primer_trim
ls demultiplexed/CO1_ReRun2_04032020/*.fastq.gz | grep -v unmatched | awk -F"/" '{print $3}' > filelist.reorient_primer_trim.txt
while read line; do
cutadapt --revcomp -g GGWACWGGWTGAACWGTWTAYCCYCC -j 36 -O 15 demultiplexed/CO1_ReRun2_04032020/${line} | cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTYGGNCAYCCNGARGTNTA -n 2 -j 36 -o demultiplexed/CO1_ReRun2_04032020/reorient_primer_trim/${line} -
done < filelist.reorient_primer_trim.txt &>>demultiplexed/fastq-multx.CO1_ReRun2_04032020.pear_merged.log

#it works!

######GenomeQuebec October (run1)######
while read line; do
        echo "######################################################### demultiplexing "$line;
        grep -w $line metadata/barcode_sheet.trawl-2018.run1.txt | awk '{print $1"\t"$2"\t"$4}' > metadata/sheets_tmp/barcode_sheet.${line}.txt;
        grep -F ${line}\. filelist.CO1_Run1_GenomeQuebec_Oct21,2019.txt > filelist_tmp/filelist.${line}.txt; #needs a period not an underscore at the end of the grep pattern
        num_files=`grep assembled filelist_tmp/filelist.${line}.txt | wc -l`;
        #num_samp=$((num_files / 2))
        #echo $num_files;
        START=1;
        END=$num_files;
        echo "processing samples for pool ${line}";
        if [ "$END" -eq "$START" ]; then
                infile=`grep assembled filelist_tmp/filelist.${line}.txt | head -n 1`; #doing this accounts for cases where we're missing an R1 or R2 file
                /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/${line}.%.fastq.gz;
        else
                for i in $(seq $START $END); do
                        echo "demultiplexing" $line "run" $i "of" $END
                        infile=`grep assembled filelist_tmp/filelist.${line}.txt| head -n ${i} | tail -1`;
                        /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/${line}.%.run${i}.fastq.gz;
                done
        fi
        echo "#########################################################";
done < metadata/pools.run1.txt &>demultiplexed/fastq-multx.CO1_Run1_GenomeQuebec_Oct21,2019.pear_merged.log

mkdir demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/reorient_primer_trim
ls demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/*.fastq.gz | grep -v unmatched | awk -F"/" '{print $3}' > filelist.reorient_primer_trim.txt
while read line; do
cutadapt --revcomp -g GGWACWGGWTGAACWGTWTAYCCYCC -j 36 -O 15 demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/${line} | cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTYGGNCAYCCNGARGTNTA -n 2 -j 36 -o demultiplexed/CO1_Run1_GenomeQuebec_Oct21,2019/reorient_primer_trim/${line} -
done < filelist.reorient_primer_trim.txt &>>demultiplexed/fastq-multx.CO1_Run1_GenomeQuebec_Oct21,2019.pear_merged.log

######GenomeQuebec September (run1)######
while read line; do
        echo "######################################################### demultiplexing "$line;
        grep -w $line metadata/barcode_sheet.trawl-2018.run1.txt | awk '{print $1"\t"$2"\t"$4}' > metadata/sheets_tmp/barcode_sheet.${line}.txt;
        grep -F ${line}\. filelist.CO1_Run1_GenomeQuebec_Sept23,2019.txt > filelist_tmp/filelist.${line}.txt; #needs a period not an underscore at the end of the grep pattern
        num_files=`grep assembled filelist_tmp/filelist.${line}.txt | wc -l`;
        #num_samp=$((num_files / 2))
        #echo $num_files;
        START=1;
        END=$num_files;
        echo "processing samples for pool ${line}";
        if [ "$END" -eq "$START" ]; then
                infile=`grep assembled filelist_tmp/filelist.${line}.txt | head -n 1`; #doing this accounts for cases where we're missing an R1 or R2 file
                /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/${line}.%.fastq.gz;
        else
                for i in $(seq $START $END); do
                        echo "demultiplexing" $line "run" $i "of" $END
                        infile=`grep assembled filelist_tmp/filelist.${line}.txt| head -n ${i} | tail -1`;
                        /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/${line}.%.run${i}.fastq.gz;
                done
        fi
        echo "#########################################################";
done < metadata/pools.run1.txt &>demultiplexed/fastq-multx.CO1_Run1_GenomeQuebec_Sept23,2019.pear_merged.log

mkdir demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/reorient_primer_trim
ls demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/*.fastq.gz | grep -v unmatched | awk -F"/" '{print $3}' > filelist.reorient_primer_trim.txt
while read line; do
cutadapt --revcomp -g GGWACWGGWTGAACWGTWTAYCCYCC -j 36 -O 15 demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/${line} | cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTYGGNCAYCCNGARGTNTA -n 2 -j 36 -o demultiplexed/CO1_Run1_GenomeQuebec_Sept23,2019/reorient_primer_trim/${line} -
done < filelist.reorient_primer_trim.txt &>>demultiplexed/fastq-multx.CO1_Run1_GenomeQuebec_Sept23,2019.pear_merged.log


######Training (run1)######
while read line; do
        echo "######################################################### demultiplexing "$line;
        grep -w $line metadata/barcode_sheet.trawl-2018.run1.txt | awk '{print $1"\t"$2"\t"$4}' > metadata/sheets_tmp/barcode_sheet.${line}.txt;
        grep -F ${line}_ filelist.CO1_Run1_Training_03062020.txt > filelist_tmp/filelist.${line}.txt; #needs an underscore after the pool number ($line) in grep pattern
        num_files=`grep assembled filelist_tmp/filelist.${line}.txt | wc -l`;
        #num_samp=$((num_files / 2))
        #echo $num_files;
        START=1;
        END=$num_files;
        echo "processing samples for pool ${line}";
        if [ "$END" -eq "$START" ]; then
                infile=`grep assembled filelist_tmp/filelist.${line}.txt | head -n 1`; #doing this accounts for cases where we're missing an R1 or R2 file
                /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_Training_03062020/${line}.%.fastq.gz;
        else
                for i in $(seq $START $END); do
                        echo "demultiplexing" $line "run" $i "of" $END
                        infile=`grep assembled filelist_tmp/filelist.${line}.txt| head -n ${i} | tail -1`;
                        /home/evan.morien/programs/ea-utils/clipper/fastq-multx -l metadata/sheets_tmp/barcode_sheet.${line}.txt $infile -o demultiplexed/CO1_Run1_Training_03062020/${line}.%.run${i}.fastq.gz;
                done
        fi
        echo "#########################################################";
done < metadata/pools.run1.txt &>demultiplexed/fastq-multx.CO1_Run1_Training_03062020.pear_merged.log

mkdir demultiplexed/CO1_Run1_Training_03062020/reorient_primer_trim
ls demultiplexed/CO1_Run1_Training_03062020/*.fastq.gz | grep -v unmatched | awk -F"/" '{print $3}' > filelist.reorient_primer_trim.txt
while read line; do
cutadapt --revcomp -g GGWACWGGWTGAACWGTWTAYCCYCC -j 36 -O 15 demultiplexed/CO1_Run1_Training_03062020/${line} | cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC -a TGRTTYTTYGGNCAYCCNGARGTNTA -n 2 -j 36 -o demultiplexed/CO1_Run1_Training_03062020/reorient_primer_trim/${line} -
done < filelist.reorient_primer_trim.txt &>>demultiplexed/fastq-multx.CO1_Run1_Training_03062020.pear_merged.log


#now creating combined fastq files for the three runs of the same data
#do it the dumb way, no more sophisticated than it needs to be
#feed this loop a list of sampleIDs, use grep to get all the files that match that sampleID, then use zcat to concatenate them into one file. if you are doing this with separate R1 and R2 files, you need to have separate grep and zcat commands for them
#samples that have no R2 file in the GQ runs: Pool3.tdna044_Run1 (left out of analysis, since merging was not possible)
ls demultiplexed/CO1_R*/reorient_primer_trim/*assembled* |grep -v Pool3.tdna044_Run1 | grep -v unmatched > filelist_tmp/samples.all.txt
awk -F'/' '{print $4}' filelist_tmp/samples.all.txt | sed 's/.fastq.gz//' | sort | uniq > sampleIDs.all.txt
while read line; do
        echo "######################################################### concatenating sample: "$line;
        grep $line filelist_tmp/samples.all.txt > filelist_tmp/filelist.${line}.txt;
        xargs -d'\n' cat <filelist_tmp/filelist.${line}.txt >demultiplexed/merged/${line}.fastq.gz
        echo "#########################################################";
done < sampleIDs.all.txt &>demultiplexed/concatenation.log


#looks good, now remove reads with Ns using vsearch
# Define variables
ls demultiplexed/merged/*.fastq.gz > filelist.merged.txt
FILELIST=filelist.merged.txt
VTHREADS=36
MINLENGTH=100
# Define binaries, temporary files and output files
VSEARCH=$(which vsearch)
TMP_FASTQ=$(mktemp)
OUTPUT=$(mktemp)
LOG="demultiplexed/merged/vsearch.log"

while read line ; do    
    FINAL_FASTA="${line}.fas"
    # Discard sequences containing Ns, discard too-short sequences, add expected error rates
    "${VSEARCH}" \
        --quiet \
        --fastq_filter "${line}" \
        --fastq_minlen "${MINLENGTH}" \
        --fastq_maxns 0 \
        --relabel_sha1 \
        --fastqout "${TMP_FASTQ}" 2>> "${LOG}"
done < "${FILELIST}"

#now these demultiplexed, primer-trimmed, combined samples will be processed in dada2 pipeline. code located here: processing.dada2.trawl2018.CO1.R