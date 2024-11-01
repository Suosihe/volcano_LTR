#!/usr/bin/bash

usageHelp="Usage: ${0##*/}\n                   Required parameters:\n               [-l reference.fa.pass.list] [-f reference.fa]\n               [-n prefix] [-s genome_size_bp_file]\n                   Optional parameters:\n               [-c sequence identity threshold for cd-hit] default 0.8\n               [-L alignment coverage for the longer sequence in cd-hit], default 0.8\n               [-T threads for cd-hit] default 0, use all CPUs\n               [-M memory limit (in MB) for cd-hit] default 0 for no limit\n               [-p parallel number for RepeatMasker] default 80\n               [-v 'div' in RepeatMasker] Masks only those repeats with (number) percent diverged from consensus, default 40\n"
pass_listHelp="* [reference.fa.pass.list] must be defined in the script. It is in your LTR_retriver path or [EDTA.raw/LTR]. \n"
ref_faHelp="* [reference.fa] must be defined in the script. It is the genome assembly used to run LTR_retriver or EDTA. \n"
prefixHelp="* [prefix] must be defined in the script. It is the prefix you defined for your output. \n"
genome_sizeHelp="* [genome_size_bp] must be defined in the script. It is a file which only contains genome size (bp) number. It can be got from index file. \n"
helpHelp="* -h, --help: print this help and exit"


printHelpAndExit() {
	echo -e "$usageHelp"
	echo -e "$pass_listHelp"
	echo -e "$ref_faHelp"
	echo -e "$prefixHelp"
	echo -e "$genome_sizeHelp"
	echo "$helpHelp"
	exit "$1"
}

while getopts "l:f:n:s:h:c:L:T:M:p:v:" opt;
do 
    case $opt in
        l) passlist=$OPTARG
            echo "found reference.fa.pass.list: $passlist";;
        h) printHelpAndExit 0;;
        f) ref_fa=$OPTARG
            echo "found reference.fa: $ref_fa";; 
        n) prefix=$OPTARG
            echo "defined prefix: $prefix";;
        s) size_file=$OPTARG
            echo "genome size file is: $size_file";;
        c) cd_c=$OPTARG
            echo "sequence identity threshold for cd-hit is reset to $cd_c";;
        L) cd_L=$OPTARG
            echo "alignment coverage for the longer sequence in cd-hit is reset to $cd_L";;
        T) cd_T=$OPTARG
            echo "threads for cd-hit is reset to $cd_T";;
        M) cd_M=$OPTARG
            echo "memory limit (in MB) for cd-hit is reset to $cd_M";;
        p) RM_p=$OPTARG
            echo "parallel number for RepeatMasker is reset to $RM_p";;
        v) RM_div=$OPTARG
            echo "Masks only those repeats with $RM_div percent diverged from consensus";;
        -) case "${OPTARG}" in
            "help")   printHelpAndExit 0;;
            *) echo "Unknown argument --${OPTARG}";
                   printHelpAndExit 1;;
           esac;;
    [?]) printHelpAndExit 1;;
    esac
done

function get_num {
	cat $size_file
}

size=$(get_num)

if [ -z "$cd_c" ]
then
    cd_c=0.8
else
    cd_c=$cd_c
fi

if [ -z "$cd_L" ]
then
    cd_L=0.8
else
    cd_L=$cd_L
fi

if [ -z "$cd_T" ]
then
    cd_T=0
else
    cd_T=$cd_T
fi

if [ -z "$cd_M" ]
then
    cd_M=0
else
    cd_M=$cd_M
fi

if [ -z "$RM_p" ]
then
    RM_p=80
else
    RM_p=$RM_p
fi

if [ -z "$RM_div" ]
then
    RM_div=40
else
    RM_div=$RM_div
fi

echo -e "\n\nParameters have been set:"
echo -e "\nThe genome size is $size bp."
echo -e "\nThe sequence identity threshold for cd-hit is $cd_c"
echo -e "\nAlignment coverage for the longer sequence in cd-hit is $cd_L"
echo -e "\nThreads for cd-hit is $cd_T"
echo -e "\nMemory limit (in MB) for cd-hit is $cd_M"
echo -e "\nParallel number for RepeatMasker is $RM_p"
echo -e "\nMasks only those repeats with $RM_div percent diverged from consensus."

echo -e "\n\nPackage checking:\n"

check_package(){
	if which $1 >/dev/null 2>&1; then
	  echo -e "$1 installed."
	else
    	  echo "error: $1 uninstall! 'cd-hit-est', 'RepeatMasker', 'blast', 'mafft', 'fasttree' must be all installed!"
	  exit
	fi
}

check_package "cd-hit-est"
check_package "RepeatMasker"
check_package "makeblastdb"
check_package "mafft"
check_package "fasttree"


echo -e "\n\nPerl scripts checking:\n"

check_perl(){
	if which $1 >/dev/null 2>&1; then
	  echo -e "$1 found."
	else
	  echo -e "error: $1 not found! Please check that /path/to/volcano/scripts should be in your \$PATH!"
	  exit
	fi
}


check_perl list_to_5ltr.pl
check_perl obtain_lib_list.pl
check_perl fam_coverage.pl
check_perl list2_ltr_seq.pl
check_perl assign_domain_based.pl
check_perl obtain_lib_list_num.pl
check_perl re_judge2.pl
check_perl add_family_info.pl

echo -e "\nAll dependent packages are installed successfully!\n"

function get_ref {
        which volcano.sh | sed 's/\/volcano.sh//g'
}

path1=$(get_ref)
echo "ref path is $path1/dataset/"


list_to_5ltr.pl $passlist $ref_fa > ${prefix}_ltr.fa

cd-hit-est -i ${prefix}_ltr.fa -o ${prefix}_clust.out -c $cd_c -aL $cd_L -T $cd_T -M $cd_M -n 5 -d 200

obtain_lib_list.pl ${prefix}_clust.out.clstr > ${prefix}_cluster.list

RepeatMasker -e ncbi -pa $RM_p -q -no_is -norna -nolow -div $RM_div -lib ${prefix}_ltr.fa -cutoff 225 $ref_fa

fam_coverage.pl ${prefix}_ltr.fa ${ref_fa}.out $size  > ${prefix}_fam_coverage

list2_ltr_seq.pl $passlist $ref_fa $prefix

makeblastdb -in ${prefix}.ltr.fa -dbtype nucl

tblastn -query $path1/dataset/Gyp_marker.fa -db ${prefix}.ltr.fa -out ${prefix}.gyp.out -max_target_seqs 1000000000 -max_hsps 1 -evalue 10e-5 
tblastn -query $path1/dataset/copia.marker.fa -db ${prefix}.ltr.fa -out ${prefix}.cop.out -max_target_seqs 1000000000 -max_hsps 1 -evalue 10e-5

extract_RT_from_blast.pl ${prefix}.gyp.out  |sort -k 1,1 -k 3,3nr -|perl -e 'while(<>){chomp;@a=split/\t/,$_;$hash{$a[0]}++;if($hash{$a[0]}==1){print ">$a[0]\n$a[4]\n";}}' -> ${prefix}_gypsy.RT.fa
extract_RT_from_blast.pl ${prefix}.cop.out  |sort -k 1,1 -k 3,3nr -|perl -e 'while(<>){chomp;@a=split/\t/,$_;$hash{$a[0]}++;if($hash{$a[0]}==1){print ">$a[0]\n$a[4]\n";}}' -> ${prefix}_copia.RT.fa

cat ${prefix}_copia.RT.fa $path1/dataset/copia.marker.fa  > ${prefix}_copia.rt.fa
cat ${prefix}_gypsy.RT.fa $path1/dataset/Gyp_marker.fa > ${prefix}_gypsy.rt.fa 


mafft ${prefix}_copia.rt.fa > ${prefix}_copia.rt.align
fasttree -quote ${prefix}_copia.rt.align > ${prefix}_copia.rt.tree

mafft ${prefix}_gypsy.rt.fa > ${prefix}_gypsy.rt.align
fasttree -quote ${prefix}_gypsy.rt.align > ${prefix}_gypsy.rt.tree

python $path1/scripts/classify_clade.py ${prefix}_copia.rt.tree ${prefix}_copia.class.csv
python $path1/scripts/classify_clade.py ${prefix}_gypsy.rt.tree ${prefix}_gypsy.class.csv

assign_domain_based.pl $prefix.ltr.acc ${prefix}_copia.rt.fa ${prefix}_gypsy.rt.fa >${prefix}_cluster_ltr_acc_domain

obtain_lib_list_num.pl ${prefix}_clust.out.clstr >${prefix}_clust.out.clstr.list

re_judge2.pl ${prefix}_clust.out.clstr.list ${prefix}_cluster_ltr_acc_domain > ${prefix}_re_judge2.out

add_family_info.pl ${prefix}_re_judge2.out ${prefix}_fam_coverage > ${prefix}_temp0

cat ${prefix}_copia.class.csv ${prefix}_gypsy.class.csv|grep -v "Digit"|perl -pe 's/\r\n/\n/g'|awk -F"," '{print $1"\t"$2}'|sort -k2,2 > temp1
join -t $'\t' -1 2 temp1 $path1/dataset/class.std.txt| awk '{print $2,$1,$3}' |sort -k1,1|sed 's/ /\t/g' >temp2
grep -v "#RepeatMasker_entry" ${prefix}_temp0|sort -k17,17 | join -t $'\t' -a 1 -e 'NULL' -1 17 - temp2 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,2.3' > temp3
sed '1s/^/#RepeatMasker_entry\tTE_family\tFull_length\tLeft_end_only\tRight_end_only\tConverted_copy_number\tTotal_entries\tTotal_length_in_bp\tWhole_genome_percentage\tClass\tSubclass\tNote\tFamily\tNUM\tType\tTREE_type\tID\tClade\n/' temp3 > ${prefix}_fam_coverage.info

cut -f 1,3 temp2 |grep "Ty3" |sed '1i ID\tClade' > ${prefix}_gypsy_clade.tsv
cut -f 1,3 temp2 |grep "Ty1" |sed '1i ID\tClade' > ${prefix}_copia_clade.tsv

rm -rf ${prefix}_temp0 temp1 temp2 temp3 ${prefix}_fam_coverage ${prefix}_copia.rt.align ${prefix}_gypsy.rt.align ${prefix}_copia.rt.fa ${prefix}_gypsy.rt.fa ${prefix}_gypsy.RT.fa ${prefix}_copia.RT.fa ${prefix}_copia.class.csv ${prefix}_gypsy.class.csv ${prefix}_clust.out ${prefix}.clust.num ${prefix}.ltr* ${prefix}_cluster_ltr_acc_domain ${ref_fa}.cat.gz ${ref_fa}.out ${ref_fa}.masked ${ref_fa}.tbl ${prefix}.gyp.out ${prefix}.cop.out ${prefix}_re_judge2.out
