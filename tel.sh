#!/usr/bin/bash

usageHelp="Usage: ${0##*/}     Required parameters:\n               [-g pass.list.gff3]  [-s sorted.sam]  [-p prefix]  [-n number]\n"
gff3Help="* [pass.list.gff3] must be defined in the script. It is in your LTR_retriver path or [EDTA.raw/LTR/]. \n"
samHelp="* [sorted.sam] must be defined in the script. It is a sorted sam file usually contains alignment results. \n"
prefixHelp="* [prefix] must be defined in the script. It is the prefix you defined for your output. \n"
numberHelp="* [number, genome_size_bp] must be defined in the script. It is a file which only contains genome size (bp) number. It can be got from index file. \n"
helpHelp="* -h, --help: print this help and exit"

function get_ref {
        which volcano0531.sh | sed 's/\/volcano0531.sh//g'
}

path1=$(get_ref)

printHelpAndExit() {
        echo -e "$usageHelp"
        echo -e "$gff3Help"
        echo -e "$samHelp"
        echo -e "$prefixHelp"
        echo -e "$numberHelp"
        echo "$helpHelp"
        exit "$1"
}

while getopts "g:s:p:h:n:" opt;
do
    case $opt in
        g) gff=$OPTARG
            echo "found pass.list.gff3: $gff";;
        h) printHelpAndExit 0;;
        s) sorted_sam=$OPTARG
            echo "found sorted.sam: $sorted_sam";;
        p) prefix=$OPTARG
            echo "prefix defined: $prefix";;
        n) num=$OPTARG
            echo "genome size file is: $num";;
        -) case "${OPTARG}" in
            "help")   printHelpAndExit 0;;
            *) echo "Unknown argument --${OPTARG}";
                   printHelpAndExit 1;;
           esac;;
    [?]) printHelpAndExit 1;;
    esac
done

grep "LTR_retrotransposon" $gff > $prefix.LTR.gff3
sed -i 's/ID/locus/g' $prefix.LTR.gff3
awk -f $path1/bin/gff2gtf.awk $prefix.LTR.gff3 >$prefix.LTR.gtf

telescope assign --exp_tag $prefix $sorted_sam $prefix.LTR.gtf 
Rscript $path1/bin/RPKM.R -i $prefix-telescope_report.tsv -n $num -o $prefix
