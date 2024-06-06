# $1 dir_name;$2 EDTA_PATH;$3 EDTAfasta

usageHelp="Usage: ${0##*/} [-d dir_prefix] [-E /path/to/EDTA e.g. ~/A_thaliana/EDTA]\n       [-f /path/to/fasta(for EDTA) e.g. ~/A_thaliana/EDTA/A_thaliana.fa.mod]\n"
dirHelp="[dir_prefix] must be defined in the script\n"
pathEDTAHelp="[/path/to/EDTA] not found\n"
famodHelp="[/path/to/fasta(for EDTA)] not found\n"
helpHelp="* -h, --help: print this help and exit"

printHelpAndExit() {
         echo -e "$usageHelp"
	 echo -e "$dirHelp"
	 echo -e "$pathEDTAHelp"
	 echo -e "$famodHelp"
	 echo "$helpHelp"
	 exit "$1"
}

while getopts "d:h:E:f:" opt;
do
    case $opt in
        d) dir_prefix=$OPTARG
            echo "dir name defined: $dir_prefix";;
        E) path_EDTA=$OPTARG
	    echo "EDTA path found: $path_EDTA";;
        f) MOD_fa=$OPTARG
	    echo "fasta for EDTA found: $MOD_fa";;
        -) case "${OPTARG}" in
            "help")   printHelpAndExit 0;;
        h) printHelpAndExit 0;;
        *) echo "Unknown argument --${OPTARG}";
             printHelpAndExit 1;;
           esac;;
    [?]) printHelpAndExit 1;;
    esac
done


# prepare the documents
mkdir $dir_prefix
cd $dir_prefix
cp $path_EDTA/*EDTA.raw/LTR/*.pass.list .; rm -rf *nmtf.pass.list;mv *.pass.list pass.list
cp $path_EDTA/$MOD_fa fasta

# find the length of fasta
samtools faidx fasta
tail -n 1 *fai|cut -f 3 > len
rm -rf *fai

