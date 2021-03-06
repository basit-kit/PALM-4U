#!/bin/ksh


echo "\$1 = " $1 

if [[ $1 == "clean" ]]
then
  cd kpp
  make distclean
  cd ..

  cd kp4
  make distclean                 
  cd ..
  exit
fi



export PATH=$PATH:`pwd`/kp4/bin

# build kpp

cd kpp
make
cd ..

# build kp4.exe

cd kp4
make install                    
cd ..

# run kp4 with default Setup
#    -d   DEFDIR=$OPTARG;;          # directory of definition files
#    -i   DE_INDEX="YES";;          # if set, deindexing
#    -f   DE_INDEX_FAST="YES";;     # if set, fast deindexing
#    -k   KEEP="YES";;              # keep Working directory
#    -o   OUTDIR=$OPTARG;;          # Output directory of Geneated Code
#    -p   PREFIX=$OPTARG;;          # Name Prefix
#    -s   KPP_SOLVER=$OPTARG;;      # Name Prefix

 
echo $PATH
#kp4.ksh -k "YES" -d `pwd`/kp4/def_small_strato
#kp4.ksh -d `pwd`/kp4/def_smog -k "YES"
#kp4.ksh -d `pwd`/kp4/def_tests -k "YES"
MECH=smog
while  getopts m:  opt   # get options
do case $opt in
      m)   MECH=$OPTARG;;          # mechanism name as part of kp4/def_[mechanism_name]

      \?)  print ${0##*/} "unknown option:" $OPTARG
           print "USAGE: ${0##*/} [ -m mechanism_name ] "
           exit 1;;
   esac
done
echo $MECH
kp4.ksh -d `pwd`/kp4/def_$MECH -k "YES" 
