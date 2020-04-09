# Shell script to run analysis_cycle.C in batch mode
# to rum:
# ./multi_analysis.sh XXX YYY
# will run analysis_cycle.C from run number XXX to run number YYY inclusive
# requires that data has been converted using mvme2root, and process_rabbit has ben run

if [ "$#" == 1 ]
then
    a=${1}
    b=${1}
elif [ "$#" == 2 ]
then
    a=${1}
    b=${2}
    if (($b < $a))
    then
        b=$a
    fi
else
    echo "Requires 1 or 2 arguments"
    exit
fi

for (( i=$a; i<=$b; i++ ))
do
    echo "Run number" $i
    root -l -q "analysis_cycle.C($i)"
done
