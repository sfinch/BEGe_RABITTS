
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
    ./process_rabbit $i
    root -l -q "plot_nmon.C($i)"
    root -l -q "run_info.C($i)"
done
