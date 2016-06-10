sed -e "s/filetag/$1/g" gp_all | gnuplot

if [ ! -e "plots_"$1 ]
then
    mkdir plots_$1
fi

mv *.eps plots_$1
