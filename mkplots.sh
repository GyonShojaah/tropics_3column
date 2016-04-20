sed -e "s/filetag/$1/g" gp_all | gnuplot
mkdir plots_$1
mv *.eps plots_$1
