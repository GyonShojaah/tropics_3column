import subprocess

cmd = 'Cl_run_cdf -s $RAD_DATA/spectra/ga7/sp_lw_ga7 -R 1 9 -I -g 2 -c -t 12 -v 13 -C 5 -B '+BASE_NAME[pool]
print cmd
subprocess.call( cmd, shell=True )
