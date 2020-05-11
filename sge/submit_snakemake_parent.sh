sbatch -J ten-sn -p genx --mem=64G --time=3-20:0:0 -n 10 --ntasks-per-node 10 --wrap "disBatch.py cmd_list.txt"
