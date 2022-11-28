#############################################################################################
# For static data
#############################################################################################

jobid="123456"
numtrials="9"

for ((j=1; j<$numtrials+1; j=j+1)) do
    echo
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp2.7/data_griewank_1000.csv 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp3/data_griewank_100000.csv 
    python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_gauss/parab_exp3/data_paraboloid_100000.csv --validsplit 0.1
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_gauss/grwk_exp3/data_griewank_100000.csv 
done

ls zz-temp-seed*.csv > allfiles.multi
python generate_ddd_figures.py allfiles.multi