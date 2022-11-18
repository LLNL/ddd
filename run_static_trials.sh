#############################################################################################
# For static data
#############################################################################################

jobid="123456"
numtrials="4"

for ((j=1; j<$numtrials; j=j+1)) do
    echo
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp2.7/data_griewank_1000.csv 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp3/data_griewank_100000.csv 
    python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/parab_exp3/data_paraboloid_1000.csv 
done

# ls zz-temp-seed*.csv > allfiles.multi
# python generate_ddd_figures.py allfiles.multi