#############################################################################################
# For static data
#############################################################################################

jobid="123456"
numtrials="10"

for ((j=1; j<$numtrials+1; j=j+1)) do
    echo
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp2.7/data_griewank_1000.csv 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/grwk_exp3/data_griewank_100000.csv 
    python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/parab_exp3/data_paraboloid_1000.csv
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/parab_exp3/data_paraboloid_100000.csv --validsplit 0.0005 --itmax 6 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/grwk_exp3/data_griewank_100000.csv --validsplit 0.0005 --initsampprop 0.002--itmax 6 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/grwk_exp3/data_griewank_100000.csv --validsplit 0.005 --initsampprop 0.02 --itmax 4 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/parab_exp3/data_paraboloid_10000.csv --validsplit 0.05 --initsampprop 0.2 --itmax 3
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/grwk_exp3/data_griewank_10000.csv --validsplit 0.005 --initsampprop 0.2 --itmax 3
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/from_ball/grwk_exp3/data_griewank_100000.csv 
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --zoomexp 3
done

ls zz-temp-seed*.csv > allfiles.multi
python generate_ddd_figures.py allfiles.multi