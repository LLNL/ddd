#############################################################################################
# For static data
#############################################################################################

jobid="123456"
numtrials="10"

for ((j=1; j<$numtrials; j=j+1)) do
    echo
    echo Starting trial with zoom exponent =  $(bc<<<"$zoomstep * $i"), seed = $j
    echo
    # python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/data_griewank_1000.csv 
    python delaunay_density_diagnostic.py --jobid ${jobid} --seed ${j} --staticdata staticdata/data_paraboloid_100000.csv 
done

ls zz-temp-seed*.csv > allfiles.multi
python generate_ddd_figures.py allfiles.multi