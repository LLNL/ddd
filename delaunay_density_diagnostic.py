################################################################################
#
# Delaunay density diangostic for MSD and grad-MSD rates
#   as described in the paper
#   Algorithm XXXX: The Delaunay Density Diagnostic
#       under review at ACM Transactions on Mathematical Software
#       (original title: Data-driven geometric scale detection via Delaunay interpolation)
#   by Andrew Gillette and Eugene Kur
#   Version 2.0, November 2023
#
# For usage information, run:
# python delaunay_density_diagnostic.py --help
#
################################################################################



#==================================================================================================#
# Load packages.  Set random state and validation split.
#==================================================================================================#

# from matplotlib.pyplot import legend
# import torch
import pandas as pd     
# from torch.autograd import Variable
# import torch.nn.functional as F
# import torch.utils.data as Data
# from torch.utils.data.sampler import SubsetRandomSampler
import numpy as np

from numpy.random import rand, default_rng
from numpy import arccos, array, degrees, absolute
from numpy.linalg import norm

from optparse import OptionParser

# import numpy.ma as ma

# import xarray as xr

from sys import exit
import os.path

import copy

import delsparse
from delsparse import delaunaysparsep as dsp


#==================================================================================================#
# Define the test function (hard coded here as the Griewank function)
#==================================================================================================#

def tf(X): # Griewank function, arbitrary dimension input
    X = X.T

    # compute and return values from the Griewank function
    term_1 = (1. / 4000.) * sum(X ** 2)
    term_2 = 1.0
    for i, x in enumerate(X):
        term_2 *= np.cos(x) / np.sqrt(i + 1)
    return 1. + term_1 - term_2

    # # instead, return values from a paraboloid:
    # return (7/20_000) *  ( X[0]**2 + 0.5*(X[1]**2) )

#==================================================================================================#
# Make query point lattice in R^dim
#==================================================================================================#
def make_test_data_grid(rng, static_data=False):

    num_samples_per_dim = options.numtestperdim

    x = np.linspace(options.queryleftbound, options.queryrightbound, num_samples_per_dim)

    print("===> Test coordinates for each dimension = ", x)
    mg_in = []
    for i in range(options.dim):
        mg_in.append(x)
    grid_pts = np.array(np.meshgrid(*mg_in))

    grid_pts = grid_pts.reshape(options.dim, num_samples_per_dim ** options.dim)
    grid_pts = grid_pts.T

    outputs_on_grid = tf(grid_pts)

    data_test_inputs  = pd.DataFrame(grid_pts)
    data_test_outputs = pd.DataFrame(outputs_on_grid)

    return data_test_inputs, data_test_outputs



#==================================================================================================#
# Collect random sample from bounding box
#==================================================================================================#
def make_random_training_in_box(rng):

    train_set_size = options.numtrainpts
    # print("==> Generating ", train_set_size, " random points.")
       
    rand_pts_n = rng.random((train_set_size, options.dim))
        
    train_box_scale_vector = np.full(options.dim, (options.bboxrightbound - options.bboxleftbound) )
    train_box_shift_vector = np.full(options.dim, options.bboxleftbound )

    # do scaling in each dim first
    for i in range(options.dim):
        rand_pts_n[:,i] *= train_box_scale_vector[i]
    # then do shifts
    for i in range(options.dim):
        rand_pts_n[:,i] += train_box_shift_vector[i]

    outputs_on_rand_n = tf(rand_pts_n)

    data_train_inputs  = pd.DataFrame(rand_pts_n)
    data_train_outputs = pd.DataFrame(outputs_on_rand_n)

    return data_train_inputs, data_train_outputs



#==================================================================================================#
# Function to compute DelaunaySparse 
#==================================================================================================#
def compute_DS_only(data_train_inputs, data_train_outputs, data_test_inputs, data_test_outputs):

    # # note: data_test_outputs is only converted to numpy if needed, transposed, 
    #         and returned as actual_test_vals

    # # WARNING: deepcopy here may be inefficient at scale
    pts_in = copy.deepcopy(data_train_inputs)
    q      = copy.deepcopy(data_test_inputs)

    interp_in = data_train_outputs
    
    actual_test_vals = data_test_outputs
    if not isinstance(actual_test_vals, np.ndarray):
        actual_test_vals = actual_test_vals.to_numpy()
    actual_test_vals = actual_test_vals.T

    actual_train_vals = data_train_outputs
    if not isinstance(actual_train_vals, np.ndarray):
        actual_train_vals = actual_train_vals.to_numpy()
    actual_train_vals = actual_train_vals.T

    interp_in_n = interp_in

    if not isinstance(interp_in_n, np.ndarray):
        interp_in_n = interp_in_n.to_numpy()
    interp_in_n = interp_in_n.T

    if not isinstance(pts_in, np.ndarray):
        pts_in = pts_in.to_numpy()
    pts_in = pts_in.T
    pts_in = np.require(pts_in, dtype=np.float64, requirements=['F'])

    if not isinstance(q, np.ndarray):
        q = q.to_numpy()
    p_in = np.asarray(q.T, dtype=np.float64, order="F")

    ir=interp_in_n.shape[0]

    interp_in_n = np.require(interp_in_n, 
                dtype=np.float64, requirements=['F'])
    simp_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                dtype=np.int32, order="F")
    weights_out = np.ones(shape=(p_in.shape[0]+1, p_in.shape[1]), 
                dtype=np.float64, order="F")
    error_out = np.ones(shape=(p_in.shape[1],), 
                dtype=np.int32, order="F")
    interp_out_n = np.zeros([interp_in_n.shape[0],p_in.shape[1]])
    interp_out_n = np.require(interp_out_n, 
                dtype=np.float64, requirements=['F'])
    rnorm_n = np.zeros(p_in.shape[1])
    rnorm_n = np.require(rnorm_n, dtype=np.float64, requirements=['F'])

    # From delsparse.py documenation:
    #   Setting EXTRAP=0 will cause all extrapolation points to be
    #   ignored without ever computing a projection. By default, EXTRAP=0.1
    #   (extrapolate by up to 10% of the diameter of PTS).
    dsp(pts_in.shape[0], pts_in.shape[1],
        pts_in, p_in.shape[1], p_in, simp_out,
        weights_out, error_out, 
        extrap=options.extrap_thresh,
        rnorm=rnorm_n,
        pmode=1,
        interp_in=interp_in_n, interp_out=interp_out_n)


    if (options.computeGrad):
            
        # # arbitrary number of outputs, as determind by interp_in_n.shape[0]
        grad_est_DS = np.zeros([interp_in_n.shape[0], simp_out.shape[1], options.dim])
        grad_est_DS.fill(999)
        
        for j in range(simp_out.shape[1]):
            # note: the value of simp_out.shape[1] should equal the number of interpolation outputs
            #       extrapolation points don't get a simp_out entry, I think?
            #
            #       mutliple test points may lie in the same simplex
            #       but that just means you might duplicate effort
            #       if you already saw a simplex and comptued the gradient(s)

            # this presumes pts_in was deep copied from data_train_inputs at start of compute_DS_only(...)
            unscaled_inputs = data_train_inputs.to_numpy().T
            # can try using scaled points instead:
            # # unscaled_inputs = pts_in

            for outputdim in range(interp_in_n.shape[0]):

                matrixA = np.zeros([options.dim+1, options.dim+1])
                for i in range(options.dim+1):
                    matrixA[i] = np.append(unscaled_inputs[:,simp_out[:,j][i]-1], interp_in_n[outputdim][simp_out[:,j][i]-1])

                coords = matrixA

                G = coords.sum(axis=0) / coords.shape[0]

                # run SVD
                u, s, vh = np.linalg.svd(coords - G)
                
                # unitary normal vector
                hyper_sfc_normal = vh[options.dim, :]

                # approx grad as normal scaled by vertical component, times -1
                grad_out = hyper_sfc_normal/hyper_sfc_normal[options.dim]
                grad_out = -grad_out[:-1]
                # print("grad out = ", grad_out)
                grad_est_DS[outputdim][j] = grad_out

            # end loop over output dimns 
    # end if computeGrad
    else:
        grad_est_DS = []


    allow_extrapolation=True
    print_errors=True
    # note: error code 1= sucessful extrap; 2 = extrap beyond threshold
    extrap_indices = np.where((error_out == 1) | (error_out == 2))
    extrap_indices = np.array(extrap_indices[0])
    # print("extrap indices = ", extrap_indices)
    # print("rnorm = ", rnorm_n)
    # print("rnorm[indices] = ", rnorm_n[extrap_indices])
    # print("type = ", type(extrap_indices))
    # print("e i [0]:",extrap_indices[0])
    
    #==============================================================================#
    # Check for errors in DelaunaySparse run
    #==============================================================================#
    if allow_extrapolation: 
        # print("Extrapolation occured at ", np.where(error_out == 1))
        # Next line replaces error code 1 (successful extrapolation) 
        #               with error code 0 (successful interpolation)
        error_out = np.where(error_out == 1, 0, error_out)
    else:
        if 1 in error_out:
            class Extrapolation(Exception): pass
            raise(Extrapolation("Encountered extrapolation point (beyond threshold) when making Delaunay prediction."))
    # Handle any errors that may have occurred.
    if (sum(error_out) != 0):
        if print_errors:
            unique_errors = sorted(np.unique(error_out))
            print(" [Delaunay errors:",end="")
            for e in unique_errors:
                if (e == 0): continue
                indices = tuple(str(i) for i in range(len(error_out))
                                if (error_out[i] == e))
                if (len(indices) > 5): indices = indices[:2] + ('...',) + indices[-2:]
                print(" %3i"%e,"at","{"+",".join(indices)+"}", end=";")
            print("] ")
        # Reset the errors to simplex of 1s (to be 0) and weights of 0s.
        bad_indices = (error_out > (1 if allow_extrapolation else 0))
        simp_out[:,bad_indices] = 1
        weights_out[:,bad_indices] = 0

    return interp_out_n, actual_test_vals, actual_train_vals, extrap_indices, grad_est_DS



#==================================================================================================#
# Main section, includes some bad input checks
#==================================================================================================#
if __name__ == '__main__':

    #==================================================================================================#
    # Provide help screen documentation. Let the user define options. Also define defaults.            #
    #==================================================================================================#

    usage = "%prog [options]"
    parser = OptionParser(usage)
    parser.add_option( "--jobid", help="Job ID.", 
        dest="jobid", type=int, default=999999)   
    parser.add_option( "--staticdatapath", help="Path to static data set from which samples will be drawn. " +
        "If no path is provided, code uses --fn option to sample data.", 
        dest="data_path", type=str, default=None)  
    parser.add_option( "--fn", help="Test function to use.  Version 1.0 of the code only supports the Griewank function.  " +
        "It is possible to code in additional functions by modifying the definition of tf(X).", 
        dest="fn_name", type=str, default="griewank")  
    parser.add_option( "--dim", dest="dim", type=int, default=2, 
        help="Dimension of input space.  Default 2.")
    parser.add_option("--extrap", dest="extrap_thresh", type=float, default=0.0,
        help="Extrapolation threshold parameter passed to DelaunaySparse.  Default 0.0.")
    parser.add_option("--maxsamp", dest="max_samp", type=int, default=20_000,
        help="Max number of samples to draw.  Default = 20,000.")
    parser.add_option("--numtrainpts", dest="numtrainpts", type=int, default=850,
        help="Initial number of samples points (n_0 in the paper).  Default = 850.")
    parser.add_option("--numtestperdim", dest="numtestperdim", type=int, default=20,
        help="Number of test points per dimension. Default = 20.")
    parser.add_option("--logbase", dest="log_base", type=float, default=1.4641,
        help="Upsampling factor b; also the base of the logarithm in rate computation.  Default 1.4641.")
    parser.add_option("--zoomctr", dest="zoom_ctr", type=float, default=0.0,
        help="Zoom modality: used only in conjunction with zoomexp option - see below. " +\
        "Default=0.0.  Use 999.0 in zoomctr or zoomexp to manually specify left/right bounds (not implemented in Version 1.0).")
    parser.add_option("--zoomexp", dest="zoom_exp", type=float, default=1.0,
        help="Zoom modality: set query bounds and bounding box such that (1) center is (x,x,...,x) where x=zoomctr"+\
        " (2) length of query grid is 10e[zoomexp] in each dimension and (3) bounding box determined from testbdsc."+\
        " Default=0.0.  Use 999.0 in zoomctr or zoomexp to manually specify left/right bounds (not implemented in Version 1.0).")
    parser.add_option("--queryleftbd", dest="queryleftbound", type=float, default=0.0,
        help="Left bound of interval used to build query point domain [a, b]^dim. Overwritten if zoom modality is used (see above). Default 0.0")
    parser.add_option("--queryrightbd", dest="queryrightbound", type=float, default=1.0,
        help="Right bound of interval used to build query point domain [a, b]^dim. Overwritten if zoom modality is used (see above). Default 1.0")
    parser.add_option("--bboxleftbd", dest="bboxleftbound", type=float, default=0.0,
        help="Left bound of interval used to build bounding box [a, b]^dim. Overwritten if zoom modality is used (see above). Default 0.0")
    parser.add_option("--bboxrightbd", dest="bboxrightbound", type=float, default=1.0,
        help="Right bound of interval used to build bounding box [a, b]^dim. Overwritten if zoom modality is used (see above). Default 1.0")
    parser.add_option("--testbdsc", dest="tb_scale", type=float, default=0.8,
        help="Query points dimension fraction (qpdf), defined as the side length of the query lattice "
        + "divided by the side length of the bounding box. Default=0.8")
    parser.add_option("--grad", dest="computeGrad", action="store_true", default=True,
	    help="Compute gradients within subroutine that calls DelaunaySparse. Default True.")
    parser.add_option("--outc", dest="out_cor", type=int, default=-1,
        help="Output coordinate to assess.  Default -1 avoids this modality and takes the first output coordinate.")
    parser.add_option("--seed", dest="spec_seed", type=int, default=0,
        help="Value passed as global seed to random number generator.  Default 0.")
    parser.add_option("--itmax", dest="it_max", type=int, default=100,
        help="Max number of iterations.  More robust to use --maxsamp to set threshold.  Default = 100.")
    parser.add_option("--save2static", dest="saveTF", action="store_true", default=False,
	    help="Alternate modality for creating static datasets from test function. If True, sample and evaluate the test function and save to csv file. Default False.")

    
    (options, args) = parser.parse_args()
    

    def echo_options(options):
        print("Selected options:")
        print()
        print("Job ID:      ", options.jobid) 
        print("Global seed for randomization: ", options.spec_seed)
        print()

        if options.data_path == None: # use test function to acquire new data
            print("Function:    ", options.fn_name)
            print("Dimension:   ", options.dim)
            print()
            print("Query points per dim:", options.numtestperdim)
            print("Total number of query points:", options.numtestperdim ** options.dim)

            # Set bounding box left/right bounds based on zoom center, zoom exponent, and scale factor qpdf
            if (options.zoom_ctr != 999 and options.zoom_exp != 999):
                options.bboxleftbound  = np.round(options.zoom_ctr - (10 ** (options.zoom_exp))/options.tb_scale,2)
                options.bboxrightbound = np.round(options.zoom_ctr + (10 ** (options.zoom_exp))/options.tb_scale,2)

            # Set query lattice left/right bounds based on bounding box bounds and scale factor qpdf
            tg_scale_fac = (1.0-options.tb_scale)/2
            interval_width = options.bboxrightbound - options.bboxleftbound
            options.queryleftbound  = options.bboxleftbound  + tg_scale_fac * interval_width
            options.queryrightbound = options.bboxrightbound - tg_scale_fac * interval_width

            print("Query point bounds in each dim: ", "[", options.queryleftbound, ", ", options.queryrightbound, "]")
            print("Query points dimension fraction (qpdf): ", options.tb_scale)
            print("Bounding box bounds in each dim: ", "[", options.bboxleftbound, ", ", options.bboxrightbound, "]")
            print()
            print("Initial sample size:", options.numtrainpts)
            print("Maximum sample size:", options.max_samp)
            print("Upsampling factor b: ", options.log_base)
            print()
            print("Using gradients? : ", options.computeGrad)
            print("Extrapolation threshold: ", options.extrap_thresh)
            # print("Output cor : ", options.out_cor)
            if (options.fn_name != 'griewank'):
                print("==> ERROR: Requested function ", options.fn_name)
                print("Only the function 'griewank' is supported by this version of the code.")
                exit()
        
        else: # static data path provided
            print("Path to static data: ", options.data_path)
            options.fn_name = 'static'
            options.zoom_ctr = 999.0
            options.zoom_exp = 999.0

        print()
        
        if (options.bboxrightbound <= options.bboxleftbound):
            print("Right bound must be larger than left bound")
            exit()
        if (options.tb_scale < 0.0001):
            print("Test bound scale must be > 0")
            exit()
        if options.log_base <= 1:
            print("Log base must be > 1.  Default is 2.0.")
            exit()
        if (options.numtestperdim ** options.dim > 10000):
            print()
            print("==> WARNING: number of query points = (query pts per dim)^(dim) =",options.numtestperdim ** options.dim,"is very large.")
            print("Exiting.")
            exit()
        if (options.extrap_thresh < 0 or options.extrap_thresh > 0.5):
            print()
            print("==> Set extrapolation threshold in [0,0.5]")
            exit()

    echo_options(options)

    globalseed = options.spec_seed
    rng = np.random.default_rng(globalseed)  


    if options.saveTF:
        # ALTERNATE MODALITY: 
        #
        # Randomly sample numtrainpts points with specified parameters
        # Evaluate the given test function on each point
        # Save into a csv file
        # Exit
        #

        data_train_inputs, data_train_outputs = make_random_training_in_box(rng)
        all_data_train = pd.concat([data_train_inputs, data_train_outputs],axis=1)
        outfname = 'temp_generated_data.csv'
        all_data_train.to_csv(outfname, index=False, header=False)    
        print(all_data_train)
        print("Saved samples from above dataframe to file", outfname)
        print("Exiting.")
        exit()

    # torch.manual_seed(globalseed)

    if options.data_path == None: # use test function to acquire new data
        data_train_inputs, data_train_outputs = make_random_training_in_box(rng)    
        data_test_inputs, data_test_outputs = make_test_data_grid(rng)
    else:
        try:
            full_dataset = pd.read_csv(options.data_path, header=None, index_col=None)
            dfrowct = full_dataset.shape[0]
            dfcolct = full_dataset.shape[1]
            print("Read in data from path.  Interpreted as",dfrowct,"points in R^",dfcolct-1,"with one output value per point.\n")
            options.dim = dfcolct-1
            options.max_samp = dfrowct
            print("Setting dimension =",dfcolct-1,"\n")
            print("Initial sample size:", options.numtrainpts)
            print("Query points per dim:", options.numtestperdim)
            print("Total number of query points:", options.numtestperdim ** options.dim)
            print("Upsampling factor b: ", options.log_base)
        except:
            print("\n Error reading in data.  To debug, check that:",
                    "\n  (1) path printed above is correct and",
                    "\n  (2) file is .csv of numerical data where each row is",
                    "the input coordinates followed by 1 output."
                    "\n  (3) sample files from ddd/staticdata/examples/ load sucessfully")
            exit()
        ########################
        # shuffle dataset and make input/output datasets
        ########################
        shuffle_seed = rng.integers(low=0, high=1_000_000, size=1)
        shuffle_seed = shuffle_seed[0] # converts array of size 1 to an integer
        full_dataset = full_dataset.sample(frac=1, random_state=shuffle_seed).reset_index(drop=True)
        print("==> Shuffled dataset based on provided seed.")
        full_dataset_inputs = full_dataset.iloc[:,0:options.dim]  # inputs  = all but last column
        full_dataset_outputs = full_dataset.iloc[:,-1:] # outputs = last column

        ########################
        # create first batch of sample points 
        ########################
        data_train_inputs  = full_dataset_inputs.iloc[0:options.numtrainpts, :]  
        data_train_outputs = full_dataset_outputs.iloc[0:options.numtrainpts, :]  
        
        ########################
        # create grid of query points by looking at percentiles of full data set
        ########################

        cols = full_dataset_inputs.columns
        num_cols = len(cols)
        
        x = np.zeros([num_cols, options.numtestperdim])

        print("==> Description of full dataset:")
        print(full_dataset_inputs.describe())
        
        print("\n==> Setting grid of query points from 25th - 75th percentile")
        i = 0
        for col in cols: 
            # set coordinate mins and maxes according to 25-75 percentiles
            cor_min = full_dataset_inputs.quantile(q=0.25)[col]
            cor_max = full_dataset_inputs.quantile(q=0.75)[col]
            x[i] = np.linspace(cor_min, cor_max, options.numtestperdim)
            i += 1

        mg_in = []
        for i in range(num_cols):
            mg_in.append(x[i])
        grid_pts = np.array(np.meshgrid(*mg_in))
        grid_pts = grid_pts.reshape(options.dim, options.numtestperdim ** options.dim)
        grid_pts = grid_pts.T

        outputs_on_grid = []
        data_test_inputs  = pd.DataFrame(grid_pts)
        data_test_outputs = pd.DataFrame(outputs_on_grid)  # intentionally empty df
        print("==> Query point grid has",data_test_inputs.shape[0],"points in R^"+str(data_test_inputs.shape[1]))

    # end else; now data_{train,test}_{inputs,outputs} are set in either case
    
    if options.data_path == None: # use test function to acquire new data
        outfname = 'zz-' + str(options.jobid) + "-" + str(options.fn_name) + "-d" + str(options.dim) + "-tpd" + str(options.numtestperdim) + "-lb" + str(options.bboxleftbound) + "-rb" + str(options.bboxrightbound) + "-tb" + str(options.tb_scale) + "-log" + str(options.log_base) +".csv"
        if (options.zoom_ctr != 999.0 and options.zoom_exp != 999.0): # add in -zoom[exponent value] before csv
            outfname = outfname[:-4] + "-zoom" + str(options.zoom_exp) + ".csv"
        if (options.spec_seed != 0): # add in -seed[seed value] before csv
            outfname = outfname[:-4] + "-seed" + str(options.spec_seed) + ".csv"
    else:
        outfname = 'zz-' + str(options.jobid) + "-" + str(options.fn_name) + "-d" + str(options.dim) + "-tpd" + str(options.numtestperdim) + "-lb25pct" + "-rb75pct" + "-log" + str(options.log_base) +".csv"
        if (options.spec_seed != 0): # add in -seed[seed value] before csv
            outfname = outfname[:-4] + "-seed" + str(options.spec_seed) + ".csv"
    print("==> Output will be stored in file ",outfname)

    results_df = []
    all_pts_in  = copy.deepcopy(data_train_inputs)
    all_pts_out = copy.deepcopy(data_train_outputs)

    if (options.out_cor == -1): # default
        out_coord = 0 # this means we will only measure error in 0th component of output; no problem if codomain is R^1
    else:
        out_coord = options.out_cor

    print("")
    print("=================================")
    print("For output coordinate ", out_coord,": ")
    print("=== results for ", options.fn_name, " ===")
    print("samples | density | prop extrap | MSD diff | MSD rate | grad diff | grad rate | analytic diff | analytic rate ") 
    print("")

    prev_error = 999999
    prev_vals_at_test = []
    prev_diff = 999999

    ########################################################################
    # create list of number of samples for each update step; 
    #    have to do in advance to avoid rounding issues
    #    can also help in future applications to know sampling rate calculation a priori
    #######################################################################

    quitloop = False
    num_samples_to_add   = np.zeros(options.it_max+1)
    total_samples_so_far = np.zeros(options.it_max+1)
    total_samples_so_far[0] = all_pts_in.shape[0]
    
    for i in range(options.it_max): # i = number of "refinements" of interpolant
        if quitloop:
            break
        #
        # upsampling rule:
        # update number of samples by replacing 2 points per unit per dimension with (logbase + 1) points per unit per dimension      
        #   ALSO: round to the nearst integer and cast as an integer - this is essential for the static data case
        #         otherwise you may add the same sample point more than once, causing an error for DS
        

        total_samples_so_far[i+1] = int(np.round(np.power((options.log_base*np.power(total_samples_so_far[i], 1/options.dim) - (options.log_base - 1)),options.dim)))
        num_samples_to_add[i]     = int(total_samples_so_far[i+1] - total_samples_so_far[i])
        if (total_samples_so_far[i+1] > options.max_samp):
            quitloop = True

    ########################################################################
    # do iterative improvement according to upsampling schedule saved in num_samples_to_add
    #######################################################################

    quitloop = False

    for i in range(options.it_max): # i = number of "refinements" of interpolant

        if quitloop:
            break

        # use the accumulated sample points as the training and the fixed test data sets as the test data
        interp_out_n, actual_test_vals, actual_train_vals, extrap_indices, grad_est_DS = compute_DS_only(all_pts_in, all_pts_out, data_test_inputs, data_test_outputs)

        prop_extrap_iterate = len(extrap_indices)/interp_out_n.shape[1]
        # print('====> proportion extrapolated = %1.2f' % prop_extrap_iterate)
        density_of_sample = all_pts_in.shape[0] ** (1/options.dim)

        if options.data_path == None: # use test function to acquire new data 
            # for analytical functions, we can compute the "actual" rate of convergence, for reference
            ds_vs_actual_at_test = np.sqrt(((interp_out_n[out_coord,:]-actual_test_vals[out_coord,:]) ** 2).mean())
            if (i == 0): 
                error_rate = 0
            else:
                error_rate =  (np.log(prev_error/ds_vs_actual_at_test))/np.log(options.log_base)
        else: # have static data; actual_test_vals = [] which causes an error below
            actual_test_vals = np.zeros_like(interp_out_n)
            ds_vs_actual_at_test = 0
            error_rate = 0

        # difference and rate computation steps for Algorithms 3.1 and 3.2
        if (i == 0):
            new_vs_prev_at_test = 0
            diff_rate = 0
            prev_diff = 0
            if (options.computeGrad):
                grad_new_vs_prev_at_test = 0
                grad_diff_rate = 0
                grad_prev_diff = 0
        elif (i == 1):
            new_vs_prev_at_test = np.sqrt(((interp_out_n[out_coord,:]-prev_vals_at_test[out_coord,:]) ** 2).mean())
            diff_rate = 0
            prev_diff = new_vs_prev_at_test
            if (options.computeGrad):
                grad_new_vs_prev_at_test = np.linalg.norm(grad_est_DS - grad_prev_vals_at_test)
                grad_diff_rate = 0
                grad_prev_diff = grad_new_vs_prev_at_test
        else: # i > 1
            new_vs_prev_at_test = np.sqrt(((interp_out_n[out_coord,:]-prev_vals_at_test[out_coord,:]) ** 2).mean())
            # computation of r_k for MSD rate
            diff_rate = np.log(prev_diff/new_vs_prev_at_test)/np.log(options.log_base) 
            prev_diff = new_vs_prev_at_test
            if (options.computeGrad):
                grad_new_vs_prev_at_test = np.linalg.norm(grad_est_DS - grad_prev_vals_at_test)
                # computation of r_k for grad-MSD rate
                grad_diff_rate = -np.log(grad_new_vs_prev_at_test/grad_prev_diff)/np.log(options.log_base)
                grad_prev_diff = grad_new_vs_prev_at_test

        if (options.computeGrad == False): 
            grad_new_vs_prev_at_test = 0.0
            grad_diff_rate = 0.0


        if (i == 0):
            print(("% 6i & %3.2f & %1.2f & -- & -- & -- & -- & %5.5f & -- \\\\") % (all_pts_in.shape[0], density_of_sample,  prop_extrap_iterate, ds_vs_actual_at_test), flush=True)
        elif (i == 1):
            print(("% 6i & %3.2f & %1.2f & %5.5f & -- & %5.5f & -- & %5.5f & %5.2f \\\\") % (all_pts_in.shape[0], density_of_sample,  prop_extrap_iterate, new_vs_prev_at_test, grad_new_vs_prev_at_test, ds_vs_actual_at_test, error_rate), flush=True)
        else:
            print(("% 6i & %3.2f & %1.2f & %5.5f & %5.2f & %5.5f & %5.2f & %5.5f & %5.2f \\\\") % (all_pts_in.shape[0], density_of_sample,  prop_extrap_iterate, new_vs_prev_at_test, diff_rate, grad_new_vs_prev_at_test, grad_diff_rate, ds_vs_actual_at_test, error_rate), flush=True)
        
        # note all_pts_in.shape[0] should be identical to total_samples_so_far[i]
        results_df.append({
            "dim of intput"    : options.dim,
            "function name"    : options.fn_name,
            "num test points"  : options.numtestperdim ** options.dim,
            "left bound"       : options.bboxleftbound,
            "right bound"      : options.bboxrightbound,
            "test grid scale"  : options.tb_scale,
            "iterate"          : i, 
            "samples"          : all_pts_in.shape[0], 
            "density"          : density_of_sample,  
            "prop extrap"      : prop_extrap_iterate, 
            "iterate diff"     : new_vs_prev_at_test,
            "iterate rate"     : diff_rate, 
            "grad diff"        : grad_new_vs_prev_at_test, 
            "grad rate"        : grad_diff_rate, 
            "actual diff"      : ds_vs_actual_at_test, 
            "actual rate"      : error_rate,
            "log base"         : options.log_base,
            "seed"             : options.spec_seed,
        })

        prev_error = ds_vs_actual_at_test
        prev_vals_at_test = interp_out_n
        if (options.computeGrad):
            grad_prev_vals_at_test = grad_est_DS


            

        ##########################################
        # get more samples for next iteration
        ##########################################
        
        # # now int(round(...)) is done at creation of num_samples_to_add
        # options.numtrainpts = int(np.round(num_samples_to_add[i]))
        options.numtrainpts = int(num_samples_to_add[i])

        # check if we will go over max samples
        if (total_samples_so_far[i] + options.numtrainpts > options.max_samp):
            print("")
            print("====> Next step would go over ", options.max_samp, " samples; breaking loop.")
            # print("====>    ALSO: setting numtrainpts to 1 to avoid generating unused points ")
            options.numtrainpts = 1
            quitloop = True

        if options.data_path == None: # use test function to acquire new data 
            new_sample_inputs, new_sample_outputs = make_random_training_in_box(rng)
        else: # have static data set
            start = int(total_samples_so_far[i])                        # index to start next sample batch
            stop  = int(total_samples_so_far[i]) + options.numtrainpts  # index to end   next sample batch
            new_sample_inputs  = full_dataset.iloc[start:stop, 0:options.dim]
            new_sample_outputs = full_dataset.iloc[start:stop,-1:] # NOTE: check selection of output coord


        ##########################################
        # update sample set for next iteration
        ##########################################

        all_pts_in = pd.concat([all_pts_in, new_sample_inputs])
        all_pts_out = pd.concat([all_pts_out, new_sample_outputs])
        
        columns = [ "dim of intput"    ,
                    "function name"    ,
                    "num test points"  ,
                    "left bound"       ,
                    "right bound"      ,
                    "test grid scale"  ,
                    "iterate"          ,
                    "samples"          ,
                    "density"          ,
                    "prop extrap"      ,
                    "iterate diff"     ,
                    "iterate rate"     ,
                    "grad diff"        ,
                    "grad rate"        ,
                    "actual diff"      ,
                    "actual rate"      ,
                    "log base",
                    "seed",
                ]
        
        df = pd.DataFrame(results_df, columns=columns)
        df.to_csv(outfname)
        
    
    # print("====> final all_pts_out shape was ", all_pts_out.shape)
    # print("===> done with output coord ", out_coord)
    print("===> results saved in ", outfname)

