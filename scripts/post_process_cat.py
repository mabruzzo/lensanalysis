import argparse
import logging
import sys
import traceback

from mpi4py import MPI

from lensanalysis.config_parsing.cosmo_config import \
    CosmoCollectionConfigBuilder
from lensanalysis.config_parsing.procedure_config import ProcedureConfig
from lensanalysis.misc.analysis_collection import AnalysisProductCollection, \
    ConvergenceMapProductCollection, ShearMapProductCollection, \
    FeatureProductCollection, set_analysis_col_value, get_analysis_col_value
from lensanalysis.misc.enum_definitions import Descriptor
from lensanalysis.misc.log import logger
from lensanalysis.misc.name_parser import SafeNameParser

from lensanalysis.procedure.building import build_procedure
from lensanalysis.procedure.io_step import LoadCollection
from lensanalysis.procedure.procedure import SimpleProcedure, DataPacket

"""
this script will post process mock shear catalogs.

Ideally, there will ultimately be a script that can post process any lensing 
product.
"""

description = ("A program that post processes Shear Catalogs and computes "
               "features to constrain cosmologies.\n\n"
               "A feature of the program is to dynamically construct the "
               "procedure that will be followed. Over the course of any "
               "procedure various analysis objects are created. To achieve the "
               "dynamic feature building, the analysis object that is created "
               "from which everything else is computed, must be specified. "
               "Then the procedure is built so that the specified analysis "
               "objects can be constructed and saved.\n\n"
               "Presently, the only features that can be computed are related "
               "to Peak Counting. \n\n"
               "The names of feature objects understood by the program are:"
               "... (fill this in later)")
               
parser = argparse.ArgumentParser(description = description)
parser.add_argument("-v", "--verbose", dest = "verbose",
                    action = "store_true", default = False,
                    help="turn output verbosity")
#parser.add_argument("-f", dest = "force",
#                    action = "store_true", default = False,
#                    help="force past initial Warning Messages.")
parser.add_argument("-b", "--begin", dest = "begin", action="store",
                    help = "Specify where to start in the procedure")
parser.add_argument("-p", "--procedure", dest = "procedure_config",
                    action = "store", required = True,
                    help = "Specify procedure configuration file")
parser.add_argument("-c", "--config", dest = "config", required = True,
                    help = ("Specify the configuration file that details how "
                            "all files are (or should be) saved."))
parser.add_argument("-s", "--save", dest = "save", action="store",
                    nargs = '+', default = None,
                    help = ("Specify a list of analysis objects to save. "
                            "This argument can NOT be called alongside -i "
                            "or --include."))
parser.add_argument("-i", "--include", dest = "include", action = "store",
                    nargs = '+', default = None,
                    help = ("Specify a list of analysis objects to save in "
                            "addition to the analysis objects specified in "
                            "the procedure configuration file. This "
                            "argument can NOT be called alongside -s or "
                            "--save."))
parser.add_argument("--fiducial",dest = "fiducial", action = "store_true",
                    default = False,
                    help = ("Specify that the analysis is being performed on "
                            "the fiducial cosmology (rather than a sampled "
                            "cosmology."))
parser.add_argument("--min",dest = "min_realization", action = "store",
                    default = None, type = int,
                    help = ("Specify the minimum realization to process."
                            "Default is the maximum allowed realization to "
                            "be processed."))
parser.add_argument("--max",dest = "max_realization", action = "store",
                    default = None, type =int,
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed."))
parser.add_argument("id",nargs=1,
                    help = "The name of cosmology to post-process.")



def setup_saved_products_helper(in_progress, analysis_storage_collection,
                                iterable):
    parser = SafeNameParser()

    for elem in iterable:
        descriptors,object_name = parser.parse_name(elem)
        if Descriptor.tomo in descriptors:
            raise NotImplementedError("Not currently equipped to handle "
                                      "tomography")
        storage = get_analysis_col_value(analysis_storage_collection,
                                         descriptors, object_name)
        if storage is None:
            raise ValueError(("{:s}, {:s} does not have any storage object "
                              "initialized").format(str(descriptors),
                                                    object_name))
        set_analysis_col_value(in_progress, descriptors, 
                               object_name, True)

def determine_saved_products(cmd_args, proc_config,
                             analysis_storage_collection):
    """
    Comes up with an instance of AnalysisProductCollection with True values 
    for every analysis product we want to save.
    """

    out = AnalysisProductCollection()
    out.conv_map = ConvergenceMapProductCollection()
    out.shear_map = ShearMapProductCollection()
    out.feature_products = FeatureProductCollection()

    parser = SafeNameParser()

    if cmd_args.save is not None:
        setup_saved_products_helper(out, analysis_storage_collection,
                                    cmd_args.save)
    else:
        if cmd_args.include is not None:
            setup_saved_products_helper(out, analysis_storage_collection,
                                        cmd_args.include)
        setup_saved_products_helper(out, analysis_storage_collection,
                                    proc_config.get_default_saved_products())
    return out


def _starting_procedure_step(cmd_args,name,analysis_storage,cosmo_storage_col):
    parser = SafeNameParser()

    if cmd_args.begin is None:
        begin = ((),"shear_cat")
    else:
        begin = parser.parse_name(cmd_args.begin)
    if begin[1] == 'shear_cat':
        assert len(begin[0]) == 0
        loader = cosmo_storage_col.get_shear_cat_loader(name)
    else:
        loader = get_analysis_col_value(analysis_storage, *begin)
        assert loader is not None

    first_step = LoadCollection(loader)
    return begin, first_step

def simple_realization_generator(min_realization,max_realization):
    """
    max_realization is maximum inclusive.
    """
    for i in range(min_realization,max_realization+1):
        logging.info("Beginning Realization {:05d}".format(i))
        yield None, DataPacket(i)

def _setup_generator(cmd_args,proc_config,nprocs=1,rank=0):
    
    if cmd_args.min_realization is not None:
        min_realization = cmd_args.min_realization
        assert isinstance(min_realization, int)
        assert min_realization >= proc_config.get_min_realization()
    else:
        min_realization = proc_config.get_min_realization()
    assert min_realization > 0
        
    if cmd_args.max_realization is not None:
        max_realization = cmd_args.max_realization
        assert isinstance(max_realization,int)
        assert max_realization <= proc_config.get_max_realization()
    else:
        max_realization = proc_config.get_max_realization()

    assert max_realization>= min_realization

    if nprocs >1:
        num_real = max_realization - min_realization +1
        num_per_task = num_real // nprocs
        cur_real = num_per_task
        if (num_real%nprocs)>0:
            if rank < (num_real%nprocs):
                cur_real +=1
                min_delta = rank*(num_per_task+1)
            else:
                min_delta = (num_real%nprocs)*(num_per_task+1)
                min_delta += (rank-(num_real%nprocs))*num_per_task
        else:
            min_delta = rank*(num_per_task)
        min_r = min_realization+min_delta
        max_r = min_r + cur_real-1

        if num_real<rank+1:
            return None
    else:
        min_r = min_realization
        max_r = max_realization

    return simple_realization_generator(min_r,max_r)

def _setup_mpi_helper(builder,comm):
    """
    This is a helper function that get's the storage collection. 

    It has separated from the rest of the setup function because the master 
    process must call this process while the other processes wait and then 
    the other processes must call this. By doing this we can avoid issues 
    where different processes try to write a directory at the same time.
    """
    if cmd_args.fiducial:
        cosmo_storage_col = builder.get_fiducial_storage_collection()
    else:
        print "Sampled"
        cosmo_storage_col = builder.get_sampled_storage_collection()
    return cosmo_storage_col

def setup(cmd_args,comm):
    """
    Parameters
    ----------
    cmd_args : dict
        The dictionary of parsed command line arguments
    comm : mpi4py.MPI.comm
        The mpi communications object or None. If the value is None, then MPI 
        is not in use.
    """

    logging.info("Setting up the procedure")
    # first let's build the cosmo collection config builder
    builder = CosmoCollectionConfigBuilder.from_config_file(cmd_args.config)

    if comm is None:
        nprocs = 1
        rank = 0 
    else:
        nprocs = comm.Get_size()
        rank = comm.Get_rank()

    if nprocs == 1:
        cosmo_storage_col = _setup_mpi_helper(builder,comm)
    else:
        rank = comm.Get_rank()
        # check an possible create the directories on rank 0.
        # The other tasks will wait around while this is going on
        if rank == 0:
            cosmo_storage_col = _setup_mpi_helper(builder,comm)
        comm.Barrier()

        # now the other tasks with locate the directories. Rank 0 is free to 
        # move on since there is no possibility of an error occuring due to 
        # creating directories in the same place at the same time.
        if rank != 0:
            cosmo_storage_col = _setup_mpi_helper(builder,comm)

    if cmd_args.fiducial:
        cosmo_storage_col = builder.get_fiducial_storage_collection()
    else:
        print "Sampled"
        cosmo_storage_col = builder.get_sampled_storage_collection()

    # if using MPI, we need to ensure that all tasks wait for the first task to
    # load the analysis storage collection for the first time (just in case the
    # storage collection needs to instantiate collections that did not already
    # exist)

    name = cmd_args.id[0]
    if name in cosmo_storage_col:
        analysis_storage = cosmo_storage_col.get_analysis_product_storage(name)
    else:
        analysis_storage = cosmo_storage_col.add_analysis_product_storage(name)

    # NOW, you can let all of the other mpi tasks catch up.


    # Now, read in the procedure configuration file
    proc_config = ProcedureConfig.from_fname(cmd_args.procedure_config)
    
    if cmd_args.save is not None and cmd_args.include is not None:
        raise ValueError("-s/--save and -i/--include cannot both be set")

    # at this point, create an instance of analysis product collection and fill
    # in the appropriate locations with True/None (whether or not we wish to
    # save the file.
    # while doing this check to ensure that such storage locations exist
    save_config = determine_saved_products(cmd_args, proc_config,
                                           analysis_storage)
    # now we do a lttle error checking

    # check to see if we are going backwards in our procedure
    # check to if we are computing unnoisy and noisy maps
    # circumvent these things if force flag was enabled


    # finally, we move to actually building our procedure.
    
    # we start with the features and work our way backwards
    begin_object, first_step = _starting_procedure_step(cmd_args,name,
                                                        analysis_storage,
                                                        cosmo_storage_col)

    # build the remaining steps
    remaining_steps = build_procedure(begin_object, proc_config,
                                      save_config, analysis_storage)
    first_step.wrapped_step = remaining_steps

    procedure = SimpleProcedure()
    procedure.add(first_step)

    logging.info("Setting up the Generator")
    # the barrier is here for debugging

    generator = _setup_generator(cmd_args, proc_config,
                                 nprocs = nprocs, rank = rank)
    return procedure,generator

def driver(cmd_args):
    # set Verbosity level
    if cmd_args.verbose:
	logging.basicConfig(level = logging.DEBUG)
    else:
        logging.basicConfig(level = logging.INFO)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    if nprocs == 1:
        comm = None
    print "CHECKPOINT A"

    if comm is not None:
        try:
            procedure,generator = setup(cmd_args,comm)
            if generator is None:
                return None
            procedure.apply_to_iterable(generator)
        except:
            print "Unexpected error on rank {:d}:".format(rank)
            print traceback.print_exception(*sys.exc_info())
            MPI.COMM_WORLD.Abort(1)
    else:
        procedure,generator = setup(cmd_args,comm)
        procedure.apply_to_iterable(generator)
    


"""
I guess there is a sequence of events
- 
- I guess I will allow for starting from an arbitrary point of procedure
- Configuration file will allow for specific sequence
- If there are conflicts about saving we will raise an error. (IE configuration 
  says yes about saving but elsewhere we prevented saving).
"""

if __name__ == '__main__':
    #Parse command arguments
    cmd_args = parser.parse_args()
    driver(cmd_args)

