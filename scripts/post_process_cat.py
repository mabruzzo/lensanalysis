import argparse
import logging

from lensanalysis.config_parsing.cosmo_config import \
    CosmoCollectionConfigBuilder
from lensanalysis.config_parsing.procedure_config import ProcedureConfig
from lensanalysis.misc.analysis_collection import AnalysisProductCollection, \
    ConvergenceMapProductCollection, ShearMapProductCollection, \
    FeatureProductCollection, set_analysis_col_value, get_analysis_col_value
from lensanalysis.misc.enum_definitions import DescriptorEnum
from lensanalysis.misc.name_parser import SafeNameParser
from lensanalysis.procedure.procedure import SimpleProcedure


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
parser.add_argument("-f", dest = "force",
                    action = "store_true", default = False,
                    help="force past initial Warning Messages.")
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
parser.add_argument("id",nargs=1,
                    help = "The name of cosmology to post-process.")
parser.add_argument("--min",dest = "min_realization", action = "store",
                    default = None,
                    help = ("Specify the minimum realization to process."
                            "Default is the maximum allowed realization to "
                            "be processed."))
parser.add_argument("--max",dest = "max_realization", action = "store",
                    default = None,
                    help = ("Specify the maximum (inclusive) realization to "
                            "process. Default is the maximum allowed "
                            "realization to be processed."))


def get_generator(min_realization,max_realizaiton,use_mpi = False):
    """
    max_realization is maximum inclusive.
    """
    pass

def setup_saved_products_helper(in_progress, analysis_storage_collection,
                                iterable):
    for elem in cmd_args.save:
        descriptors,object_name = parser.parse_name(elem)
        if DescriptorEnum.tomo in descriptors:
            raise NotImplementedError("Not currently equipped to handle "
                                      "tomography")
        storage = get_analysis_col_value(analysis_storage_collection,
                                         descriptors, object_name)
        if storage is None:
            raise ValueError(("{:s}, {:s} does not have any storage object "
                              "initialized").format(str(descriptors),
                                                    object_name))
        set_analysis_col_value(out, descriptors, object_name)

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

def setup_procedure(cmd_args):

    logging.info("Setting up the procedure")
    # first let's build the cosmo collection config builder
    builder = CosmoCollectionConfigBuilder.from_config_file(cmd_args.config)

    if cmd_args.fiducial:
        cosmo_storage_col = builder.get_fiducial_storage_collection()
    else:
        cosmo_storage_col = builder.get_sampled_storage_collection()

    # if using MPI, we need to ensure that all tasks wait for the first task to
    # load the analysis storage collection for the first time (just in case the
    # storage collection needs to instantiate collections that did not already
    # exist)

    name = cmd_args.id
    if name in cosmo_storage_col:
        analysis_storage = cosmo_storage_col.add_analysis_product_storage(name)
    else:
        analysis_storage = cosmo_storage_col.get_analysis_product_storage(name)

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
                                           analysis_storage_collection)

    # now we do a lttle error checking

    # check to see if we are going backwards in our procedure
    # check to if we are computing unnoisy and noisy maps
    # circumvent these things if force flag was enabled


    # finally, we move to actually building our procedure.
    # start at the end with the features we compute and work our way backwards
    
    # we start with the features and work our way backwards
    

def driver(cmd_args):
    # create logger
    logger = logging.getLogger("log")

    conh = logging.StreamHandler()
    #Verbosity level
    if cmd_args.verbose:
	logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    formatter_console = logging.Formatter('%(levelname)s - %(message)s')
    conh.setFormatter(formatter_console)
    logger.addHandler(conh)

    procedure = setup_procedure(cmd_args)
    


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
    print cmd_args.save

