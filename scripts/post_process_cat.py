from lensanalysis.procedure import SimpleProcedure

"""
this script will post process mock shear catalogs.

Ideally, there will ultimately be a script that can post process any lensing 
product.
"""

def get_generator(min_realization,max_realizaiton,use_mpi = False):
    """
    max_realization is maximum inclusive.
    """
    pass

def setup_procedure():
    pass

def driver(procedure):
    pass


"""
I guess there is a sequence of events
- 
- I guess I will allow for starting from an arbitrary point of procedure
- Configuration file will allow for specific sequence
- If there are conflicts about saving we will raise an error. (IE configuration 
  says yes about saving but elsewhere we prevented saving).
"""

if __name__ == '__main__':
    pass
