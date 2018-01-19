import logging

formatter_console = logging.Formatter('%(levelname)s - %(message)s')

conh = logging.StreamHandler()
conh.setFormatter(formatter_console)

logger = logging.getLogger("lensanalysis.driver")
logprocedure = logging.getLogger("lensanalysis.procedure")

logger.addHandler(conh)
