import logging



conh = logging.StreamHandler()
conh.setFormatter(formatter_console)

formatter_console = logging.Formatter('%(levelname)s - %(message)s')

logger = logging.getLogger("lensanalysis.driver")
logprocedure = logging.getLogger("lensanalysis.procedure")

logger.addHandler(conh)
