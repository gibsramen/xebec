import logging


def get_logger(logfile):
    xebec_logger = logging.getLogger("xebec")
    xebec_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile, mode="w")
    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] :: %(message)s"
    )
    fh.setFormatter(formatter)
    xebec_logger.addHandler(fh)

    return xebec_logger
