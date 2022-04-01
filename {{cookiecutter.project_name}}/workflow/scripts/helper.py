import logging


def get_logger(logfile, rulename):
    xebec_logger = logging.getLogger("xebec")
    xebec_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile, mode="w")
    formatter = logging.Formatter(
        f"[%(asctime)s - %(name)s - {rulename}] :: %(message)s"
    )
    fh.setFormatter(formatter)
    xebec_logger.addHandler(fh)

    return xebec_logger
