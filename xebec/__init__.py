import logging
import os


__version__ = "0.4.2"
COOKIE_DIR = os.path.dirname(__file__)


def get_logger(logfile, rulename):
    xebec_logger = logging.getLogger("xebec")
    xebec_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile, mode="w")
    formatter = logging.Formatter(
        f"[%(asctime)s - %(name)s - {rulename}] :: %(message)s"
    )
    fh.setFormatter(formatter)
    xebec_logger.addHandler(fh)

    logging.captureWarnings(True)
    logging.getLogger("py.warnings").addHandler(fh)

    return xebec_logger
