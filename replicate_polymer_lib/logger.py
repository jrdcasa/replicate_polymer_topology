import logging


# =============================================================================
class CustomFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    FORMATS = {
        logging.ERROR: "\n\tERROR: %(asctime)s: %(msg)s",
        logging.WARNING: "\n\tWARNING: %(msg)s",
        logging.DEBUG: "%(asctime)s: %(msg)s",
        "DEFAULT": "%(msg)s",
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        date_fmt = '%d-%m-%Y %d %H:%M:%S'
        formatter = logging.Formatter(log_fmt, date_fmt)
        return formatter.format(record)


# =============================================================================
def init_logger(name, fileoutput="output.log", append=True, inscreen=False):
    """

    Args:
        name:
        fileoutput:
        append:
        inscreen:

    Returns:

    """

    loggerName = logging.getLogger(name=name)
    loggerName.setLevel(logging.DEBUG)

    # create a file handler
    if append:
        h1 = logging.FileHandler(fileoutput, 'a')
    else:
        h1 = logging.FileHandler(fileoutput, 'w')

    # create a logging format
    h1.setFormatter(CustomFormatter())

    # add the handlers to the logger
    loggerName.addHandler(h1)

    # Output also in screen
    if inscreen:
        f1 = logging.StreamHandler()
        f1.setFormatter(CustomFormatter())
        loggerName.addHandler(f1)

    return loggerName
