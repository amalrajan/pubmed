import logging
from typing import Optional


def get_logger(name: Optional[str] = None, debug: bool = False) -> logging.Logger:
    """Set up and return a configured logger.

    :param name: Name of the logger, defaults to None
    :type name: Optional[str]
    :param debug: If True, set logging level to DEBUG
    :type debug: bool
    :return: Configured logger
    :rtype: logging.Logger
    """
    logger_name = name or "pubmed_tool"
    logger = logging.getLogger(logger_name)

    # Remove all existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Set logging level
    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)

    # Create console handler and set level
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # Create formatter and attach to handler
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)

    # Add new handler to logger
    logger.addHandler(ch)

    return logger
