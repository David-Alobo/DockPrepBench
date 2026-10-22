"""Logging configuration."""

import logging


def setup_logging(log_level: int = logging.INFO) -> None:
    """
    Set up logging configuration.

    Args:
        log_level (int): Logging level (e.g., logging.DEBUG, logging.INFO).
    """
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s — %(name)s — %(levelname)s — %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.getLogger().setLevel(log_level)
