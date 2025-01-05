#!/usr/bin/env python
# Created on: Jan 3, 2025 at 7:30:53 PM
__author__ = "Michael Cuoco"

# create logger for package
import logging

# Root Logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# Console Handler
formatter = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(name)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

__all__ = ["logger"]
