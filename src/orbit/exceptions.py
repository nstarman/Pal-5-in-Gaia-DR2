# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : exceptions
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Exceptions."""

__author__ = "Nathaniel Starkman"


#############################################################################
# IMPORTS

# GENERAL
import warnings
from warnings import warn


#############################################################################
# CODE
#############################################################################


class integrationWarning(UserWarning):
    """integrationWarning"""
    pass


warnings.simplefilter("always", integrationWarning)


#############################################################################
# END
