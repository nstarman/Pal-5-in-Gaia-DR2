# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : paths initialization file
# AUTHOR  : Nathaniel Starkman
# PROJECT : starkython
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""initialization file for util/paths.

Routine Listings
----------------
current_file_directory

"""

__author__ = "Nathaniel Starkman"

__all__ = [
    "current_file_directory"
]


##############################################################################
# CODE
##############################################################################


def current_file_directory(__file__):
    """current_file_directory.

    Parameters
    ----------
    __file__ : builtin
        the furrent file directory
        pops the file from the directory, leaving the the folder

    Returns
    -------
    cur_dir : str
        the current directory

    """
    cur_dir = "/".join(__file__.split("/")[:-1])

    return cur_dir


# /def


##############################################################################
# END
