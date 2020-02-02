# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : pickle
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""pickle."""

__author__ = "Nathaniel Starkman"

#############################################################################
# IMPORTS

# GENERAL
import pickle


############################################################################
# CODE
############################################################################


def dump(obj, fname, protocol=None, *, fopt="b", fix_imports=True):
    """Wrap pickle.dump.

    *fname* replaces *file* and is a string for the filename
    this file is auto opened and closed

    """
    with open(fname, "w" + fopt) as file:
        pickle.dump(obj, file, protocol=protocol, fix_imports=fix_imports)

    return


# /def


# ----------------------------------------------------------------------------


def load(
    fname, *, fopt="b", fix_imports=True, encoding="ASCII", errors="strict"
):
    """Wrap pickle.load.

    *fname* replaces *file* and is a string for the filename
    this file is auto opened and closed

    """
    with open(fname, "r" + fopt) as file:
        res = pickle.load(
            file, fix_imports=fix_imports, encoding=encoding, errors=errors
        )

    return res


# /def


############################################################################
# END
