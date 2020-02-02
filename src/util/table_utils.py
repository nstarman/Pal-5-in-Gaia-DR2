# -*- coding: utf-8 -*-

# ----------------------------------------------------------------------------
#
# TITLE   : Table utilities
# AUTHOR  : Nathaniel Starkman
# PROJECT : AST1501
#
# ----------------------------------------------------------------------------

# Docstring and Metadata
"""Table utilities.

Routine Listings
----------------
neg_to_nan
add_units_to_Table
add_color_col
add_calculated_col
rename_columns
drop_colnames

"""

__author__ = "Nathaniel Starkman"
__credits__ = ["Astropy"]


#############################################################################
# IMPORTS

# GENERAL
import itertools
import numpy as np

from astropy.table import Table, QTable


##############################################################################
# CODE
##############################################################################


def neg_to_nan(df, col):
    """Set negative values in `col` to NaN.

    This edits in-place

    Parameters
    ----------
    df: array_like
        array
    col: index, index str, slicer
        the column selector

    """
    df[col][df[col] < 0] = np.NaN


# /def

# ----------------------------------------------------------------------------


def add_units_to_Table(df, udict=None):
    """Add units to an astropy Table.

    Takes Table and returns QTable

    Parameters
    ----------
    df : Table or QTable
    udict: dict, optional
        dictionary of column names and corresponding units

    Returns
    -------
    qdf: QTable
        same table as `df`, with units

    """
    # Adding Units, if corresponding column in Table
    for key, unit in udict.items():
        if key in df.columns and df[key].unit != unit:
            setattr(df[key], "unit", unit)

    qdf = QTable(df)

    return qdf


# /def


# ----------------------------------------------------------------------------


def add_color_col(df, c1, c2, **kw):
    """Add color column.

    Parameters
    ----------
    df: Table, QTable
    c1: str
    c2: str

    Other Parameters
    ----------------
    color: str
        name of color

    Notes
    -----
    name of color column is `{c1}-{c2}` if not provided as `color`
    adds error column

    TODO
    ----
    make this a function of add_calculated_col

    """
    color = kw.get("color", c1 + "-" + c2)
    if color in df.colnames:
        print("{} already in table".format(color))
        return None

    try:
        c1ind = df.index_column(c1 + "_err")
    except KeyError:
        c1ind = df.index_column(c1)
        noerr = True
    except ValueError:
        c1ind = df.index_column(c1)
        noerr = True
    else:
        noerr = False

    try:
        c2ind = df.index_column(c2 + "_err")
    except KeyError:
        c2ind = df.index_column(c1)
        noerr = True
    except ValueError:
        c2ind = df.index_column(c2)
        noerr = True
    else:
        noerr = False if noerr is False else True

    colindex = max(c1ind, c2ind) + 1
    df.add_column(df[c1] - df[c2], index=colindex, name=color)
    df[color].info.description = color + " color"

    # Adding Errors
    if noerr is False:
        colindex = df.index_column(color) + 1
        df.add_column(
            np.sqrt(df[c1 + "_err"] ** 2 + df[c2 + "_err"] ** 2),
            index=colindex,
            name=color + "_err",
        )
        df[color + "_err"].info.description = "error in {} color [mag]"

    return

# /def


# ----------------------------------------------------------------------------


def add_calculated_col(df, func, *cols, funcargs=[], funckw={}, **kw):
    """Add a calculated column in two column variables.

    Parameters
    ----------
    func: function
        function over df cols
        form function(*columns, *extra_arguments, **keyword_argument)
    cols: list of str
        the names of the columns in df
    funcargs: list
        list of extra arguments for func
    funckw: dict
        dictionary of extra keyword arguments for func

    Other Parameters
    ----------------
    name: str
        function name
        defaults to 'func({},{})'
    index: int
        index to put column at
    description: str
        column description
    return: bool
        whether to return modified df

    Returns
    -------
    if kw['return'] is True:
        returns modified df
    else:
        it modifies in place

    """
    name = kw.get("name", "func{}".format(str(tuple(cols))))

    if name in df.colnames:
        print("{} already in table".format(name))
        return None

    colindex = kw.get("index", None)

    df.add_column(
        func(*[df[c] for c in cols], *funcargs, **funckw),
        index=colindex,
        name=name,
    )

    df[name].info.description = kw.get("description", name)

    if kw.get("return", False) is True:
        return df


# /def


# ----------------------------------------------------------------------------


def rename_columns(df, *args, **kw):
    """Rename columns.

    Parameters
    ----------
    df: Table, QTable
    args: tuple
        [(name, rename), (name, rename), ...]
    kw:
        {name: rename, name: rename, name: rename, ...}
        chains together args, kw
        kw takes precedence

    """
    for n, rn in itertools.chain(args, kw.items()):
        df.rename_column(n, rn)

    return

# /def


# ----------------------------------------------------------------------------


# def drop_colnames(colnames, *args):
#     """helper function for making a table from another table, dropping some names

#     Parameters
#     ----------
#     colnames: list
#         list of names in the original table
#     args: list
#         list of strings of names to drop

#     Returns
#     -------

#     """
#     names = np.array(colnames[:])  # shallow copy just in case
#     inds = ~np.in1d(names, args)

#     return list(names[inds])


# # /def


##############################################################################
# END
