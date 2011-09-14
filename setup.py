#!/usr/bin/env python
#
# Distutils setup script for munkres_sparse
#
# $Id$
# ---------------------------------------------------------------------------

from distutils.core import setup
import re
import os
import sys
import imp

# Load the data.

here = os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path = [here] + sys.path
mf = os.path.join(here, 'munkres_sparse.py')
munkres_sparse = imp.load_module('munkres_sparse', open(mf), mf,
                          ('__init__.py', 'r', imp.PY_SOURCE))
long_description = munkres_sparse.__doc__
version = str(munkres_sparse.__version__)
author = munkres_sparse.__author__
email = munkres_sparse.__email__

url = munkres_sparse.__url__
license = munkres_sparse.__license__

# Run setup

setup(
    name="munkres_sparse",
    version=version,
    description="kuhn-munkres algorithm for the Assignment Problem with an incomplete cost matrix",
    long_description=long_description,
    url=url,
    license=license,
    author=author,
    author_email=email,
    py_modules=["munkres_sparse"],
    classifiers = [
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GPLv2',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Mathematics', 
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
