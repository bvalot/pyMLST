#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
python Mlst Local Search Tool

.. currentmodule:: pymlst
.. moduleauthor:: benoit_valot <benoit.valot@univ-fcomte.fr>
"""

from . import config
from .version import __version__, __release__  # noqa
from .wg.core import open_wg
from .cla.core import open_cla
