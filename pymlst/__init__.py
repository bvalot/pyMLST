#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
python Mlst Local Search Tool

.. currentmodule:: pymlst
.. moduleauthor:: benoit_valot <benoit.valot@univ-fcomte.fr>
"""
import configparser

from .version import __version__, __release__  # noqa
from pkg_resources import Requirement, resource_filename
import os
import shutil


config_path = resource_filename(Requirement.parse("pymlst"), "config.conf")


def update_binary_paths(paths):
    config = get_config()
    if config is None:
        config = configparser.ConfigParser()

    if not config.has_section('BINARIES'):
        config.add_section('BINARIES')

    for k, v in paths.items():
        config['BINARIES'][k] = os.path.abspath(v)

    with open(config_path, 'w') as f:
        config.write(f)


def reset_config():
    if os.path.exists(config_path):
        os.remove(config_path)


def get_config():
    if os.path.exists(config_path):
        config = configparser.ConfigParser()
        config.read(config_path)
        return config
    else:
        return None


def get_binary_path(bin_name):
    config = get_config()
    if config is not None:
        if config.has_option('BINARIES', bin_name):
            return config.get('BINARIES', bin_name)
    return shutil.which(bin_name)  # path research


def list_binary_paths():
    config = get_config()
    if config is not None:
        if config.has_section('BINARIES'):
            return config.items('BINARIES')
    return []
