"""Binary dependencies management utility functions.

This module provides tools to manage the binary dependencies
within the application. Those binaries are external dependencies
and are left to the user to provide.

The following functions are used within the application to retrieve
binaries location that might have been defined by the user.
"""

import os

import configparser
import shutil

from pkg_resources import Requirement, resource_filename


config_path = resource_filename(Requirement.parse("pymlst"), "pymlst.conf")


def update_binary_paths(paths):
    """Updates the paths stored in the configuration file."""
    config = get_config()
    if config is None:
        config = configparser.ConfigParser()

    if not config.has_section('BINARIES'):
        config.add_section('BINARIES')

    for key, value in paths.items():
        config['BINARIES'][key] = os.path.abspath(value)

    with open(config_path, 'w') as file:
        config.write(file)


def reset_config():
    """Removes the configuration file."""
    if os.path.exists(config_path):
        os.remove(config_path)


def get_config():
    """Retrieves the configuration file."""
    if os.path.exists(config_path):
        config = configparser.ConfigParser()
        config.read(config_path)
        return config
    return None


def get_binary_path(bin_name):
    """Retrieves a binary path."""
    config = get_config()
    if config is not None:
        if config.has_option('BINARIES', bin_name):
            return config.get('BINARIES', bin_name)
    return shutil.which(bin_name)  # path research


def list_binary_paths():
    """Lists the binary paths stored in the configuration file."""
    config = get_config()
    if config is not None:
        if config.has_section('BINARIES'):
            return config.items('BINARIES')
    return []
