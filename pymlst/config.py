import configparser
import os
import shutil

_ROOT = os.path.abspath(os.path.dirname(__file__))
_CONF_PATH = os.path.join(_ROOT, 'pymlst.conf')

_BIN_SECTION = 'BINARIES'


def get_data(path):
    return os.path.join(_ROOT, 'data', path)


def get_config():
    """Retrieves the configuration file."""
    conf = configparser.ConfigParser()
    if os.path.exists(_CONF_PATH):
        conf.read(_CONF_PATH)
        return conf
    return conf


def update_binary_paths(paths):
    """Updates the paths stored in the configuration file."""
    conf = get_config()

    if not conf.has_section(_BIN_SECTION):
        conf.add_section(_BIN_SECTION)

    for key, value in paths.items():
        conf[_BIN_SECTION][key] = os.path.abspath(value)

    with open(_CONF_PATH, 'w') as file:
        conf.write(file)


def reset_binary_paths():
    """Removes the configuration file."""
    conf = get_config()
    if conf.has_section(_BIN_SECTION):
        conf.remove_section(_BIN_SECTION)


def get_binary_path(bin_name):
    """Retrieves a binary path."""
    conf = get_config()
    if conf.has_option(_BIN_SECTION, bin_name):
        return conf.get(_BIN_SECTION, bin_name)
    return shutil.which(bin_name)  # path research


def list_binary_paths():
    """Lists the binary paths stored in the configuration file."""
    conf = get_config()
    if conf.has_section(_BIN_SECTION):
        return conf.items(_BIN_SECTION)
    return []
