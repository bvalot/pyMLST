import configparser
import os
import shutil

_ROOT = os.path.abspath(os.path.dirname(__file__))
_CONF_PATH = os.path.join(_ROOT, 'pymlst.conf')

_BIN_SECTION = 'BINARIES'
_LOG_SECTION = 'LOGGING'
_LOG_LEVEL = 'Log Level'

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

def write_config(conf):
    with open(_CONF_PATH, 'w') as file:
        conf.write(file)

def get_config():
    """Retrieves the configuration file."""
    conf = configparser.ConfigParser()
    if os.path.exists(_CONF_PATH):
        conf.read(_CONF_PATH)
    return conf


def update_binary_paths(paths):
    """Updates the paths stored in the configuration file."""
    conf = get_config()

    if not conf.has_section(_BIN_SECTION):
        conf.add_section(_BIN_SECTION)

    for key, value in paths.items():
        conf[_BIN_SECTION][key] = os.path.abspath(value)
    write_config(conf)


def reset_binary_paths():
    """Removes the configuration file."""
    conf = get_config()
    if conf.has_section(_BIN_SECTION):
        conf.remove_section(_BIN_SECTION)
    write_config(conf)


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

def get_logging_level():
    """Return log level"""
    conf = get_config()
    if conf.has_section(_LOG_SECTION):
        return conf.get(_LOG_SECTION, _LOG_LEVEL)
    return("INFO")
    
def set_logging_level(levelname):
    """Defined level of logging"""
    conf = get_config()

    if not conf.has_section(_LOG_SECTION):
        conf.add_section(_LOG_SECTION)
    conf[_LOG_SECTION][_LOG_LEVEL] = levelname  
    write_config(conf)
