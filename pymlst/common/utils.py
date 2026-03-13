import logging
import os
import gzip
from pathlib import Path
from io import TextIOWrapper

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from alembic import command
from alembic.config import Config
from sqlalchemy import create_engine, inspect, select

from pymlst import config
from pymlst.common import flag, exceptions



def records_to_dict(records):
    """Convert SeqIO records to dictionary."""
    seq_dict = {}
    for seq in records:
        seq_dict[seq.id] = seq.seq.upper()
    return seq_dict


# ==================== FASTA FORMAT DETECTION ====================

def is_fasta_handle(handle):
    """
    Check if a file handle contains FASTA format data.
    
    Handles:
    - Standard FASTA headers (starting with '>')
    - Comment lines (starting with ';' or '#')
    - Empty lines
    - Works with non-seekable handles
    
    Args:
        handle: A file handle
    
    Returns:
        bool: True if the handle appears to be in FASTA format
    """
    try:
        # Save current position if seekable
        can_seek = hasattr(handle, 'seekable') and handle.seekable()
        if can_seek:
            current_pos = handle.tell()
        
        # Read lines until we find a FASTA header or run out of lines
        found_header = False
        line_count = 0
        max_lines = 20  # Check up to 20 lines to handle long comment blocks
        
        while line_count < max_lines:
            try:
                line = handle.readline()
                if not line:  # EOF
                    break
                
                line = line.strip()
                line_count += 1
                
                # Skip empty lines
                if not line:
                    continue
                
                # Skip comment lines (common in FASTA files)
                if line.startswith(';') or line.startswith('#'):
                    continue
                
                # Found a FASTA header
                if line.startswith('>'):
                    found_header = True
                    break
                
                # If we find any non-empty, non-comment line that's not '>',
                # and we haven't found a header yet, it's probably not FASTA
                if line_count > 5:  # After 5 non-comment lines without '>', give up
                    break
                
            except Exception:
                break
        
        # Restore position if seekable
        if can_seek:
            try:
                handle.seek(current_pos)
            except Exception:
                pass
        
        return found_header
        
    except Exception:
        # If anything fails, assume it's FASTA (backward compatible)
        return True


# ==================== ENHANCED READ_GENOME ====================

def read_genome(handle):
    """
    Read genome from an open file handle with format validation.
    
    This function:
    1. Detects if the file is gzip compressed (magic bytes)
    2. Creates appropriate streaming handle
    3. Validates FASTA format
    4. Parses sequences
    
    Args:
        handle: An open file handle (text or binary mode)
    
    Returns:
        dict: Dictionary with sequence IDs as keys and sequences as values
    
    Raises:
        ValueError: If the file is not in FASTA format
        IOError: If there are issues reading the file
    """
    original_handle = handle
    parse_handle = None
    file_path = None
    
    # Get file path if available
    if hasattr(handle, 'name'):
        file_path = handle.name
    
    try:
        # ===== STEP 1: Detect file type and create appropriate handle =====
        
        # Check if we can read magic bytes
        can_read_magic = False
        magic = None
        
        try:
            if hasattr(handle, 'seekable') and handle.seekable():
                pos = handle.tell()
                magic = handle.read(2)
                handle.seek(pos)
                can_read_magic = True
        except:
            pass
        
        # Decision tree for handle creation
        if can_read_magic and magic == b'\x1f\x8b':
            # Gzip detected by magic bytes
            logging.debug("Gzip file detected by magic bytes")
            if file_path:
                parse_handle = gzip.open(file_path, 'rt')
            else:
                parse_handle = gzip.open(handle, 'rt')
        
        elif file_path and file_path.endswith('.gz'):
            # Gzip detected by extension
            logging.debug("Gzip file detected by .gz extension")
            parse_handle = gzip.open(file_path, 'rt')
        
        elif hasattr(handle, 'mode') and 'b' in handle.mode:
            # Binary mode but not gzip
            logging.debug("Binary mode handle (non-gzip), converting to text")
            parse_handle = TextIOWrapper(handle, encoding='utf-8')
        
        else:
            # Text mode handle
            logging.debug("Using handle directly")
            parse_handle = handle
        
        # ===== STEP 2: Validate FASTA format =====
        
        if not is_fasta_handle(parse_handle):
            raise ValueError(
                "The input file does not appear to be in FASTA format. "
                "FASTA files should contain headers starting with '>'. "
                "Comments starting with ';' or '#' are ignored."
            )
        
        # ===== STEP 3: Parse sequences =====
        
        # Reset position to beginning
        if hasattr(parse_handle, 'seekable') and parse_handle.seekable():
            parse_handle.seek(0)
        
        # Parse FASTA records (streaming)
        logging.debug("Starting FASTA parsing")
        records = SeqIO.parse(parse_handle, 'fasta')
        result = records_to_dict(records)
        logging.debug(f"FASTA parsing complete: {len(result)} sequences")
        
        return result
        
    except ValueError:
        # Re-raise ValueError as is (for FASTA format errors)
        raise
    except Exception as e:
        logging.error(f"Error in read_genome: {e}")
        # Fallback to original behavior
        try:
            if hasattr(original_handle, 'seekable') and original_handle.seekable():
                original_handle.seek(0)
            records = SeqIO.parse(original_handle, 'fasta')
            return records_to_dict(records)
        except Exception as e2:
            raise IOError(f"Error reading genome: {e2}")
    
    finally:
        # Clean up the parse handle only if we created it
        if parse_handle is not None and parse_handle != original_handle:
            try:
                parse_handle.close()
                logging.debug("Closed parse handle")
            except:
                pass




def write_genome(genome_dict, handle):
    """Write a genome dictionary to a file handle in FASTA format."""
    for seq_id, seq in genome_dict.items():
        handle.write('>' + str(seq_id) + '\n' + str(seq) + '\n')


def file_name(handle):
    """Extract file name without extension from a file handle."""
    filename = os.path.basename(handle.name)
    if filename.endswith(".fasta"):
        return filename.rstrip(".fasta")
    if filename.endswith(".fna"):
        return filename.rstrip(".fna")
    else:
        return filename.split('.')[0]
    

def strip_file(file):
    """Read lines from a file and strip newline characters."""
    found = []  
    if file is not None:
        for line in file:
            found.append(line.rstrip('\n'))
    return found


def compar_seqs(seqs):
    """Compare sequences and count polymorphic sites."""
    count = 0
    for index in range(0, len(seqs[0])):
        seqs_char = {s[index] for s in seqs}
        if '-' in seqs_char:
            seqs_char.remove('-')
        if len(seqs_char) > 1:
            count += 1
    return count


def write_count(count, texte):
    """Write text to a count file if it exists."""
    if count:
        count.write(texte)


def validate_sequence(sequence):
    """Validate that a sequence is a valid coding sequence."""
    try:
        sequence.translate(cds=True, table=11)
    except TranslationError:
        return False
    else:
        return True


def create_logger():
    """Create and configure logger based on configuration."""
    log = config.get_logging_level()
    if log == "DEBUG":
        level = logging.DEBUG
    elif log == "INFO":
        level = logging.INFO
    elif log == "WARNING":
        level = logging.WARNING
    else:
        level = logging.ERROR
    logging.basicConfig(
        level=level,
        format='[%(levelname)s: %(asctime)s] %(message)s')


def clean_kwargs(kwargs):
    """Removes kwargs with None values produced by Click."""
    for key, value in kwargs.copy().items():
        if value is None:
            kwargs.pop(key)
    return kwargs


def get_output(kwargs):
    """Extract output from kwargs for extractor."""
    if 'output' in kwargs:
        out_kwargs = {'output': kwargs['output']}
        kwargs.pop('output')
    else:
        out_kwargs = {}
    return kwargs, out_kwargs


def check_type(conn, mlst_type):
    """Check if database matches expected MLST type."""
    inspector = inspect(conn)
    tables = inspector.get_table_names()
    if len(tables) == 0:
        return
    elif 'mlst_type' not in tables:
        logging.warning('The base missing mlst_type metadata, continue with %s', mlst_type)
        return
    m_t = conn.execute(
        select(flag.mlst_type.c.name)
    ).fetchone()
    if m_t.name != mlst_type:
        raise exceptions.WrongBaseType(
            'The base you are attempting to perform '
            'on belongs to the wrong MLST type')


def get_updated_engine(path, module):
    """Get a SQLAlchemy engine with database updated to latest schema."""
    env_path = config.get_data(os.path.join('alembic', module))
    alembic_cfg = Config()
    alembic_cfg.set_main_option('script_location', env_path)
    logging.getLogger('alembic').setLevel(logging.CRITICAL)

    if path is None:
        engine = create_engine('sqlite://')
    else:
        engine = create_engine('sqlite:///' + os.path.abspath(path))

    with engine.begin() as conn:
        check_type(conn, module)
        alembic_cfg.attributes['connection'] = conn
        command.upgrade(alembic_cfg, 'head')

    return engine


def clean_geneid(geneid):
    """Remove '_' on geneid to be compatible with kma search."""
    return geneid.replace("_", "-")
