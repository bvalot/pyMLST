"""add CLI command file."""

import os
import subprocess
import tempfile
from contextlib import contextmanager
import click

import pymlst
from pymlst.common import utils, exceptions


@contextmanager
def open_genome_file(genome_path):
    """
    Context manager for handling compressed genome files.
    Creates a temporary file for decompressed content.
    
    Args:
        genome_path: Path to the genome file (.fasta or .fasta.gz)
        
    Yields:
        Path to the temporary file containing genome data
    """
    temp_file = None
    
    try:
        if genome_path.endswith('.gz'):
            # Create temporary file
            temp_file = tempfile.NamedTemporaryFile(mode='w+t', suffix='.fasta', delete=False)
            temp_file.close()
            
            # Decompress using zcat
            with open(temp_file.name, 'w') as f_out:
                result = subprocess.run(
                    ['zcat', genome_path],
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True
                )
                f_out.write(result.stdout)
            
            # Yield the temporary file path
            yield temp_file.name
            
        else:
            # Regular file - yield the original path
            yield genome_path
                
    finally:
        # Note: We don't delete the temp file here because it will be used
        # by BLAT. The caller is responsible for cleaning up.
        pass


@click.command(name='add')
@click.option('--strain', '-s',
              help='Name of the strain (default: genome name).')
@click.option('--identity', '-i',
              type=float, 
              default=0.95,
              help='Minimum identity to search gene (default=0.95).')
@click.option('--coverage', '-c',
              type=float,
              default=0.90,
              help='Minimum coverage to search gene (default=0.9).')
@click.argument('database',
                type=click.Path(exists=True))
@click.argument('genome',
                type=click.Path(exists=True))

def cli(database, genome, strain, identity, coverage):
    """Adds a strain GENOME to the wgMLST DATABASE.
    
    GENOME can be a FASTA file or a compressed FASTA file (.gz)
    """
    temp_file_path = None
    
    try:
        # Prepare arguments
        kwargs = {
            'strain': strain,
            'identity': identity,
            'coverage': coverage
        }
        
        # Handle compressed file if needed
        with open_genome_file(genome) as genome_file_path:
            temp_file_path = genome_file_path if genome_file_path != genome else None
            
            # Open the genome file for reading
            with open(genome_file_path, 'r') as genome_handle:
                with pymlst.open_wg(os.path.abspath(database)) as mlst:
                    mlst.add_strain(genome_handle, **utils.clean_kwargs(kwargs))
        
    except exceptions.PyMLSTError as err:
        raise click.ClickException(str(err))
    except subprocess.CalledProcessError as err:
        raise click.ClickException(f"Error decompressing file: {err}")
    except Exception as err:
        raise click.ClickException(f"Error: {err}")
    finally:
        # Clean up temporary file
        if temp_file_path and os.path.exists(temp_file_path):
            os.unlink(temp_file_path)
