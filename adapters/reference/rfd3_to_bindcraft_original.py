#!/usr/bin/env python3
"""
Convert RFD3 .cif.gz outputs to bindcraft PDB format.

Steps:
1. Convert CIF to PDB
2. Swap chain labels (B <-> A)
3. Relabel residues in chain A (after swap) to account for chain breaks
   - Each geometric gap in chain gets +1 added to residue numbering
4. Save with format: l{chain_B_length}_s{5_digit_timestamp}.pdb
"""

import gzip
import os
import sys
from pathlib import Path
from datetime import datetime
from Bio import PDB
from Bio.PDB import PPBuilder
import time


def check_chain_breaks(chain):
    """
    Detect geometric chain breaks in a chain.
    Returns list of (before_residue_number, after_residue_number) tuples.
    """
    residues = list(chain)
    if len(residues) < 2:
        return []
    
    breaks = []
    for i in range(len(residues) - 1):
        res1 = residues[i]
        res2 = residues[i + 1]
        
        # Check if both residues have CA atoms
        if "CA" not in res1 or "CA" not in res2:
            continue
        
        # Calculate distance between CA atoms
        ca1 = res1["CA"]
        ca2 = res2["CA"]
        distance = ca1 - ca2
        
        # Chain break if distance > 4 Angstroms (typical peptide bond ~3.8)
        if distance > 4.0:
            breaks.append((res1.get_id()[1], res2.get_id()[1]))
    
    return breaks


def relabel_residues_for_chain_breaks(chain, breaks):
    """
    Relabel residues to account for chain breaks.
    For each break before a residue, increment that residue's number by 1.
    This creates a gap of 1 in residue numbering (residues separated by 2 indices).
    
    Args:
        chain: BioPython chain object
        breaks: List of (before_resnum, after_resnum) tuples indicating chain breaks
    """
    if not breaks:
        return  # No breaks, no relabeling needed
    
    residues = list(chain)
    
    # Create a mapping of old residue numbers to new ones
    renumbering = {}
    for residue in residues:
        current_resnum = residue.get_id()[1]
        
        # Count how many breaks occur BEFORE or AT this residue
        # A break (before_num, after_num) is before/at this residue if after_num <= current_resnum
        offset = sum(1 for before_num, after_num in breaks if after_num <= current_resnum)
        
        renumbering[current_resnum] = current_resnum + offset
    
    # Renumber in descending order of new ID to avoid conflicts
    residues_sorted = sorted(residues, key=lambda r: -renumbering[r.get_id()[1]])
    
    for residue in residues_sorted:
        old_resnum = residue.get_id()[1]
        new_resnum = renumbering[old_resnum]
        residue.id = (residue.id[0], new_resnum, residue.id[2])


def convert_rfd3_to_bindcraft(cif_gz_path, output_dir=None, prefix=""):
    """
    Convert a single RFD3 CIF.GZ file to bindcraft PDB format.
    
    Args:
        cif_gz_path: Path to .cif.gz file
        output_dir: Output directory (default: same as input)
        prefix: Optional prefix to add to filename
    
    Returns:
        Path to output PDB file
    """
    cif_gz_path = Path(cif_gz_path)
    
    if not cif_gz_path.exists():
        raise FileNotFoundError(f"Input file not found: {cif_gz_path}")
    
    if output_dir is None:
        output_dir = cif_gz_path.parent
    else:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse CIF from gzip
    print(f"Processing: {cif_gz_path.name}")
    
    with gzip.open(cif_gz_path, 'rt') as f:
        cif_parser = PDB.MMCIFParser(QUIET=True)
        structure = cif_parser.get_structure(cif_gz_path.stem, f)
    
    # Get the first model
    model = structure[0]
    chains = list(model)
    
    # Assert exactly 2 chains
    if len(chains) != 2:
        raise ValueError(f"Expected 2 chains, found {len(chains)}: {[c.id for c in chains]}")
    
    # Get chain A and B
    chain_a = None
    chain_b = None
    for chain in chains:
        if chain.id == 'A':
            chain_a = chain
        elif chain.id == 'B':
            chain_b = chain
    
    if chain_a is None or chain_b is None:
        raise ValueError(f"Chains A and B not found. Found: {[c.id for c in chains]}")
    
    # Record length of original chain A (will become new chain B after swap)
    chain_b_length = len(chain_a)
    
    # Swap chain IDs: A becomes B, B becomes A
    chain_a.id = 'X'  # Temporary
    chain_b.id = 'A'
    chain_a.id = 'B'
    
    # Now re-get the chain references after swap
    # New chain A is the one with id='A' (was original B)
    # New chain B is the one with id='B' (was original A)
    new_chain_a = None
    new_chain_b = None
    for chain in model:
        if chain.id == 'A':
            new_chain_a = chain
        elif chain.id == 'B':
            new_chain_b = chain
    
    # Relabel residues in new chain A (originally chain B) for chain breaks
    breaks = check_chain_breaks(new_chain_a)
    if breaks:
        print(f"  Chain breaks found: {len(breaks)} breaks")
        relabel_residues_for_chain_breaks(new_chain_a, breaks)
        print(f"  Residues relabeled for breaks")
    else:
        print(f"  No chain breaks detected")
    
    # Generate output filename
    # Generate 5-digit timestamp
    timestamp_ms = int((time.time() * 1000) % 100000)
    timestamp_str = f"{timestamp_ms:05d}"
    
    # Build filename with optional prefix
    if prefix:
        output_filename = f"{prefix}_l{chain_b_length}_s{timestamp_str}.pdb"
    else:
        output_filename = f"_l{chain_b_length}_s{timestamp_str}.pdb"
    output_path = output_dir / output_filename
    
    # Write PDB
    pdb_writer = PDB.PPBuilder()
    pdb_io = PDB.PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(str(output_path))
    
    print(f"  Saved to: {output_filename}")
    
    return output_path


def batch_convert(input_dir, output_dir=None, prefix="", pattern="*.cif.gz"):
    """
    Convert all CIF.GZ files in a directory.
    
    Args:
        input_dir: Directory containing .cif.gz files
        output_dir: Output directory (default: same as input)
        prefix: Optional prefix to add to filenames
        pattern: File pattern to match
    
    Returns:
        List of output file paths
    """
    input_dir = Path(input_dir)
    
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    cif_files = sorted(input_dir.glob(pattern))
    
    if not cif_files:
        print(f"No files matching pattern '{pattern}' found in {input_dir}")
        return []
    
    print(f"Found {len(cif_files)} files to convert\n")
    
    output_paths = []
    errors = []
    
    for cif_file in cif_files:
        try:
            output_path = convert_rfd3_to_bindcraft(cif_file, output_dir, prefix=prefix)
            output_paths.append(output_path)
        except Exception as e:
            error_msg = f"ERROR: {cif_file.name} - {str(e)}"
            print(f"  {error_msg}")
            errors.append(error_msg)
    
    print(f"\n{'='*60}")
    print(f"Conversion complete!")
    print(f"Successfully converted: {len(output_paths)}/{len(cif_files)}")
    
    if errors:
        print(f"Errors ({len(errors)}):")
        for error in errors:
            print(f"  {error}")
    
    return output_paths


if __name__ == "__main__":
    # Example usage: script.py <input_path> [output_dir] [prefix]
    if len(sys.argv) > 1:
        input_path = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else None
        prefix = sys.argv[3] if len(sys.argv) > 3 else ""
    else:
        # Default: use the RFD3_Kinesin_hotspot outputs
        input_path = "/home/tadas/code/startup_lanternfish_binders/binder_runs/RFD3_Kinesin_hotspot/outputs_non_loopy"
        output_dir = None
        prefix = ""
    
    input_path = Path(input_path)
    
    if input_path.is_file():
        # Single file
        convert_rfd3_to_bindcraft(input_path, output_dir, prefix=prefix)
    elif input_path.is_dir():
        # Batch convert directory
        batch_convert(input_path, output_dir, prefix=prefix)
    else:
        print(f"Path not found: {input_path}")
        sys.exit(1)
