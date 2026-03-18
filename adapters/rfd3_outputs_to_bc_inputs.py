#!/usr/bin/env python3
import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resseq: int
    icode: str


def parse_atom_line(line: str) -> Optional[Tuple[ResidueKey, str]]:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    chain = line[21]
    resseq = int(line[22:26])
    icode = line[26]
    return ResidueKey(chain=chain, resseq=resseq, icode=icode), line


def parse_xyz(line: str) -> Tuple[float, float, float]:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def format_atom_line(line: str, chain: str, resseq: int, icode: str) -> str:
    return f"{line[:21]}{chain}{resseq:4d}{icode}{line[27:]}"


def collect_residue_order(lines: List[str]) -> List[ResidueKey]:
    seen = set()
    order: List[ResidueKey] = []
    for line in lines:
        parsed = parse_atom_line(line)
        if parsed is None:
            continue
        key, _ = parsed
        if key not in seen:
            seen.add(key)
            order.append(key)
    return order


def load_target_map(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
    if path is None:
        return {}
    data = json.loads(path.read_text())
    return data.get("residue_map", {})


def collect_ca_coords(lines: List[str]) -> Dict[ResidueKey, Tuple[float, float, float]]:
    coords: Dict[ResidueKey, Tuple[float, float, float]] = {}
    for line in lines:
        parsed = parse_atom_line(line)
        if parsed is None:
            continue
        key, _ = parsed
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        if key not in coords:
            coords[key] = parse_xyz(line)
    return coords


def find_chain_break_offsets(
    residue_order: List[ResidueKey],
    ca_coords: Dict[ResidueKey, Tuple[float, float, float]],
    chain_id: str,
    break_threshold: float = 4.0,
) -> Dict[ResidueKey, int]:
    chain_residues = [r for r in residue_order if r.chain == chain_id]
    offsets: Dict[ResidueKey, int] = {}
    breaks_seen = 0
    prev_key: Optional[ResidueKey] = None
    for key in chain_residues:
        if prev_key is not None and prev_key in ca_coords and key in ca_coords:
            x1, y1, z1 = ca_coords[prev_key]
            x2, y2, z2 = ca_coords[key]
            dist = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
            if dist > break_threshold:
                breaks_seen += 1
        offsets[key] = breaks_seen
        prev_key = key
    return offsets


def build_new_numbering(
    residue_order: List[ResidueKey],
    target_numbering_map: Dict[str, Dict[str, str]],
    target_break_offsets: Dict[ResidueKey, int],
) -> Dict[ResidueKey, Tuple[str, int, str]]:
    mapping: Dict[ResidueKey, Tuple[str, int, str]] = {}

    binder_counter = 1
    for key in residue_order:
        if key.chain == "A":
            mapping[key] = ("B", binder_counter, " ")
            binder_counter += 1
            continue

        if key.chain == "B":
            target_key = f"{key.resseq}{key.icode.strip()}"
            mapped = target_numbering_map.get("B", {}).get(target_key)
            if mapped:
                num = "".join(ch for ch in mapped if ch.isdigit())
                icode = "".join(ch for ch in mapped if not ch.isdigit()) or " "
                resnum = int(num)
                mapping[key] = ("A", resnum + target_break_offsets.get(key, 0), icode[0])
            else:
                mapping[key] = ("A", key.resseq + target_break_offsets.get(key, 0), key.icode)
            continue

        raise ValueError(f"Unexpected chain in RFD3 output: {key.chain!r}")
    return mapping


def rewrite_pdb(
    lines: List[str], mapping: Dict[ResidueKey, Tuple[str, int, str]]
) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
    rewritten: List[str] = []
    residue_map: Dict[str, Dict[str, str]] = {"A_to_B": {}, "B_to_A": {}}

    for line in lines:
        parsed = parse_atom_line(line)
        if parsed is None:
            rewritten.append(line)
            continue
        old_key, atom_line = parsed
        new_chain, new_resseq, new_icode = mapping[old_key]
        old_repr = f"{old_key.resseq}{old_key.icode.strip()}"
        new_repr = f"{new_resseq}{new_icode.strip()}"

        if old_key.chain == "A":
            residue_map["A_to_B"][old_repr] = new_repr
        elif old_key.chain == "B":
            residue_map["B_to_A"][old_repr] = new_repr

        rewritten.append(format_atom_line(atom_line, new_chain, new_resseq, new_icode))
    return rewritten, residue_map


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Adapter for RFD3 outputs to standardized BC inputs."
    )
    parser.add_argument("--input-pdb", required=True, type=Path)
    parser.add_argument("--target-numbering-map", type=Path, default=None)
    parser.add_argument("--output-pdb", required=True, type=Path)
    parser.add_argument("--output-chain-map", required=True, type=Path)
    args = parser.parse_args()

    lines = args.input_pdb.read_text().splitlines(keepends=True)
    residue_order = collect_residue_order(lines)
    target_map = load_target_map(args.target_numbering_map)
    ca_coords = collect_ca_coords(lines)
    target_break_offsets = find_chain_break_offsets(residue_order, ca_coords, chain_id="B")
    mapping = build_new_numbering(residue_order, target_map, target_break_offsets)
    rewritten, residue_map = rewrite_pdb(lines, mapping)

    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    args.output_chain_map.parent.mkdir(parents=True, exist_ok=True)

    args.output_pdb.write_text("".join(rewritten))
    chain_meta = {
        "input_chain_semantics": {"A": "binder", "B": "target"},
        "output_chain_semantics": {"A": "target", "B": "binder"},
        "residue_mapping": residue_map,
        "target_break_offsets_applied": True,
    }
    args.output_chain_map.write_text(json.dumps(chain_meta, indent=2))


if __name__ == "__main__":
    main()
