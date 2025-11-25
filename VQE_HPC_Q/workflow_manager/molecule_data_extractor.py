import ast
import re

def normalize_atom_list(atom_data):
    """
    Convert atom entries into PySCF-compatible strings.
    Accepts:
    - ["H 0 0 1"]
    - [["Li", (0,0,0)], ["H", (0,0,1.56)]]
    """
    normalized = []

    for item in atom_data:
        # Case 1: already a string
        if isinstance(item, str):
            normalized.append(item.strip())
            continue

        # Case 2: list/tuple like ["Li", (0,0,0)]
        if isinstance(item, (list, tuple)) and len(item) == 2:
            symbol = item[0]
            coords = item[1]

            if isinstance(coords, (list, tuple)) and len(coords) == 3:
                x, y, z = coords
                normalized.append(f"{symbol} {x} {y} {z}")
                continue

        raise ValueError(f"Invalid atom format: {item}")

    return normalized


def molecule_data_parser(filepath):
    """
    Parses molecule file and returns:
    { verbose, basis, spin, charge, symmetry, atom }
    With ATOM automatically normalized into strings.
    """
    keys = ["verbose", "basis", "spin", "charge", "symmetry", "atom"]
    data = {k: None for k in keys}

    with open(filepath, 'r') as f:
        content = f.read()

    for key in keys:
        pattern = rf"^\s*{key}\s*=\s*(.*?)(?=^\s*\w+\s*=|\Z)"
        match = re.search(pattern, content, flags=re.DOTALL | re.MULTILINE)

        if not match:
            continue

        raw_value = match.group(1).strip()
        if raw_value.endswith(","):
            raw_value = raw_value[:-1].strip()

        try:
            value = ast.literal_eval(raw_value)
        except Exception:
            value = raw_value.strip()

        
        if key == "atom" and value is not None:
            value = normalize_atom_list(value)

        data[key] = value

    return data


