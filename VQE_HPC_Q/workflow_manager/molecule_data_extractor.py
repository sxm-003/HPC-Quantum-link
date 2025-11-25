import ast
import re

def molecule_data_parser(filepath):
    keys = ["verbose", "basis", "spin", "charge", "symmetry", "atom"]
    data = {}

    with open(filepath, 'r') as f:
        content = f.read()

    for key in keys:
        pattern = rf"^\s*{key}\s*=\s*(.*?)(?=^\s*\w+\s*=|\Z)"
        match = re.search(pattern, content, flags=re.DOTALL | re.MULTILINE)

        if not match:
            # No key found — return None or raise error?
            data[key] = None
            continue

        raw_value = match.group(1).strip()

        try:
            # Try python literal parsing
            value = ast.literal_eval(raw_value)
        except Exception:
            # fallback: treat as a raw string
            value = raw_value

        data[key] = value

    return tuple(data[k] for k in keys)


print(molecule_data_parser("test.txt"))
