import ast
import re

def molecule_data_parser(filepath):
    keys = [ "verbose", "basis", "spin","charge", "symmetry", "atom"]
    data= {}
    
    with open (filepath, 'r') as f:
        lines = f.read()
        for key in keys:
            if key in lines:
                try:
                    pattern = rf"{key}\s*=\s*(.*?)(?=\n\w+\s*=|\Z)"
                    match = re.search(pattern, lines, re.DOTALL)
                    if match:
                        value = match.group(1).strip()
                        data[key] = ast.literal_eval(value)
                       

                except Exception as e:
                    raise ValueError(
                        f"FIle format error - check {key}, unable to parse value {value}\n"
                        f"{e}"
                        )
    return (
            data.get("verbose"),
            data.get("basis"),
            data.get("spin"),
            data.get("charge"),
            data.get("symmetry"),
            data.get("atom")
)



