# Load two YAML files and deep-merge them
import yaml
import os

def deep_merge(dest, src):
    for k, v in src.items():
        if k in dest and isinstance(dest[k], dict) and isinstance(v, dict):
            deep_merge(dest[k], v)
        else:
            dest[k] = v

def load_many(files):
    out = {}
    for p in files:
        with open(p, 'r') as fh:
            data = yaml.safe_load(fh) or {}
            deep_merge(out, data)
    return out

files = ['config/2022_v12_samples.yml', 'config/2023_v12_samples.yml']
merged = load_many(files)
with open('config/merged_samples.yml', 'w') as fh:
    yaml.safe_dump(merged, fh)
print("Merged YAML files into:", os.path.abspath('config/merged_samples.yml'))