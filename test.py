import os
def load_config(config_file="config.txt"):
    config = {}
    with open(config_file, "r") as file:
        for line in file:
            if line.strip() and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                config[key.strip()] = value.strip()
    return config

# Load configuration
config = load_config()

catalog_dir = config.get("catalog_dir")
ilc_map_dir = config.get("ilc_map_dir")

# Ensure output directory exists
#os.makedirs(catalog_dir, exist_ok=True)

print(f"Reading from: {catalog_dir}")
print(f"Saving results to: {ilc_map_dir}")