import requests
import pandas as pd
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed

# Get list of all universal metabolites
base_url = "http://bigg.ucsd.edu/api/v2/"
list_url = base_url + "universal/metabolites"
response = requests.get(list_url)

# check if request is going through
if response.status_code != 200:
    raise Exception("Failed to fetch metabolite list")

metabolites = response.json()["results"]
print(f"Found {len(metabolites)} metabolites. Fetching details...")

# function that fetches specific information for one metabolite, i.e. BIGG ID, name, formulas and charges
def fetch_metabolite_details(met):
    bigg_id = met.get("bigg_id", "")
    name = met.get("name", "")
    url = f"{base_url}universal/metabolites/{bigg_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()  # converts JSON response to a dictionary (data)
            formulae = data.get("formulae", [])  # if formula not available, use empty list []
            charges = data.get("charges", [])

            # Safe quoting for CSV
            safe_name = str(name).replace('"', "'")
            name = f'"{safe_name}"'

            return {
                "bigg_id": bigg_id,
                "name": name,
                "formulas": str(formulae),
                "charges": str(charges)
            }
    except Exception as e:
        print(f"Error with {bigg_id}: {e}")

    safe_name = str(name).replace('"', "'")
    name = f'"{safe_name}"'
    return {
        "bigg_id": bigg_id,
        "name": name,
        "formulas": "[]",
        "charges": "[]"
    }


# Use ThreadPoolExecutor to parallelise requests (faster runtime)
results = []
with ThreadPoolExecutor(max_workers=25) as executor:
    futures = [executor.submit(fetch_metabolite_details, met) for met in metabolites]
    for i, future in enumerate(as_completed(futures)):
        results.append(future.result())
        if i % 500 == 0:
            print(f"{i}/{len(metabolites)} done...")

# Save to CSV
df = pd.DataFrame(results)
df.to_csv("../bigg_metabolites_complete.csv", index=False, quoting=csv.QUOTE_MINIMAL)
print("Saved to bigg_metabolites_complete.csv")


# Read previously created CSV with all BIGG metabolites
df_bigg_met = pd.read_csv("../bigg_metabolites_complete.csv", quotechar='"')

# Convert stringified lists back to real lists
df_bigg_met["formulas"] = df_bigg_met["formulas"].apply(ast.literal_eval)
df_bigg_met["charges"] = df_bigg_met["charges"].apply(ast.literal_eval)