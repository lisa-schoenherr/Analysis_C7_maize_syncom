# IMPORTS
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed


# FUNCTION
def fetch_metabolite_details(met):
    bigg_id = met[:-2]
    url = f"{base_url}/{bigg_id}"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()  # converts JSON response to a dictionary (data)
            bigg_model = data.get("compartments_in_models", [])

            bigg_model_id = []
            for model in bigg_model:
                bigg_model_id.append(model["model_bigg_id"])

            return {
                "bigg_id": bigg_id,
                "name": met,
                "bigg_models": list(set(bigg_model_id))
            }

    except Exception as e:
        print(f"Error with {bigg_id}: {e}")

    return {
        "bigg_id": bigg_id,
    }


# MAIN
# Get list of all universal metabolites
base_url = "http://bigg.ucsd.edu/api/v2/universal/metabolites"
response = requests.get(base_url)

# list of metabolites where you wanna get Bigg Model ID from
metabolis = []

# check if request is going through
if response.status_code != 200:
    raise Exception("Failed to fetch metabolite list")

# Use ThreadPoolExecutor to parallelise requests
# Set max_workers to 25 for faster execution
results = []
with ThreadPoolExecutor(max_workers=25) as executor:
    futures = [executor.submit(fetch_metabolite_details, met) for met in metabolis]
    for i, future in enumerate(as_completed(futures)):
        results.append(future.result())

df = pd.DataFrame(results)
df