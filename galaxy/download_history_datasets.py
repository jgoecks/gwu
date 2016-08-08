#
# Download datasets from some Galaxy histories.
#
# TODO: generalize this to download a history's datasets from any history.
#

import os
import subprocess
from bioblend.galaxy import GalaxyInstance

API_KEY = ""
GALAXY_URL = "https://usegalaxy.org"

# Get Galaxy instance.
gi = GalaxyInstance(url=GALAXY_URL, key=API_KEY)

# Download datasets from histories.
for history in gi.histories.get_histories():
    if history['name'].startswith("Plate"):
        print "Downloading datasets from %s" % history['name']
        os.mkdir(history['name'])
        history_metadata = gi.histories.show_history(history['id'], contents=False)
        for dataset_id in history_metadata['state_ids']['ok']:
            dataset_metadata = gi.datasets.show_dataset(dataset_id)
            # print dataset_metadata['name'], dataset_metadata['file_size']
            print "Downloading dataset %s" % dataset_metadata['name']
            subprocess.call(["wget", "-O", os.path.join(history["name"], dataset_metadata['name']),
                            GALAXY_URL + dataset_metadata["download_url"]])
