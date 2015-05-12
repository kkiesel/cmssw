#!/bin/bash

# $1 should be 'fast' or 'full'

# write all output in this file:
exec > >(tee $1FileList.py)

echo "filelist = []"
for dataset in $(das_client.py --limit=0 --query="dataset=/Single*/*kiesel*$1*/* instance=prod/phys03"); do
    for file in $(das_client.py --limit=0 --query="file dataset=$dataset instance=prod/phys03"); do
        echo "filelist.append('"$file"')"
    done
done

