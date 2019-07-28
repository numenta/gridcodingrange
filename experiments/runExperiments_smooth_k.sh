#!/bin/bash

set -e -x

cd "$(dirname "$0")"
pwd

DATA_DIR=data/smooth_k/$(python -c "import time; import uuid; print('{}-{}'.format(time.strftime('%Y%m%d-%H%M%S'), uuid.uuid1()))")
mkdir -p ./$DATA_DIR

python -u generate_bases.py $DATA_DIR --phaseResolution 0.2 --m 3 4 5 6 --k 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 --numTrials 100
python -u measure_unique_sidelength.py $DATA_DIR
