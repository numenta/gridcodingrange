#!/bin/bash

set -e -x

cd "$(dirname "$0")"
pwd

DATA_DIR=data/smooth_k_filtered/$(python -c "import time; import uuid; print('{}-{}'.format(time.strftime('%Y%m%d-%H%M%S'), uuid.uuid1()))")
mkdir -p ./$DATA_DIR

python -u generate_bases.py $DATA_DIR --phaseResolution 0.2 --m 5 --k 2.0 2.0001 2.00014142136 2.0002 2.00028284271 2.0004 2.00056568542 2.0008 2.00113137085 2.0016 2.0022627417 2.0032 2.0045254834 2.0064 2.0090509668 2.0128 2.0181019336 2.0256 2.0362038672 2.0512 2.07240773439 2.1024 2.14481546879 2.2048 2.28963093757 2.4096 2.57926187515 2.8192 3.0 3.0001 3.00014142136 3.0002 3.00028284271 3.0004 3.00056568542 3.0008 3.00113137085 3.0016 3.0022627417 3.0032 3.0045254834 3.0064 3.0090509668 3.0128 3.0181019336 3.0256 3.0362038672 3.0512 3.07240773439 3.1024 3.14481546879 3.2048 3.28963093757 3.4096 3.57926187515 3.8192 4.0 --normalizeScales --filtered --reuseBases --numTrials 100
python -u measure_unique_sidelength.py $DATA_DIR