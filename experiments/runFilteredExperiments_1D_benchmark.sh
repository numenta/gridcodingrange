#!/bin/bash

set -e -x

cd "$(dirname "$0")"
pwd

DATA_DIR=data/1D_benchmark_filtered/$(python -c "import time; import uuid; print('{}-{}'.format(time.strftime('%Y%m%d-%H%M%S'), uuid.uuid1()))")
mkdir -p ./$DATA_DIR

python -u generate_bases.py $DATA_DIR --phaseResolution 0.4 0.2 0.1 0.05 0.025 --m 1 2 3 --k 1 --filtered --imposeScales --hardImpose --numTrials 2000
python -u measure_unique_sidelength.py --scaleMinimalBox $DATA_DIR
