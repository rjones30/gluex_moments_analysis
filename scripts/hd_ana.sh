#!/bin/bash

batch=0
nthreads=4

JANA_GEOMETRY_URL=ccdb://GEOMETRY/main_HDDS.xml
JANA_CALIB_URL=sqlite:////group/halld/www/halldweb/html/dist/ccdb.sqlite
JANA_CALIB_CONTEXT=variation=mc

hd_root --config=hd_ana.config \
        --nthreads=$nthreads \
        -PJANA:BATCH_MODE=$batch \
        -PNTHREADS=$nthreads \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        -PTRK:SAVE_TRUNCATED_DEDX=1 \
       $*
