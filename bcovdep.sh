#!/bin/bash

TOOLS_DIR=/work/user/yanll/Tools/ZH_pipelines/tools

${TOOLS_DIR}/covdepstat \
        <(cat sample_info2.txt) \
        /work/user/yanll/Tools/ZH_pipelines/data/regions \
        > bcovdep-summary2.txt

