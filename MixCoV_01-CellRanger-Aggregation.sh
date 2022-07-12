#!/bin/bash

# aggregate data with Cell Ranger
cellranger aggr --id=PanCoV_scRNA_Agg --csv=PanCoV_CellRanger_AggregationFile.csv --localcores 20  --mempercore 5
