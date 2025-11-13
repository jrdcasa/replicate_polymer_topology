#!/bin/bash

#conda deactivate
source ~/Programacion/sandboxes/sandbox_common/bin/activate

# TOPOLOGY
#topology_cmd -i C1_methyl.xsd   -r C1_methyl_headtail.dat   -a C1_methyl_residues.dat   -p C1_methyl
topology_cmd -i PE_model.xsd  -r PE_model_headtail.dat  -a PE_model_residues.dat  -p PE_model
topology_cmd -i mPP_model.xsd -r mPP_model_headtail.dat -a mPP_model_residues.dat -p mPP_model
topology_cmd -i rPP_model.xsd -r rPP_model_headtail.dat -a rPP_model_residues.dat -p rPP_model
