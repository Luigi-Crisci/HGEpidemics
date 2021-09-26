#!/bin/bash
julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/em/sa_params_em.json &> logs/em.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/isolation/sa_params_isolation.json &> logs/isolation.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/lockdown/sa_params_lockdown.json &> logs/lockdown.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/nsga/sa_params.json &> logs/nsga.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/ppm/sa_params_ppm.json &> logs/ppm.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/quarantine/sa_params_quarantine.json &> logs/quarantine.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/tracing/sa_params_tracing.json &> logs/tracing.log &
echo "$!" >> logs/ids.txt