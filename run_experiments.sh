#!/bin/bash
truncate --size 0 logs/ids.txt

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

# ### Immunization ###

# julia src/experiments/immunization/spreading_experiment_with_immunization.jl src/experiments/immunization/SA/configs/age_immunization_uniform/sa_params.json &> logs/age_immunization_uniform.log &
# echo "$!" >> logs/ids.txt

# julia src/experiments/immunization/spreading_experiment_with_immunization.jl src/experiments/immunization/SA/configs/age_immunization_mixed/sa_params.json &> logs/age_immunization_mixed.log &
# echo "$!" >> logs/ids.txt

### Timed NPIs ###

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/mixed/lockdown--ppm-em-tracing--immunize/sa_params.json &> logs/lockdown--ppm-em-tracing.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/mixed/ppm-em-quarantine/sa_params.json &> logs/ppm-em-quarantine.log &
echo "$!" >> logs/ids.txt

julia src/experiments/NPIs/spreading_experiment_with_NPIs.jl src/experiments/NPIs/SA/configs/npis/mixed/lockdown--nothing/sa_params.json &> logs/lockdown--nothing.log &
echo "$!" >> logs/ids.txt