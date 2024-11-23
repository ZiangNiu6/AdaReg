# write activate renv
Rscript -e 'renv::activate(); renv::restore()'

# 1. Extract the gRNA names from data
eps_list=(0.1 0.2 0.4)
directory=code/formal-simulation

# 2. null simulation
# Submit the job for the current eps
output_null_dir="$directory/results/null"
mkdir -p "$output_null_dir"
echo "Rscript $directory/simulation-null.R $output_null_dir" | qsub -N "null_simulation" -l m_mem_free=12G


# 3. alternative simulation
# Loop over each epsilon parameter in eps_list
for eps in "${eps_list[@]}"
do
    # Create the output directory for the current epsilon
    output_alternative_dir="$directory/results/alternative_eps_${eps}"
    mkdir -p "$output_alternative_dir"
    
    # jump to next iteration if the file already exists
    if [ -f "$output_alternative_dir/results.rds" ]; then
        echo "Results already exist for eps=$eps. Skipping."
        continue
    fi
    
    # Submit the job for the current eps
    echo "Rscript $directory/simulation-alternative.R $eps $output_alternative_dir" | qsub -N "alternative_simulation_eps_${eps}" -l m_mem_free=12G
    
    echo "Submitted jobs for $eps greedy algorithm!"
done
    
