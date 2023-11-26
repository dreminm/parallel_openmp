import subprocess
import os

parameters = {
    'L': [1.0, 3.0],
    'OMP_NUM_THREADS': [1, 2, 4, 8, 16, 32],
    'grid_size': [128, 256]
}

output_file = 'output.txt'

with open(output_file, 'w') as file:
    for L in parameters['L']:
        for grid in parameters['grid_size']:
            for threads in parameters['OMP_NUM_THREADS']:
                command = ["./main", str(L), str(grid)]

                modified_env = os.environ.copy()
                modified_env["OMP_NUM_THREADS"] = str(threads)

                result = subprocess.run(command, env=modified_env, text=True, capture_output=True)
                
                file.write(f"Output for L={L}, OMP_NUM_THREADS={threads}, grid_size={grid}:\n")
                file.write(result.stdout)
                file.write('\n')

print(f"Outputs saved to {output_file}")