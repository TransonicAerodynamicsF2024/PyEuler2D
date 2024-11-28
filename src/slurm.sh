#!/encs/bin/tcsh

#SBATCH --job-name=Euler2DC1   ## Give the job a name
#SBATCH --mem=64Gb              ## Memory allocation
#SBATCH --time=30:00:00
#SBATCH --main-type=BEGIN,END
#SBATCH --main-user=paramvir.lobana@mail.concordia.ca

# Load Python module
module load python/3.9.7

# Display the loaded modules for verification <OPTIONAL>
module list

# Python script parameters
# Replace the values of the flags below with your specific arguments

srun python your_script.py \
    -i 1000 \           # Number of iterations
    -C 0.8 \            # CFL number
    -m /path/to/mesh \  # Path to mesh file
    -A 5 \              # Angle of attack in degrees
    -M 0.9              # Mach number
