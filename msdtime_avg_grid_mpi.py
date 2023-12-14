import MDAnalysis as mda
from numba import jit
import os
import numpy as np
import logging
import argparse
from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

log_filename = f"msd_log_rank_{rank}.txt"
logging.basicConfig(filename=log_filename, level=logging.INFO, format='%(asctime)s - %(message)s')

parser = argparse.ArgumentParser(description='Compute Mean Square Displacement (MSD) from MD simulation data.')
parser.add_argument('-t', '--tpr', required=True, help='Path to the input .tpr file.')
parser.add_argument('-r', '--trr', required=True, help='Path to the input .trr file.')
parser.add_argument('-s', '--select', default="all", help='Atom selection string. Default is "all".')
args = parser.parse_args()

trr_filename = os.path.basename(args.trr).split('.')[0]

u = mda.Universe(args.tpr, args.trr)
lipids = u.select_atoms(args.select)

box = u.dimensions[:2]
grid_size = 10
grid_spacing = box / grid_size

@jit(nopython=True)
def compute_euclidean_norm(diff_array):
    return np.sqrt(np.sum(diff_array**2, axis=1))

@jit(nopython=True)
def compute_msd_for_dt(initial_positions, final_positions):
    displacement = compute_euclidean_norm(final_positions - initial_positions)
    msd = np.mean(displacement**2)
    return msd

msds_sets = {i: [] for i in range(grid_size * grid_size)}
total_grids = grid_size * grid_size

# Divide grids among available cores
grids_per_core = total_grids // size
start_grid = rank * grids_per_core
end_grid = (rank + 1) * grids_per_core

# For the last rank, handle potential leftovers
if rank == size - 1:
    end_grid = total_grids

for grid in range(start_grid, end_grid):
    i = grid // grid_size
    j = grid % grid_size

    logging.info(f"Processing grid {grid + 1} out of {total_grids} by Rank {rank}...")
    
    lower_bound = np.array([i * grid_spacing[0], j * grid_spacing[1], -np.inf])
    upper_bound = np.array([(i + 1) * grid_spacing[0], (j + 1) * grid_spacing[1], np.inf])

    lipids_in_current_grid = lipids.select_atoms(
        "prop x >= {} and prop x < {} and prop y >= {} and prop y < {}".format(
            lower_bound[0], upper_bound[0], lower_bound[1], upper_bound[1]
        )
    )

    nframes = len(u.trajectory)
    for dt in range(1, nframes//2 + 1):
        msds_for_dt = []

        for start_frame in range(0, nframes - dt):
            u.trajectory[start_frame]
            initial_positions = lipids_in_current_grid.positions[:, :2].copy()

            u.trajectory[start_frame + dt]
            final_positions = lipids_in_current_grid.positions[:, :2].copy()

            msd = compute_msd_for_dt(initial_positions, final_positions)
            msds_for_dt.append(msd)

        msds_sets[grid].append(np.mean(msds_for_dt))

logging.info(f"Rank {rank}: MSD calculations complete!")

all_msds_sets = comm.gather(msds_sets, root=0)

if rank == 0:
    combined_msds_sets = {}
    for msd_set in all_msds_sets:
        combined_msds_sets.update(msd_set)

    directory = "msd_20fs_grid_time_avg"
    if not os.path.exists(directory):
        os.makedirs(directory)

    for key, msd_values in combined_msds_sets.items():
        output_file_name = f"{trr_filename}_grid_{key}_MSD.txt"
        file_path = os.path.join(directory, output_file_name)
        taus = np.arange(len(msd_values))
        data_to_save = np.column_stack((taus, msd_values))
        np.savetxt(file_path, data_to_save, delimiter='\t', header="Tau\tMSD", comments="")
        logging.info(f"Results for grid {key} saved to {file_path}")
