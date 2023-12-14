import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import pdist, squareform
import os
import logging
import argparse
from mpi4py import MPI
import xdrlib

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

log_filename = f"msd_log_rank_{rank}.txt"
logging.basicConfig(filename=log_filename, level=logging.INFO, format='%(asctime)s - %(message)s')

parser = argparse.ArgumentParser(description='Compute correlation from MD simulation data.')
parser.add_argument('-t', '--tpr', required=True, help='Path to the input .tpr file.')
parser.add_argument('-r', '--trr', required=True, help='Path to the input .trr file.')
parser.add_argument('-s', '--select', default="name PO4", help='Atom selection string. Default is "name PO4".')
args = parser.parse_args()

#@jit(nopython=True)
def distance(xi,yi,xt,yt,L):
    dx=xt-xi
    dy=yt-yi
    corrx=np.round(dx/L)
    corry=np.round(dy/L)
    dx=dx-corrx*L
    dy=dy-corry*L
    dr = np.zeros((len(dx),2))
    dr[:,0] = dx
    dr[:,1] = dy
    mag = np.linalg.norm(dr, axis=1)
    dr = dr / mag[:, np.newaxis]
    return dr

#@jit(nopython=True)
def corr_func(i, dist_nd, unit_disp):
    dim = 2
    cr_i = np.zeros(dim)
    for j in range(N):
        if i == j:
            continue
        dist = dist_nd[i,j]
        if r <= dist <= (r + dr):
            cr_i += unit_disp[j]
    return cr_i

u = mda.Universe(args.tpr, args.trr)
lipids = u.select_atoms(args.select)
trr_filename = os.path.basename(args.trr).split('.')[0]

L = 300
dt = 100
dr = 0.2
d = 0.2
R = np.arange(0.0001, L/2, dr)
dim = 2
N = len(lipids)
rho = N/(L**2)
ensembles = np.arange(0, 800, 20)
N_ensembles = len(ensembles)

# Split the ensembles among the MPI processes
local_ensembles = np.array_split(ensembles, size)[rank]

local_gr_en = np.zeros(len(R))

logging.info(f"Rank {rank} processing frames: {local_ensembles[0]} to {local_ensembles[-1]}")
local_ensembles = np.array_split(ensembles, size)[rank]
local_gr_en = np.zeros(len(R))

for ti in local_ensembles:
    u.trajectory[ti+dt]
    final_pos = lipids.positions
    rt = lipids.positions[:, :2]
    xt = lipids.positions[:, 0]
    yt = lipids.positions[:, 1]

    u.trajectory[ti]
    initial_pos = lipids.positions
    ri = lipids.positions[:, :2]
    xi = lipids.positions[:, 0]
    yi = lipids.positions[:, 1]

    dist_nd_sq = np.zeros(N*(N-1)//2)
    for di in range(dim):
        pos_1d = ri[:,di][:,np.newaxis]
        dist_1d = pdist(pos_1d)
        dist_1d[dist_1d > L*0.5] -= L
        dist_nd_sq += dist_1d**2
    dist_nd = np.sqrt(dist_nd_sq)
    dist_nd = squareform(dist_nd)

    unit_disp = distance(xi,yi,xt,yt,L)

    gr = np.array([])
    for r in R:
        cr = 0
        for i in range(N):
            cr_dot = corr_func(i, dist_nd, unit_disp)
            cr_i = np.dot(cr_dot, unit_disp[i])/(2*np.pi*r*dr*rho)
            cr += cr_i
        cr = cr/N
        gr = np.append(gr, cr)

    local_gr_en += gr

# Gather results from all MPI processes
all_gr_en = np.zeros_like(local_gr_en)
comm.Reduce(local_gr_en, all_gr_en, op=MPI.SUM, root=0)

if rank == 0:
    gr_en = all_gr_en / N_ensembles
    
    # Make sure the directory exists
    directory = "correlation_data"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Assuming `key` is already defined somewhere in your code
    output_file_name = f"{trr_filename}_cor.txt"
    file_path = os.path.join(directory, output_file_name)
    
    # Stacking data and saving
    data_to_save = np.column_stack((R, gr_en))
    np.savetxt(file_path, data_to_save, delimiter='\t', header="Radial_Range\tCorrelation", comments="")

    logging.info(f"Results saved to {file_path}")
