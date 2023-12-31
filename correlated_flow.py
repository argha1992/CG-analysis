import os
import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import pdist, squareform
from numba import jit, prange
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process .tpr and .trr files for correlation analysis.")
    parser.add_argument("tpr_file", help="Input .tpr file path.")
    parser.add_argument("trr_file", help="Input .trr file path.")
    args = parser.parse_args()
    return args.tpr_file, args.trr_file

tpr_file, trr_file = parse_arguments()

#@jit(nopython=True)
def distance(xi, yi, xt, yt, L):
    dx = xt - xi
    dy = yt - yi
    corrx = np.round(dx / L)
    corry = np.round(dy / L)
    dx = dx - corrx * L
    dy = dy - corry * L
    dr = np.zeros((len(dx), 2))
    dr[:, 0] = dx
    dr[:, 1] = dy
    mag = np.linalg.norm(dr, axis=1)
    dr = dr / mag[:, np.newaxis]
    return dr

#@jit(nopython=True)
def corr_func(i, dist_nd, unit_disp, r, dr, rho):
    cr_disp_i = np.zeros(2)
    cr_i = 0
    for j in prange(dist_nd.shape[0]):
        if i == j:
            continue
        dist = dist_nd[i, j]
        if r <= dist <= (r + dr):
            cr_disp_i += unit_disp[j]
            cr_i += 1
    return cr_disp_i, cr_i

u = mda.Universe(tpr_file, trr_file)
lipids = u.select_atoms("name PO4")
L = 300
dt = 100
dr = 0.2
d = 0.2
R = np.arange(0.0001, L/2, dr)
dim = 2
gr_en = np.zeros(len(R))
gr_disp_en = np.zeros(len(R))
ensembles = np.arange(0, 1000, 10)
N_ensembles = len(ensembles)
N = len(lipids)
rho = N / (L**2)

for ti in ensembles:
    u.trajectory[ti+dt]
    final_pos = lipids.positions
    rt = final_pos[:,:2]
    xt = final_pos[:,0]
    yt = final_pos[:,1]

    u.trajectory[ti]
    initial_pos = lipids.positions
    ri = initial_pos[:,:2]
    xi = initial_pos[:,0]
    yi = initial_pos[:,1]

    dist_nd_sq = np.zeros(N*(N-1)//2)
    for di in range(dim):
        pos_1d = ri[:, di][:, np.newaxis]
        dist_1d = pdist(pos_1d)
        dist_1d[dist_1d > L * 0.5] -= L
        dist_nd_sq += dist_1d**2
    dist_nd = np.sqrt(dist_nd_sq)
    dist_nd = squareform(dist_nd)
    unit_disp = distance(xi, yi, xt, yt, L)

    gr = np.array([])
    gr_disp = np.array([])
    for r in R:
        cr_disp = 0
        cr = 0
        for i in range(N):
            cr_dot, cr_i = corr_func(i, dist_nd, unit_disp, r, dr, rho)
            cr_disp_i = np.dot(cr_dot, unit_disp[i]) / (2 * np.pi * r * dr * rho)
            cr_disp += cr_disp_i
            cr += cr_i / (2 * np.pi * r * dr * rho)
        cr = cr / N
        cr_disp = cr_disp / N
        gr = np.append(gr, cr)
        gr_disp = np.append(gr_disp, cr_disp)

    gr_en += gr
    gr_disp_en += gr_disp

gr_en = gr_en / N_ensembles
gr_disp_en = gr_disp_en / N_ensembles

result = np.column_stack((R, gr_en, gr_disp_en))

base_filename = os.path.splitext(os.path.basename(trr_file))[0]
dirname = 'correlation_data_1'

if not os.path.exists(dirname):
    os.makedirs(dirname)

result_file_path = os.path.join(dirname, f'{base_filename}_10_cor.txt')

np.savetxt(result_file_path, result, header='R gr_en, gr_disp_en', comments='', delimiter='\t', fmt='%10.5f')
