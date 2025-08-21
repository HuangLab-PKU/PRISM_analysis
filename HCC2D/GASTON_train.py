# basic packages
import warnings
warnings.filterwarnings('ignore')

import os
from pathlib import Path
import shutil

# data processing packages
import numpy as np
import scanpy as sc
from glmpca import glmpca
from gaston import neural_net

# visualization packages
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "xelatex",      # 使用 XeLaTeX，如果不需要 LaTeX 公式渲染，可以省略
    'font.family': 'serif',          # 字体设置为衬线字体
    'text.usetex': False,            # 禁用 LaTeX，使用 Matplotlib 内置文字渲染
    'pgf.rcfonts': False,            # 禁用 pgf 的默认字体管理
    'pdf.fonttype': 42,              # 确保字体为 TrueType 格式，可被 Illustrator 编辑
    'ps.fonttype': 42,               # EPS 文件也使用 TrueType 格式
    'figure.dpi': 300,               # 设置图形分辨率
    'savefig.dpi': 300,              # 保存的图形文件分辨率
    'axes.unicode_minus': False,     # 避免负号问题
})

import torch
if torch.cuda.is_available(): torch.set_default_device('cuda')
else: torch.set_default_device('cpu')
print("Default device:", torch.tensor([1.0]).device)


# workdir
BASE_DIR = Path(r'G:\spatial_data\analysis')
RUN_ID = '20230523_HCC_PRISM_probe_refined'
base_path = BASE_DIR / RUN_ID
typ_path = base_path / "cell_typing"
output_path = base_path / "GASTON" / 'demo_pc5'
output_path.mkdir(exist_ok=True)

# copy current file to output_path, existing files will be overwritten
shutil.copy(__file__, output_path / Path(__file__).name)

## preprocess
# Load the data
combine_adata_st = sc.read_h5ad(typ_path / 'combine_adata_st.h5ad')
adata_direct = sc.read_h5ad(typ_path / 'adata_leiden_res_1.h5ad')
adata = adata_direct[adata_direct.obs.index.isin(combine_adata_st.obs.index)]
adata.obs = combine_adata_st.obs
adata.obs = adata.obs.rename(columns={'X_pos':'X', 'Y_pos':'Y'})
adata.obsm['spatial'] = adata.obs.loc[:, ['X', 'Y']].values
adata.obs['batch'] = adata.obs['dataset']
adata = adata[adata.obs['type'] != 'other']

# select a subset of marker genes
# marker_gene_list = [
#     'HBV', 'AFP', 'GPC3', 
# ]
# adata = adata[:, marker_gene_list]
print(adata)

counts_mat = np.array(adata.raw.X)
coords_mat = adata.obsm['spatial']
cell_types = adata.obs.subtype.values
gene_labels = np.array(adata.var.index)

## run GLM-PCA
# GLM-PCA parameters
num_dims = 5
penalty = 1.5

# CHANGE THESE PARAMETERS TO REDUCE RUNTIME
num_iters = 30
eps = 1e-4
glmpca_res = glmpca.glmpca(counts_mat.T, num_dims, fam="poi", penalty=penalty, verbose=True,
                        ctl = {"maxIter":num_iters, "eps":eps, "optimizeTheta":True})
A = glmpca_res['factors'] # should be of size N x num_dims, where each column is a PC

# visualize top GLM-PCs
R, C = -(-num_dims//5), 5
fig,axs = plt.subplots(R, C, figsize=(C*5, R*5))
for r in range(R):
    for c in range(C):
        i = r*C+c
        ax_tmp = axs[r,c] if R > 1 else axs[c]
        ax_tmp.scatter(-1*coords_mat[:, 0], coords_mat[:, 1], c=A[:, i], cmap='Reds', s=0.5)
        ax_tmp.set_title(f'GLM-PC{i}')
plt.savefig(output_path / 'glmpca.png')
plt.close()

np.save(output_path / 'glmpca.npy', A)
np.save(output_path / 'coords_mat.npy', coords_mat)
np.save(output_path / 'counts_mat.npy', counts_mat)
np.save(output_path / 'cell_labels.npy', cell_types)
np.save(output_path / 'gene_labels.npy', gene_labels)


## Train GASTON network
# Load data generated above
path_to_coords = output_path / 'coords_mat.npy'
path_to_glmpca = output_path / 'glmpca.npy'
S = np.load(path_to_coords)
A = np.load(path_to_glmpca)

# z-score normalize S and A
S_torch, A_torch = neural_net.load_rescale_input_data(S, A)

###################################### 
# NEURAL NET PARAMETERS (USER CAN CHANGE)
# architectures are encoded as list, eg [20,20] means two hidden layers of size 20 hidden neurons
isodepth_arch = [20, 20]  # architecture for isodepth neural network d(x,y) : R^2 -> R 
expression_fn_arch = [20, 20]  # architecture for 1-D expression function h(w) : R -> R^G
num_epochs = 200000  # number of epochs to train NN (NOTE: it is sometimes beneficial to train longer)
checkpoint = 500  # save model after number of epochs = multiple of checkpoint
out_dir = output_path / 'model'  # folder to save model runs
optimizer = "adam"
num_restarts = 30
######################################

seed_list = range(num_restarts)
for seed in seed_list:
    print(f'training neural network for seed {seed}')
    out_dir_seed = f"{out_dir}/rep{seed}"
    os.makedirs(out_dir_seed, exist_ok=True)    
    # Move model to GPU
    mod, loss_list = neural_net.train(
        S_torch, A_torch, S_hidden_list=isodepth_arch, A_hidden_list=expression_fn_arch, 
        epochs=num_epochs, checkpoint=checkpoint, 
        save_dir=out_dir_seed, optim=optimizer, seed=seed, save_final=True)