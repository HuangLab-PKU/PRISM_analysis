# STAGATE analysis
# Loading the Packages
import os
from pathlib import Path
import pickle
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.signal import argrelextrema
from scipy.signal import find_peaks

import seaborn as sns
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

# Load one slide exp
base_path = Path(r'F:\spatial_data\analysis\20250222_combined_analysis_of_pseudo_HCC3D')
STAGATE_path = base_path / "STAGATE"
STAGATE_path.mkdir(exist_ok=True)

## proprocess of data
combine_adata = sc.read_h5ad(r'F:\spatial_data\processed\PRISM_publication\PRISM30_pseudo3D_HCC_20_slides\cell_typing\combine_adata_st.h5ad')
adata_direct = sc.read_h5ad(r'F:\spatial_data\processed\PRISM_publication\PRISM30_pseudo3D_HCC_20_slides\cell_typing\adata.h5ad')
adata = adata_direct[adata_direct.obs.index.isin(combine_adata.obs.index)]
adata.obs = combine_adata.obs

# format for later analysis
adata.obs = adata.obs.rename(columns={'X_pos':'X', 'Y_pos':'Y'})
adata.obsm['spatial'] = adata.obs.loc[:, ['X', 'Y']].values

## clustering of ROIs
type_colormap = {
    'Liver':(1,0.392,0),
    'Tumor':(0.751,0.491,0),
    'Endo':(1,0,1),
    'Ep':(0,1,0),
    'CAF':(0,0,1),
    'DC':(1,0.259,0),
    'Mait':(1,0,0.434),
    'Mast':(1,0,0),
    'Monocyte':(0,0.471,1),
    'Neutrophil':(0.6,0.6,0),
    'Macrophage':(0.2,0.6,0),
    'CD4+':(0.5,0.5,0.5),
    'CD8+':(1,0.8,0),
    'T_reg':(0,1,0.672),
    'T_proliferation':(0,1,0.636),
    'B':(0,1,1),
    'NK':(1,0,0),
    'other':(0.9,0.9,0.9),
}
colors = [_ for _ in type_colormap.values()][:-1]
colors[5] = (0,0,0)

## training of STAGATE and build graph
import STAGATE
import STAGATE_pyG

def Stats_Spatial_Net_manual(adata, outpath=None):
    import matplotlib.pyplot as plt
    Num_edge = adata.uns['Spatial_Net']['Cell1'].shape[0]
    Mean_edge = Num_edge/adata.shape[0]
    plot_df = pd.value_counts(pd.value_counts(adata.uns['Spatial_Net']['Cell1']))
    plot_df = plot_df/adata.shape[0]
    fig, ax = plt.subplots(figsize=[3,2])
    plt.ylabel('Percentage')
    plt.xlabel('')
    plt.title('Number of Neighbors (Mean=%.2f)'%Mean_edge)
    ax.bar(plot_df.index, plot_df)
    if not outpath is None:
        plt.savefig(outpath)
        plt.close()

# for rad_cutoff in tqdm([50,100,150,200,250,300,400]):
# for rad_cutoff in tqdm([150,200,250,300,400]):
for rad_cutoff in tqdm([250]):
    current_path = STAGATE_path / f'rad_cutoff_{str(rad_cutoff)}'
    current_path.mkdir(exist_ok=True)

    section_order = [f'slice{_}' for _ in range(20)]
    STAGATE.Cal_Spatial_Net_3D(adata, rad_cutoff_2D=rad_cutoff, rad_cutoff_Zaxis=rad_cutoff, key_section='slice', section_order=section_order)
    Stats_Spatial_Net_manual(adata, outpath=current_path / 'stats_spatial_net.png')
    adata.write(current_path / 'adata_STAGATE.h5ad')

    adata_STAGATE = sc.read_h5ad(current_path / 'adata_STAGATE.h5ad')
    adata_STAGATE = STAGATE_pyG.train_STAGATE(adata, hidden_dims=[64, 24])
    adata_STAGATE.write(current_path / 'adata_STAGATE.h5ad')

    adata_STAGATE = sc.read_h5ad(current_path / 'adata_STAGATE.h5ad')
    sc.pp.neighbors(adata_STAGATE, n_neighbors=50, use_rep='STAGATE')
    sc.tl.umap(adata_STAGATE)
    adata_STAGATE.write_h5ad(current_path / 'adata_STAGATE.h5ad')

    mclust_dir = current_path / 'mClust_test'
    mclust_dir.mkdir(exist_ok=True)
    for cluster in range(2, 20):
        adata_STAGATE = sc.read_h5ad(current_path / 'adata_STAGATE.h5ad')
        adata_STAGATE = STAGATE.mclust_R(adata_STAGATE, used_obsm='STAGATE', used_savename=f'mclust_{cluster}', num_cluster=cluster)
        adata_STAGATE.obs[f'mclust_{cluster}'] = adata_STAGATE.obs[f'mclust_{cluster}'].astype(str)
        adata_STAGATE.write(current_path / 'adata_STAGATE.h5ad')
        
        fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
        sc.pl.umap(adata_STAGATE, color=f'mclust_{cluster}', palette=colors, s=1, show=False, ax=ax[0])
        sc.pl.embedding(adata_STAGATE, basis="spatial", color=f'mclust_{cluster}', palette=colors, s=1, show=False, ax=ax[1]) #, legend_loc=False)
        plt.tight_layout()
        plt.savefig(mclust_dir / f'cluster_{cluster}.png')
        plt.close()