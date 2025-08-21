# STAGATE analysis
# Loading the Packages
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc

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

# workdir 
BASE_DIR = Path(r'G:\spatial_data\analysis')
RUN_ID = '20250222_combined_analysis_of_pseudo_HCC3D'

# Load one slide exp
base_path = BASE_DIR / f'{RUN_ID}'
data_path = base_path / "segmented"
typ_path = base_path / "cell_typing"
output_path = base_path / "STAGATE"
output_path.mkdir(exist_ok=True)

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

rad_cutoff = 250
# rad_cutoff_path = output_path / f'rad_cutoff_{str(rad_cutoff)}'
rad_cutoff_path = output_path
rad_cutoff_path.mkdir(exist_ok=True)

import STAGATE
mclust_dir = rad_cutoff_path / 'mClust_test'
mclust_dir.mkdir(exist_ok=True)
adata_STAGATE = sc.read_h5ad(rad_cutoff_path / 'adata_STAGATE_light.h5ad')
# for cluster in range(2,15):
cluster = 42
adata_STAGATE = STAGATE.mclust_R(adata_STAGATE, used_obsm='STAGATE', used_savename=f'mclust_{cluster}', num_cluster=cluster)
adata_STAGATE.obs[f'mclust_{cluster}'] = adata_STAGATE.obs[f'mclust_{cluster}'].astype(str)

fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(20, 10))
sc.pl.umap(adata_STAGATE, color=f'mclust_{cluster}', palette=colors, s=0.1, show=False, ax=ax[0])
sc.pl.embedding(adata_STAGATE, basis="spatial", color=f'mclust_{cluster}', palette=colors, s=0.1, show=False, ax=ax[1]) #, legend_loc=False)
plt.tight_layout()
plt.savefig(mclust_dir / f'cluster_{cluster}.png')
plt.close()

adata_STAGATE.write_h5ad(rad_cutoff_path / 'adata_STAGATE_mclust.h5ad')