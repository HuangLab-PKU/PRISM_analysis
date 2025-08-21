# STAGATE analysis
# Loading the Packages
import os
from pathlib import Path
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Data Processing
import scanpy as sc

# Visualization
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

## proprocess of data
combine_adata = sc.read_h5ad(typ_path/'adata_PRISM_HCC_3D_pseudo_st_category_reordered.h5ad')
combine_adata.obs['X'] = combine_adata.obsm['spatial3d'][:,0]
combine_adata.obs['Y'] = combine_adata.obsm['spatial3d'][:,1]
combine_adata.obs['Z'] = combine_adata.obsm['spatial3d'][:,2]
combine_adata.obsm['spatial'] = combine_adata.obs.loc[:, ['X', 'Y']].values


## training of STAGATE and build graph
import STAGATE_pyG as STAGATE
rad_cutoff = 100
rad_cutoff_3d = 100
rad_cutoff_path = output_path / f'rad_cutoff_{str(rad_cutoff)}'
rad_cutoff_path.mkdir(exist_ok=True)

section_order = [f'layer{_}' for _ in range(20)]
STAGATE.Cal_Spatial_Net_3D(combine_adata, rad_cutoff_2D=rad_cutoff, rad_cutoff_Zaxis=rad_cutoff_3d, key_section='layer', section_order=section_order)
combine_adata_STAGATE = STAGATE.train_STAGATE(combine_adata, hidden_dims=[64, 30])
combine_adata_STAGATE.write(rad_cutoff_path / 'combine_adata_STAGATE.h5ad')

combine_adata_STAGATE = sc.read_h5ad(rad_cutoff_path / 'combine_adata_STAGATE.h5ad')
sc.pp.neighbors(combine_adata_STAGATE, n_neighbors=50, use_rep='STAGATE')
sc.tl.umap(combine_adata_STAGATE)
combine_adata_STAGATE.write_h5ad(rad_cutoff_path / 'combine_adata_STAGATE.h5ad')