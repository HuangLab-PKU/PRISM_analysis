import os
from pathlib import Path
from tqdm import tqdm
from pprint import pprint
import pickle
import yaml
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scanpy as sc

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
# workdir 
BASE_DIR = Path(r'G:\spatial_data')
RUN_ID_LIST = []

# analysis dir
RUN_ID = '20250222_combined_analysis_of_pseudo_HCC3D'
analysis_dir = BASE_DIR / 'analysis' / RUN_ID
typ_path = analysis_dir / "cell_typing"

output_path = analysis_dir / "STAGATE"
output_path.mkdir(exist_ok=True)

## define region based on cell composition and gene expression
adata_STAGATE = sc.read_h5ad(output_path / 'adata_STAGATE.h5ad')
with open(os.path.join(analysis_dir / "STAGATE_analysis_params.yaml"), "r") as f:
    STAGATE_params = yaml.load(f, Loader=yaml.FullLoader)
STAGATE_annotate = STAGATE_params['STAGATE_mclust_11_annotate']
pprint(STAGATE_annotate, sort_dicts=False)
map_dict = {str(mclust): str(num+1) for num, (key, mclust_list) in enumerate(STAGATE_annotate.items()) for mclust in mclust_list}
pprint(map_dict)
adata_STAGATE.obs['mclust_annotated'] = adata_STAGATE.obs['mclust_11'].astype(str).map(map_dict).astype(str)
adata_STAGATE.write_h5ad(output_path / 'adata_STAGATE_annotated.h5ad')