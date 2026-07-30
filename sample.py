from re import S
import numpy as np
import sys
from io import StringIO
import tifffile

sys.path[-1], sys.path[0] = sys.path[0], sys.path[-1]
# from pypne import pnextract
import pypne

image = np.load("./image_100_300_300_pore.npz")["arr_0"]
image_VElems, pn = pypne.pnextract(
    image,
    0.1,
    config_settings={"write_all": True, "name": "pn", "output_path": "./pn_output"},
    verbose=True,
    n_workers=10,
)

image = image_VElems[1:-1, 1:-1, 1:-1]
# tifffile.imwrite("./image_VElems_temp.tif", image)
