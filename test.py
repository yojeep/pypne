from re import S
import numpy as np
import sys
from io import StringIO

sys.path[-1], sys.path[0] = sys.path[0], sys.path[-1]
# from pypne import pnextract
import scipy.io as sio
import pypne

image = np.fromfile("./image_100_300_300_pore.raw", dtype=np.uint8).reshape(
    300, 300, 100
)
# image = image[:100, :100, :100]
image_VElems, pn = pypne.pnextract(
    image,
    0.1,
    config_settings={"write_all": True, "name": "pn", "output_path": "./pn_output"},
    verbose=True,
    n_workers=64,
)



for key in pn:
    print(len(pn[key]))
image = image_VElems[1:-1, 1:-1, 1:-1]
# image = image.reshape(300,300,100).tofile("./image_100_300_300_pore_VElems.raw")
