
# pypne – Python wrapper for pore-network extraction

**pypne** is a Python interface to the high-performance C++ pore-network extraction code originally developed in [pnextract](https://github.com/aliraeini/pnextract). It enables seamless integration of pore-network modeling into scientific Python workflows, allowing users to extract topologically and geometrically accurate pore networks directly from 3D binary images of porous media.

This package wraps the core `pnextract` algorithm using pybind11, providing:
- Easy-to-use Python API with NumPy compatibility
- Full control over algorithm parameters and output options
- Parallel execution support via configurable worker threads

The extracted pore network can be used for downstream simulations (e.g., with [pnflow](https://github.com/aliraeini/pnflow)) or analysis of transport properties in porous materials.

---

## Installation

Precompiled wheels may be available for Windows and Linux

```bash
pip install pypne
```

Install from source, ensure you have a compatible C++ compiler and Python ≥3.8. Then

```bash
pip install git+https://github.com/iPMLab/pypne.git
pip install .
```


---

## Quick Start

```python
import numpy as np
from pypne import pnextract

# Load or generate a 3D binary image (0 = pore, 1 = solid)
image = np.load("./image_100_300_300_pore.npz")["arr_0"]

# Extract pore network
velems, pn = pnextract(
    image,
    resolution=1.0,
    verbose=True,
    n_workers=4
)

print("Number of pores:", len(pn["pore._id"]))
print("Pore radii:", pn["pore.radius"])
```

---

## API Reference

### `pnextract(image, resolution=1.0, config_settings=None, verbose=False, n_workers=1)`

Extracts a pore network from a 3D binary image.

#### Parameters:
- **`image`** (`np.ndarray`, shape `(nz, ny, nx)`, dtype `uint8`):  
  Binary 3D image where `0` denotes pore space and `1` denotes solid.
- **`resolution`** (`float`, optional):  
  Physical voxel size (e.g., in micrometers). Default: `1.0`.
- **`config_settings`** (`dict`, optional):  
  Configuration dictionary controlling output files and algorithm behavior. See below.
- **`verbose`** (`bool`, optional):  
  If `True`, prints C++ progress logs. Default: `False`.
- **`n_workers`** (`int`, optional):  
  Number of threads for parallel processing. Use `≤0` to auto-detect CPU count. Default: `1`.

#### Returns:
- **`velems`** (`np.ndarray`, shape `(nz+2, ny+2, nx+2)`):  
  Labeled void elements image (extended by 1 voxel on each side).
- **`pn`** (`dict`):  
  Pore network data with the following keys:

##### Pore properties:
- `'pore._id'`
- `'pore.x'`, `'pore.y'`, `'pore.z'`
- `'pore.connection_number'`
- `'pore.volume'`
- `'pore.radius'`
- `'pore.shape_factor'`
- `'pore.clay_volume'`

##### Throat properties:
- `'throat._id'`
- `'throat.pore_1_index'`, `'throat.pore_2_index'`
- `'throat.radius'`
- `'throat.shape_factor'`
- `'throat.total_length'`
- `'throat.conduit_lengths_pore1'`, `'throat.conduit_lengths_pore2'`
- `'throat.length'`
- `'throat.volume'`
- `'throat.clay_volume'`

---

### Configuration Settings (`config_settings`)

You can customize behavior and output via the `config_settings` dict:

```python
config = {
    "write_Statoil": False,
    "write_radius": False,
    "write_elements": False,
    "write_hierarchy": False,
    "write_throatHierarchy": False,
    "write_vtkNetwork": False,
    "write_throats": False,
    "write_poreMaxBalls": False,
    "write_throatMaxBalls": False,
    "write_all": False,          # if True, enables all write_* flags
    "output_path": "./results",
    "name": "my_network",        # base filename for outputs
    "minRPore": 1.5,             # minimum pore radius (optional)
    "medialSurfaceSettings": "_clipROutx _clipROutyz _midRf _MSNoise _lenNf _vmvRadRelNf _nRSmoothing _RCorsnf _RCorsn"
}
```

> **Note**: If any `write_*` option is enabled, `output_path` and `name` must be provided (or will default to current directory and `"pn"`).

For `medialSurfaceSettings`, provide a space-separated string of 9 float values corresponding to:
```
_clipROutx _clipROutyz _midRf _MSNoise _lenNf _vmvRadRelNf _nRSmoothing _RCorsnf _RCorsn
```

---

## Output Files (when enabled)

When writing is enabled, the following files may be generated (prefix = `output_path/name`):
- `{prefix}.dat` – Statoil format network
- `{prefix}_radius.raw` – Radius map
- `{prefix}_VElems.raw` – Void element labels (int32, Fortran order)
- `{prefix}_hierarchy.txt` – Pore hierarchy
- `{prefix}_throatHierarchy.txt` – Throat hierarchy
- `{prefix}.vtp` – VTK network for visualization
- `{prefix}_throats.dat` – Throat list
- `{prefix}_poreMaxBalls.txt`, `{prefix}_throatMaxBalls.txt` – Maximal balls

---

## References

This implementation is based on the algorithms described in:

- [Pore-network extraction from micro-computerized-tomography images](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.80.036307)  
  *Phys. Rev. E 80, 036307 (2009)*

- [Generalized network modeling: Network extraction as a coarse-scale discretization of the void space of porous media](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.96.013312)  
  *Phys. Rev. E 96, 013312 (2017)*

Original C++ code: [pnextract](https://github.com/aliraeini/pnextract)

Please cite these papers if you use this code in your research.

---

## Contact

For questions or support, contact:

- Jinping Yang – [12227053@zju.edu.cn](mailto:12227053@zju.edu.cn)
- Qingyang Lin – [qingyan_lin@zju.edu.cn](mailto:qingyan_lin@zju.edu.cn)
- Mingliang Qu – [mingliangqu@zju.edu.cn](mailto:mingliangqu@zju.edu.cn)
