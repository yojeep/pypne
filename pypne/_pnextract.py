from dataclasses import asdict, dataclass, field
from typing import Optional
import numpy as np
from pathlib import Path
import os
import sys
from io import StringIO
from contextlib import contextmanager, nullcontext, redirect_stdout
from .libcpp import pypne_cpp


@dataclass
class PnConfig:
    # 结构定义与默认值在同一个地方，一目了然
    write_Statoil: bool = False
    write_elements: bool = False
    write_all: bool = False
    output_path: Path = field(default_factory=lambda: Path.cwd())
    name: str = "pn"
    minRp: Optional[float] = None
    nBP6: int = 2
    clipROutx: Optional[int] = None
    clipROutyz: Optional[int] = None
    midRf: Optional[float] = None
    MSNoise: Optional[float] = None
    lenNf: Optional[int] = None
    vmvRadRelNf: Optional[float] = None
    nRSmoothing: Optional[int] = None
    RCorsnf: Optional[float] = None
    RCorsn: Optional[float] = None


# 使用时极其清爽，不需要重复写一堆键名
cfg = PnConfig()


@contextmanager
def suppress_stdout():
    sys.stdout.flush()
    saved = os.dup(1)
    try:
        with Path(os.devnull).open("w") as devnull:
            os.dup2(devnull.fileno(), 1)
            with redirect_stdout(devnull):
                yield
    finally:
        os.dup2(saved, 1)
        os.close(saved)


def pnextract(
    image,
    resolution=1.0,
    config_settings=None,
    verbose=False,
    n_workers=1,
):
    """
    image : 3D numpy array of binary image data,0 is the value to be extracted
    resolution : resolution of the image, default is 1.0
    verbose : whether to print the progress of the algorithm, default is False

    return :
    1. extracted image, which shape is (nz+2, ny+2, nx+2)
    2. pore network(pn).
    pn is a dict containing the following keys:
    'pore._id'
    'pore.x'
    'pore.y'
    'pore.z'
    'pore.connection_number'
    'pore.volume'
    'pore.radius'
    'pore.shape_factor'

    'throat._id'
    'throat.pore_1_index'
    'throat.pore_2_index'
    'throat.radius'
    'throat.shape_factor'
    'throat.total_length'
    'throat.conduit_lengths_pore1'
    'throat.conduit_lengths_pore2'
    'throat.length'
    'throat.volume'
    'throat.clay_volume'

    config_settings: a dictionary containing the following keys:
    write_Statoil:false,
    write_elements:false,
    write_all:false,
    minRp: minimum radius of pore, using default value _minRp=min(1.25, avgR*0.25)+0.5
    clipROutx=0.05;
    clipROutyz=0.98;
    midRf=0.7;
    MSNoise=1.*abs(_minRp)+1.;
    lenNf=0.6;
    vmvRadRelNf=1.1;
    nRSmoothing=3;
    RCorsnf=0.15;
    RCorsn=abs(_minRp);
    output_path : path to output file, using default value "./pn(with desired suffix)"
    """
    cfg = asdict(PnConfig())

    if config_settings:
        unsupported_keys = config_settings.keys() - cfg.keys()
        if unsupported_keys:
            raise ValueError(
                f"Unsupported keys: {unsupported_keys}, supported keys are {set(cfg.keys())}"
            )
        cfg.update(config_settings)

    if cfg["write_all"]:
        for k in cfg:
            if k.startswith("write_"):
                cfg[k] = True

    if any(v for k, v in cfg.items() if k.startswith("write_")):
        cfg["output_path"] = Path(cfg["output_path"]).resolve()
        Path(cfg["output_path"]).mkdir(parents=True, exist_ok=True)
        cfg["output_path"] = Path(cfg["output_path"]) / cfg["name"]
    cfg["output_path"] = str(cfg["output_path"])

    image = image.astype(bool, copy=False)
    nz, ny, nx = image.shape
    # 直接根据 verbose 决定是否使用 suppress_stdout
    n_workers = (
        os.cpu_count() + 1 + n_workers
        if n_workers < 0
        else max(1, min(n_workers, os.cpu_count()))
    )
    with suppress_stdout() if not verbose else nullcontext():
        res = pypne_cpp.pnextract(
            nz, ny, nx, resolution, image.reshape(-1), cfg, n_workers
        )

    image_VElems = res["VElems"].reshape(nz + 2, ny + 2, nx + 2)
    pn = res["pn"]
    link1 = pn["link1"]
    link2 = pn["link2"]
    node1 = pn["node1"]
    node2 = pn["node2"]
    link1_arr = np.genfromtxt(
        StringIO(link1),
        delimiter=None,
        skip_header=1,
        usecols=(0, 1, 2, 3, 4, 5),
        dtype=[
            ("throat__id", "int32"),
            ("throat_pore_1_index", "int32"),
            ("throat_pore_2_index", "int32"),
            ("throat_radius", "float32"),
            ("throat_shape_factor", "float32"),
            ("throat_total_length", "float32"),
        ],
    )

    link2_arr = np.genfromtxt(
        StringIO(link2),
        delimiter=None,
        usecols=(0, 1, 2, 3, 4, 5, 6, 7),
        dtype=[
            ("throat__id", "int32"),
            ("throat_pore_1_index", "int32"),
            ("throat_pore_2_index", "int32"),
            ("throat_conduit_lengths_pore1", "float32"),
            ("throat_conduit_lengths_pore2", "float32"),
            ("throat_length", "float32"),
            ("throat_volume", "float32"),
            ("throat_clay_volume", "float32"),
        ],
    )

    node1_arr = np.genfromtxt(
        StringIO(node1),
        delimiter=None,
        skip_header=1,
        usecols=(0, 1, 2, 3, 4),
        dtype=[
            ("pore__id", "int32"),
            ("pore_x", "float32"),
            ("pore_y", "float32"),
            ("pore_z", "float32"),
            ("pore_connection_number", "int32"),
        ],
    )

    node2_arr = np.genfromtxt(
        StringIO(node2),
        delimiter=None,
        usecols=(0, 1, 2, 3, 4),
        dtype=[
            ("pore__id", "int32"),
            ("pore_volume", "float32"),
            ("pore_radius", "float32"),
            ("pore_shape_factor", "float32"),
            ("pore_clay_volume", "float32"),
        ],
    )

    pn["pore._id"] = node1_arr["pore__id"]
    pn["pore.x"] = node1_arr["pore_x"]
    pn["pore.y"] = node1_arr["pore_y"]
    pn["pore.z"] = node1_arr["pore_z"]
    pn["pore.connection_number"] = node1_arr["pore_connection_number"]
    pn["pore.volume"] = node2_arr["pore_volume"]
    pn["pore.radius"] = node2_arr["pore_radius"]
    pn["pore.shape_factor"] = node2_arr["pore_shape_factor"]
    pn["pore.clay_volume"] = node2_arr["pore_clay_volume"]
    pn["throat._id"] = link1_arr["throat__id"]
    pn["throat.pore_1_index"] = link1_arr["throat_pore_1_index"]
    pn["throat.pore_2_index"] = link1_arr["throat_pore_2_index"]
    pn["throat.radius"] = link1_arr["throat_radius"]
    pn["throat.shape_factor"] = link1_arr["throat_shape_factor"]
    pn["throat.total_length"] = link1_arr["throat_total_length"]
    pn["throat.conduit_lengths_pore1"] = link2_arr["throat_conduit_lengths_pore1"]
    pn["throat.conduit_lengths_pore2"] = link2_arr["throat_conduit_lengths_pore2"]
    pn["throat.length"] = link2_arr["throat_length"]
    pn["throat.volume"] = link2_arr["throat_volume"]
    pn["throat.clay_volume"] = link2_arr["throat_clay_volume"]

    if cfg["write_elements"]:
        image_VElems.astype(np.int32, copy=False).tofile(
            f"{cfg['output_path']}_VElems_{image_VElems.shape[2]}x{image_VElems.shape[1]}x{image_VElems.shape[0]}_s32_le.raw"
        )
    return image_VElems, pn
