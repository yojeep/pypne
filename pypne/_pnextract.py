import numpy as np
from pathlib import Path
import os
import sys
from io import StringIO
from contextlib import contextmanager, nullcontext, redirect_stdout
from .libcpp import pypne_cpp


@contextmanager
def suppress_stdout():
    """
    mute stdout(Python + C/C++)
    """
    original_stdout_fd = sys.stdout.fileno()
    original_stdout = os.dup(original_stdout_fd)
    null_fd = os.open(os.devnull, os.O_WRONLY)

    try:
        os.dup2(null_fd, original_stdout_fd)
        with open(os.devnull, "w") as f, redirect_stdout(f):
            yield  # 执行代码
    finally:
        os.dup2(original_stdout, original_stdout_fd)
        os.close(null_fd)
        os.close(original_stdout)


_true_set = {"yes", "true", "t", "y", "1"}
_false_set = {"no", "false", "f", "n", "0"}


def str2bool(value, raise_exc=False):
    if isinstance(value, bool):
        return value
    if (
        isinstance(value, str)
        or sys.version_info[0] < 3
        and isinstance(value, basestring)
    ):
        value = value.lower()
        if value in _true_set:
            return True
        if value in _false_set:
            return False

    if raise_exc:
        raise ValueError('Expected "%s"' % '", "'.join(_true_set | _false_set))
    return None


def pnextract(image, resolution=1.0, config_settings=None, verbose=False,n_workers=1):
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


    config_settings: a dictionary containing the following keys:
    write_Statoil:false,
    write_radius:false,
    write_elements:false,
    write_hierarchy:false,
    write_throatHierarchy:false,
    write_vtkNetwork:false,
    write_throats:false,
    write_poreMaxBalls:false,
    write_throatMaxBalls:false,


    output_path : path to output file, using default value "./pn(with desired suffix)"

    minRPore: minimum radius of pore, using default value _minRp=min(1.25, avgR*0.25)+0.5

    medialSurfaceSettings: medial surface settings, using the following default values:
    _clipROutx=0.05;
        _clipROutyz=0.98;
        _midRf=0.7;
        _MSNoise=1.*abs(_minRp)+1.;
        _lenNf=0.6;
        _vmvRadRelNf=1.1;
        _nRSmoothing=3;
        _RCorsnf=0.15;
        _RCorsn=abs(_minRp);

    If you wants to set medialSurfaceSettings, you should use config_settings like this:
    config_settings['medialSurfaceSettings'] = "_clipROutx _clipROutyz _midRf _MSNoise _lenNf _vmvRadRelNf _nRSmoothing _RCorsnf _RCorsn"
    change the arguments to values you want.
    """
    Path_cwd = Path.cwd()
    default_config = {
        "write_Statoil": False,
        "write_radius": False,
        "write_elements": False,
        "write_hierarchy": False,
        "write_throatHierarchy": False,
        "write_vtkNetwork": False,
        "write_throats": False,
        "write_poreMaxBalls": False,
        "write_throatMaxBalls": False,
        "write_all": False,
        "output_path": (Path_cwd).resolve(),
        "name": "pn",
        "minRPore": None,
        "medialSurfaceSettings": None,
    }

    if config_settings is not None:
        default_config.update(config_settings)
    default_config = {k: v for k, v in default_config.items() if v is not None}
    default_config = {str(k): str(v) for k, v in default_config.items()}
    if str2bool(default_config["write_all"]):
        for k in default_config:
            if k.startswith("write_"):
                default_config[k] = "true"
    anything2write = any(
        str2bool(default_config[k]) for k in default_config if k.startswith("write_")
    )
    if anything2write:
        os.makedirs(default_config["output_path"], exist_ok=True)
        default_config["output_path"] = os.path.join(
            default_config["output_path"], default_config["name"]
        )
        default_config.pop("name")
    else:
        default_config.pop("output_path")
        default_config.pop("name")
    image = image.astype(np.uint8)
    nz, ny, nx = image.shape
    # 直接根据 verbose 决定是否使用 suppress_stdout
    with suppress_stdout() if not verbose else nullcontext():
        res = pypne_cpp.pnextract(
            nx, ny, nz, resolution, image.reshape(-1), default_config.copy(), n_workers
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

    if str2bool(default_config["write_elements"]):
        image_VElems.astype(np.int32, copy=False).tofile(
            default_config["output_path"] + "_VElems.raw"
        )
    return image_VElems, pn
