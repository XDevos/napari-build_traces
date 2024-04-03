#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import uuid
import matplotlib

matplotlib.use("Agg")

from datetime import datetime
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.patches import Polygon
from modules.module import Module
from scipy.spatial import KDTree
from stardist import random_label_cmap

sys.path.append("../../..")

from chromapylot.core.core_types import AnalysisType, DataType
from chromapylot.core.data_manager import load_json, save_ecsv
from chromapylot.parameters.pipeline_params import PipelineParams
from chromapylot.parameters.acquisition_params import AcquisitionParams
from chromapylot.parameters.matrix_params import MatrixParams
from chromapylot.core.run_args import RunArgs
from chromapylot.modules.build_trace import BuildTrace3DModule

import json
import pathlib
from typing import TYPE_CHECKING
from napari.layers import Layer, Image
from napari.types import LayerDataTuple
from magicgui import magic_factory
from magicgui.widgets import FunctionGui
import numpy as np
import scipy.optimize as spo

if TYPE_CHECKING:
    import napari


def _method_choices(wdg):
    # Update operations
    opts = ["masking", "clustering"]
    return opts


def on_init(widget: FunctionGui):
    """called each time widget_factory creates a new widget."""

    @widget.save.changed.connect
    def on_save_changed(save_val: bool):
        if save_val:
            widget.save_path.visible = True
        else:
            widget.save_path.visible = False


@magic_factory(
    widget_init=on_init,
    pixel_size_xy={"value": 0.1},
    pixel_size_z={"value": 0.25},
    z_binning={"value": 2},
    method={"choices": _method_choices},
    mask_to_process={"value": None},
    kdtree_distance_threshold={"value": 1.0},
    call_button="Build traces",
    save={
        "widget_type": "CheckBox",
        "value": False,
        "name": "save",
        "text": "Save configuration",
    },
    save_path={"mode": "d", "visible": False},
)
def do_build_trace(
    localization_path: pathlib.Path,
    output_folder: pathlib.Path,
    pixel_size_xy: float,
    pixel_size_z: float,
    z_binning: int,
    method: str,
    mask_to_process: str,
    kdtree_distance_threshold: float,
    save: bool,
    save_path: pathlib.Path,
) -> LayerDataTuple:


    # load a parameters file from local directory
    param_path = pathlib.Path("/home/xdevos/Repositories/pyHiM/src/toolbox/parameter_file/parameters.json")
    with open(param_path) as f:
        param = json.load(f)

    # update the parameters
    param["common"]["acquisition"]["pixel_size_xy"] = pixel_size_xy    
    param["common"]["acquisition"]["pixel_size_z"] = pixel_size_z
    param["common"]["acquisition"]["zBinning"] = z_binning
    param["common"]["buildsPWDmatrix"]["tracing_method"] = list(method)
    param["common"]["buildsPWDmatrix"]["KDtree_distance_threshold_mum"] = kdtree_distance_threshold    

    
    # INITIALIZATION
    pipe_params = PipelineParams(param, AnalysisType.TRACE)
    mod = BuildTrace3DModule(pipe_params.acquisition, pipe_params.matrix)
    ref_file = None
    ref_files = ref_file.split(",") if ref_file else []

    # MODULE EXECUTION
    mod.load_reference_data(ref_files)
    input_data = mod.load_data(localization_path, None)
    output_data = mod.run(input_data)
    mod.save_data(output_data, output_folder, localization_path)
    #############################################################################################
    #############################################################################################

    # if save:
    #     data_to_print = {"operationTODO": operationTODO, "image": layer0_name}
    #     with open("json_data_todelete.json", "w") as outfile:
    #         json.dump(data_to_print, outfile)
    #     print("Configuration save as json file")

    return
