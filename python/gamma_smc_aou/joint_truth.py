"""Small dependency-free native accelerator for exact pairwise tree ages."""
import ctypes
from functools import lru_cache
from pathlib import Path
import numpy as np


@lru_cache(maxsize=1)
def library():
    path=Path(__file__).resolve().parents[2]/"bin/joint_pair_tmrca.so"
    native=ctypes.CDLL(str(path))
    function=native.joint_pair_tmrca
    function.argtypes=[np.ctypeslib.ndpointer(dtype=np.int32,flags="C_CONTIGUOUS"),
                       np.ctypeslib.ndpointer(dtype=np.float64,flags="C_CONTIGUOUS"),
                       np.ctypeslib.ndpointer(dtype=np.int32,flags="C_CONTIGUOUS"),
                       ctypes.c_int64,ctypes.c_int64,
                       np.ctypeslib.ndpointer(dtype=np.float64,flags="C_CONTIGUOUS")]
    function.restype=ctypes.c_int
    return function


def pair_ages(tree,node_pairs,node_times):
    pairs=np.ascontiguousarray(node_pairs,dtype=np.int32)
    times=np.ascontiguousarray(node_times,dtype=np.float64)
    parents=np.ascontiguousarray(tree.parent_array,dtype=np.int32)
    result=np.empty(len(pairs),dtype=np.float64)
    code=library()(parents,times,pairs,len(pairs),len(times),result)
    if code:raise ValueError(f"Invalid pair/tree for exact MRCA extraction: {code}")
    return result
