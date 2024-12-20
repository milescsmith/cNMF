from importlib.metadata import PackageNotFoundError, version

from loguru import logger

from cnmf.cnmf import cNMF, load_df_from_npz, main, save_df_to_npz
from cnmf.logging import init_logger
from cnmf.preprocess import Preprocess

logger.disable(__package__)
init_logger(verbose=3, save_log=True)
try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"

__all__ = ["Preprocess", "cNMF", "load_df_from_npz", "main", "save_df_to_npz"]
