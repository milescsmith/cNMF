from importlib.metadata import PackageNotFoundError, version

from cnmf.cnmf import cNMF, load_df_from_npz, main, save_df_to_npz
from cnmf.preprocess import Preprocess

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"

__all__ = ["Preprocess", "cNMF", "load_df_from_npz", "main", "save_df_to_npz"]
