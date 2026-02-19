import pathlib
import logging
import ROOT as rt

rt.gROOT.SetBatch(True)
rt.EnableThreadSafety()
rt.TH1.AddDirectory(False)

logger = logging.getLogger(__name__)
if not logger.hasHandlers():
    logging.basicConfig(level=logging.DEBUG)


def GetObjectFromFile(file, path):
    """
    Get an object from an already-opened TFile.

    Args:
        file (TFile): An already-opened TFile object.
        path (str): Path to the object in the TFile, delimited by '/'.

    Returns:
        The requested ROOT object.

    Raises:
        ValueError: If any part of the path cannot be found.
    """
    dirs = path.strip().split("/")
    current = file
    for name in dirs:
        next_obj = current.Get(name) or current.FindObject(name)
        if not next_obj:
            raise ValueError(
                f"Object '{name}' not found in '{current.GetName()}' (full path: '{path}')"
            )
        current = next_obj
    return current


def GetObject(filename, path):
    """
    Get a ROOT object from a file by filename.

    Wrapper around GetObjectFromFile that handles opening and closing the file.

    Args:
        filename (str): Full path to the ROOT file.
        path (str): Path to the object in the TFile, delimited by '/'.

    Returns:
        A clone of the requested ROOT object.

    Raises:
        ValueError: If the file cannot be opened or the object is not found.
    """
    file = rt.TFile(filename, "READ")
    if not file or file.IsZombie():
        raise ValueError(f"Could not open ROOT file: {filename}")
    obj = GetObjectFromFile(file, path).Clone()
    file.Close()
    return obj


def CreateOutputDir(path, rename_old=False):
    """
    Ensure the parent directory of the given path exists.

    Args:
        path (str or Path): Path whose parent directory should be created.
        rename_old (bool): If True and a file already exists at path,
                           rename it with a numeric suffix before returning.
    """
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if rename_old and path.exists():
        counter = 1
        while True:
            new_path = path.parent / f"{path.stem}_{counter:04d}{path.suffix}"
            if not new_path.exists():
                path.rename(new_path)
                logger.debug("Renamed existing file %s -> %s", path, new_path)
                break
            counter += 1
