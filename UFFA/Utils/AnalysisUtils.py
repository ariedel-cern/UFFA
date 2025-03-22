import pathlib
import ROOT


def CheckDictEntry(dict, key, type):
    if key not in dict:
        raise ValueError(f"{key} not configured")
    if not isinstance(dict[key], type):
        raise TypeError(f"{key} has incorrect type")


def GetObjectFromFile(file, path):
    """
    Get an object from file
    Args:
        file (TFile): TFile object that it already opend
        path (string): path to the object in the TFile, delimited by /

    Returns:
       Searched object
    """

    # split path into a list
    dirs = path.split("/")

    #
    current_dir = file
    next_dir = None

    for dir in dirs:
        # try to find object in TDirectory
        try:
            next_dir = current_dir.Get(dir)
        except:
            pass
        # if this failed, try to find object in TList
        if not next_dir:
            try:
                next_dir = current_dir.findobject(dir)
            except:
                pass
        if not next_dir:
            raise ValueError(f"Object at {path} not found")
        # reset the pointer
        current_dir = next_dir
    return current_dir


def GetObject(filename, path):
    """
    Get object from file with filename
    Wrapper function for GetObjectFromFile which takes care of opening and closing the file
    Args:
        filename (string): filename with full path
        path (string): path to the object in the TFile, delimited by /

    Returns:
       Searched object
    """
    file = ROOT.TFile(filename, "READ")
    object = GetObjectFromFile(file, path).Clone()
    file.Close()
    return object


def CreateOutputDir(path, rename_old=False):
    """
    Make sure that the parent directory of path exists
    """
    # create parent dir of path
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # if set, check if file at path already exists
    # if so, rename
    if rename_old:
        if path.exists():
            base = path.stem  # File name without extension
            suffix = path.suffix  # File extension (if any)
            parent = path.parent
        else:
            return

        counter = 1
        while True:
            new_name = f"{base}_{counter:04d}{suffix}"
            new_path = parent / new_name
            if not new_path.exists():  # Found a free name
                path.rename(new_path)  # Rename the existing file
                break
            counter += 1
