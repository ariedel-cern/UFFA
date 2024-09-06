import ROOT as rt


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
