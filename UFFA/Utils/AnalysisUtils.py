def CheckDictEntry(dict, key, type):
    if key not in dict:
        raise ValueError(f"{key} not configured")
    if not isinstance(dict[key], type):
        raise TypeError(f"{key} has incorrect type")
