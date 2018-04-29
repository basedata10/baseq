def dynamicbin_reader(path):
    import pandas as pd
    return pd.read_table(path, names=["chr", "start", 'absstart', 'end', 'range', 'length', 'GC'], sep=",")