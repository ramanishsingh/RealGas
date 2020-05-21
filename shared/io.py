
def read_csv(f_name: str, delimiter: str):
    """
    :param f_name: name of file
    :param delimiter: delimiter separating columns
    """
    with open(f_name) as f:
        keys = next(f).rstrip('\n').split(delimiter)
        data = {key: [] for key in keys}
        for line in f:
            vals = line.rstrip('\n').split(delimiter)
            for i, j in zip(keys, vals):
                data[i].append(j)

    return data
