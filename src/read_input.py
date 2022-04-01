import os


DEFAULT_DICT = {'MAGNETIC_ATOM': 'Fe', 'MAX_T': 1600}


def try_float(item : str):
    try:
        return float(item)
    except ValueError:
        return item


def read_input(input_path: str) -> dict:
    read_path = os.path.join(input_path, 'INPUT')
    with open(read_path, 'r') as f:
        in_list = [st.strip().replace(" ", "").split(':')
                   for st in f.readlines()]
    # transform input dictionary into list
    out_dict = {item[0]: try_float(item[1]) for item in in_list}
    return out_dict


def update_defaults(input_path: str, default_dict: dict) -> dict:
    input_dict = read_input(input_path)
    return default_dict.update(input_dict)


if __name__ == '__main__':
    read_input(os.getcwd())
