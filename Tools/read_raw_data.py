from glob import glob
import numpy as np


def read_data(plt_file):
    '''

    This function reads the raw (i.e. not averaged to cell centers) data 
    from a WarpX plt file. The plt file must have been written with the
    raw_data option turned on, so that it contains a raw_data sub-directory.
    This is only really useful for single-level data.

    Arguments:

        plt_file : An AMReX plt_file file. Must contain a raw_data directory.
    
    Returns:

        A data dictionary where the keys are field name strings and the values 
        are numpy arrays. The shape of the numpy arrays will cover the entire 
        simulation domain.

    Example:

        >>> data = read_data("plt00016")
        >>> print(data.keys())
        >>> print(data['Ex'].shape)

    '''
    raw_file = plt_file + "/raw_fields/Level_0/"
    field_names = _get_field_names(raw_file)
    
    data = {}
    for field in field_names:
        data[field] = _read_field(raw_file, field)

    return data


def _get_field_names(raw_file):
    header_files = glob(raw_file + "*_H")
    return [hf.split("/")[-1][:-2] for hf in header_files]


def _string_to_numpy_array(s):
    return np.array([int(v) for v in s[1:-1].split(",")], dtype=np.int64)


def _line_to_numpy_arrays(line):
    lo_corner = _string_to_numpy_array(line[0][1:])
    hi_corner = _string_to_numpy_array(line[1][:])
    node_type = _string_to_numpy_array(line[2][:-1]) 
    return lo_corner, hi_corner, node_type


def _read_header(raw_file, field):
    header_file = raw_file + field + "_H"
    with open(header_file, "r") as f:

        # skip the first five lines
        for _ in range(5):
            f.readline()

        # read boxes
        boxes = []
        for line in f:
            clean_line = line.strip().split()
            if clean_line == [')']:
                break
            lo_corner, hi_corner, node_type = _line_to_numpy_arrays(clean_line)
            boxes.append((lo_corner, hi_corner, node_type))

        # read the file and offset position for the corresponding box
        file_names = []
        offsets = []
        for line in f:
            if line.startswith("FabOnDisk:"):
                clean_line = line.strip().split()
                file_names.append(clean_line[1])
                offsets.append(int(clean_line[2]))

    return boxes, file_names, offsets


def _combine_boxes(boxes):
    lo_corners, hi_corners = zip(*[(box[0], box[1]) for box in boxes])
    domain_lo = np.min(lo_corners, axis=0)
    domain_hi = np.max(hi_corners, axis=0)
    return domain_lo, domain_hi


def _read_field(raw_file, field_name):

    boxes, file_names, offsets = _read_header(raw_file, field_name)

    lo, hi = _combine_boxes(boxes)
    data = np.zeros(hi - lo + 1)

    for box, fn, offset in zip(boxes, file_names, offsets):
        lo = box[0]
        hi = box[1]
        shape = hi - lo + 1
        with open(raw_file + fn, "rb") as f:
            f.seek(offset)
            f.readline()  # always skip the first line
            arr = np.fromfile(f, 'float64', np.product(shape))
            arr = arr.reshape(shape, order='F')
            data[lo[0]:hi[0]+1, lo[1]:hi[1]+1, lo[2]:hi[2]+1] = arr

    return data


if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    data = read_data("plt00016")

    # print the shapes of all the fields
    for name, vals in data.items():
        print(name, vals.shape)
    
    # make a projection along the z-axis of the 'By' field
    plt.pcolormesh(data['By'].sum(axis=2))
    plt.savefig('By')
