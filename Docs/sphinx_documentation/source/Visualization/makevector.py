from paraview.util.vtkAlgorithm import *

@smproxy.filter()
@smproperty.input(name="InputDataset", port_index=0)
@smdomain.datatype(dataTypes=["vtkOverlappingAMR"])
class MakeVector(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType="vtkOverlappingAMR", inputType="vtkOverlappingAMR")

    def _find_names(self, fields):
        # This method uses regular expressions to figure out
        # vectors that can be merged. It assumes that the component
        # scalars are named as [xyz]name or name[xyz]
        import re
        array_names = fields.keys()
        to_merge = {}
        matches = {}
        for a in array_names:
            grp = re.search(r'^([xyz])[_-]?([\w]*)', a)
            if grp:
                grp = grp.groups()
                if grp[0] and grp[1]:
                    try:
                        matches[grp[1]][grp[0]] = a
                    except KeyError:
                        matches[grp[1]] = { grp[0] : a }
        for a in array_names:
            grp = re.search(r'([\w]*)[_-]?([xyz])$', a)
            if grp:
                grp = grp.groups()
                if grp[0] and grp[1]:
                    try:
                        matches[grp[0]][grp[1]] = a
                    except KeyError:
                        matches[grp[0]] = { grp[1] : a }
        for key, value in matches.items():
            keys = set(value.keys())
            if (len(keys) == 2 or len(keys) == 3) and 'x' in keys \
                and 'y' in keys:
                try:
                    to_merge[key] = (value['x'], value['y'], value['z'])
                except KeyError:
                    to_merge[key] = (value['x'], value['y'])
        return to_merge

    def RequestData(self, request, inInfoVec, outInfoVec):
        # These modules are needed for the data model and numpy convenience
        # methods.
        from vtkmodules.vtkCommonDataModel import vtkOverlappingAMR
        from vtkmodules.numpy_interface import dataset_adapter as dsa
        from vtkmodules.numpy_interface import algorithms as algs
        # Boilerplate code to get the input and output
        inp = vtkOverlappingAMR.GetData(inInfoVec[0], 0)
        output = vtkOverlappingAMR.GetData(outInfoVec, 0)

        # Pass all of the AMR structure and fields to the output.
        output.ShallowCopy(inp)

        # Create objects that offer convenience APIs that interact
        # with numpy.
        inp = dsa.WrapDataObject(inp)
        output = dsa.WrapDataObject(output)

        # Workhorse. Merge all arrays of for [xyz]name and name[xyz] in
        # point and cell data
        for inp_field, output_field in [(inp.CellData, output.CellData), \
                                        (inp.PointData, output.PointData)]:
            # This method find the arrays names that match the criteria
            to_merge = self._find_names(inp_field)
            for name, arrays in to_merge.items():
                if len(arrays) == 2:
                    # make_vector is really all we need.
                    output_field.append(
                        algs.make_vector(inp_field[arrays[0]],
                                         inp_field[arrays[1]]), name)
                else:
                    output_field.append(
                        algs.make_vector(inp_field[arrays[0]],
                                         inp_field[arrays[1]],
                                         inp_field[arrays[2]]), name)

        return 1
