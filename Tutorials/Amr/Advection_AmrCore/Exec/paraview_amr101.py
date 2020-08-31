# based on traces generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *

import subprocess
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--frame_rate', type=int, default=15, help="Frame rate for generating movies, i.e. number of plots per second in the movie.")
parser.add_argument('-r', '--resolution', type=int, default=1024, help="(Square) resolution of output movie.")
parser.add_argument('-d', '--spacedim', type=int, default=3, help="Dimensionality of the problem: 2 or 3")
args = parser.parse_args()

def generate_movie_3D(AllPlotFiles):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
    plt00.CellArrayStatus = []

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Properties modified on plt00
    plt00.CellArrayStatus = ['phi']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    renderView1.ViewSize = [1200, 1200]

    # get layout
    layout1 = GetLayout()

    # show data in view
    plt00Display = Show(plt00, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00Display.Representation = 'Outline'
    plt00Display.ColorArrayName = [None, '']
    plt00Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00Display.SelectOrientationVectors = 'None'
    plt00Display.ScaleFactor = 0.1
    plt00Display.SelectScaleArray = 'None'
    plt00Display.GlyphType = 'Arrow'
    plt00Display.GlyphTableIndexArray = 'None'
    plt00Display.GaussianRadius = 0.005
    plt00Display.SetScaleArray = [None, '']
    plt00Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00Display.OpacityArray = [None, '']
    plt00Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00Display.PolarAxes = 'PolarAxesRepresentation'
    plt00Display.ScalarOpacityUnitDistance = 0.030761993184097912

    # reset view to fit data
    renderView1.ResetCamera()

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # change solid color
    plt00Display.AmbientColor = [0.0, 1.0, 0.0]
    plt00Display.DiffuseColor = [0.0, 1.0, 0.0]

    # create a new 'Slice'
    slice1 = Slice(Input=plt00)
    slice1.SliceType = 'Plane'
    slice1.HyperTreeGridSlicer = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.5, 0.5, 0.0625]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    slice1.HyperTreeGridSlicer.Origin = [0.5, 0.5, 0.0625]

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1.SliceType)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = [None, '']
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.1
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 0.005
    slice1Display.SetScaleArray = [None, '']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = [None, '']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(slice1Display, ('FIELD', 'vtkBlockColors'))

    # get color transfer function/color map for 'vtkBlockColors'
    vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')

    # get opacity transfer function/opacity map for 'vtkBlockColors'
    vtkBlockColorsPWF = GetOpacityTransferFunction('vtkBlockColors')

    # set scalar coloring
    ColorBy(slice1Display, ('CELLS', 'phi'))

    # rescale color and/or opacity maps used to include current data range
    slice1Display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'phi'
    phiLUT = GetColorTransferFunction('phi')

    # get opacity transfer function/opacity map for 'phi'
    phiPWF = GetOpacityTransferFunction('phi')

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # get color legend/bar for phiLUT in view renderView1
    phiLUTColorBar = GetScalarBar(phiLUT, renderView1)

    # change scalar bar placement
    phiLUTColorBar.WindowLocation = 'AnyLocation'
    phiLUTColorBar.Position = [0, 0.75]
    phiLUTColorBar.ScalarBarLength = 0.2

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.5, 0.5, 2.3291959654285184]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0625]
    renderView1.CameraParallelScale = 0.7098635432250342

    # save animation
    output_movie_base = "amr101_3D"
    output_movie = output_movie_base + ".avi"
    SaveAnimation(output_movie,
                  renderView1,
                  ImageResolution=[1200, 1200],
                  FrameRate=args.frame_rate,
                  FrameWindow=[0, len(AllPlotFiles)-1])

    return output_movie_base, output_movie

def generate_movie_2D(AllPlotFiles):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
    plt00.CellArrayStatus = []

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # Properties modified on plt00
    plt00.CellArrayStatus = ['phi']

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1309, 923]
    renderView1.ViewSize = [1200, 1200]

    # get layout
    layout1 = GetLayout()

    # show data in view
    plt00Display = Show(plt00, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00Display.Representation = 'Outline'
    plt00Display.ColorArrayName = [None, '']
    plt00Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00Display.SelectOrientationVectors = 'None'
    plt00Display.ScaleFactor = 0.1
    plt00Display.SelectScaleArray = 'None'
    plt00Display.GlyphType = 'Arrow'
    plt00Display.GlyphTableIndexArray = 'None'
    plt00Display.GaussianRadius = 0.005
    plt00Display.SetScaleArray = [None, '']
    plt00Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00Display.OpacityArray = [None, '']
    plt00Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00Display.PolarAxes = 'PolarAxesRepresentation'
    plt00Display.ScalarOpacityUnitDistance = 0.0701538780193358

    # reset view to fit data
    renderView1.ResetCamera()

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

    # get the material library
    materialLibrary1 = GetMaterialLibrary()

    # update the view to ensure updated data information
    renderView1.Update()

    # change representation type
    plt00Display.SetRepresentationType('Surface')

    # set scalar coloring
    ColorBy(plt00Display, ('CELLS', 'phi'))

    # rescale color and/or opacity maps used to include current data range
    plt00Display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'phi'
    phiLUT = GetColorTransferFunction('phi')

    # get opacity transfer function/opacity map for 'phi'
    phiPWF = GetOpacityTransferFunction('phi')

    # show color bar/color legend
    plt00Display.SetScalarBarVisibility(renderView1, True)

    # create a new 'AMReX/BoxLib Grid Reader'
    plt00_1 = AMReXBoxLibGridReader(FileNames=AllPlotFiles)
    plt00_1.CellArrayStatus = []

    # show data in view
    plt00_1Display = Show(plt00_1, renderView1, 'AMRRepresentation')

    # trace defaults for the display properties.
    plt00_1Display.Representation = 'Outline'
    plt00_1Display.ColorArrayName = [None, '']
    plt00_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    plt00_1Display.SelectOrientationVectors = 'None'
    plt00_1Display.ScaleFactor = 0.1
    plt00_1Display.SelectScaleArray = 'None'
    plt00_1Display.GlyphType = 'Arrow'
    plt00_1Display.GlyphTableIndexArray = 'None'
    plt00_1Display.GaussianRadius = 0.005
    plt00_1Display.SetScaleArray = [None, '']
    plt00_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    plt00_1Display.OpacityArray = [None, '']
    plt00_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    plt00_1Display.DataAxesGrid = 'GridAxesRepresentation'
    plt00_1Display.PolarAxes = 'PolarAxesRepresentation'
    plt00_1Display.ScalarOpacityUnitDistance = 0.0701538780193358

    # update the view to ensure updated data information
    renderView1.Update()

    # change solid color
    plt00_1Display.AmbientColor = [0.0, 1.0, 0.0]
    plt00_1Display.DiffuseColor = [0.0, 1.0, 0.0]

    # get color legend/bar for phiLUT in view renderView1
    phiLUTColorBar = GetScalarBar(phiLUT, renderView1)

    # change scalar bar placement
    phiLUTColorBar.WindowLocation = 'AnyLocation'
    phiLUTColorBar.Position = [0, 0.75]
    phiLUTColorBar.ScalarBarLength = 0.2

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [0.5, 0.5, 10000.0]
    renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
    renderView1.CameraParallelScale = 0.5843857695756589

    # save animation
    output_movie_base = "amr101_2D"
    output_movie = output_movie_base + ".avi"
    SaveAnimation(output_movie,
                  renderView1,
                  ImageResolution=[1200, 1200],
                  FrameRate=args.frame_rate,
                  FrameWindow=[0, len(AllPlotFiles)-1])

    return output_movie_base, output_movie

def convert_avi_to_gif(output_movie_base, output_movie):
    # use ffmpeg to convert the avi movie into an animated gif
    ffmpeg_convert_to_gif = 'ffmpeg -y -i {} -vf "fps=35,scale={}:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 {}.gif'.format(output_movie, args.resolution, output_movie_base)
    subprocess.run(ffmpeg_convert_to_gif, shell=True)

if __name__ == "__main__":
    if not (args.spacedim == 2 or args.spacedim == 3):
        print("Please specify --spacedim D (with D=2 or D=3)")
        exit()

    if args.frame_rate <= 0:
        print("Please specify --frame_rate F (with F > 0)")
        exit()

    if args.resolution <= 0:
        print("Please specify --resolution R (with R > 0)")
        exit()

    # get all the plotfiles
    PlotFiles = sorted(glob.glob("plt" + "[0-9]"*5))

    # call the 2D or 3D vis function
    output_movie_base = None
    output_movie = None

    if args.spacedim == 3:
        output_movie_base, output_movie = generate_movie_3D(PlotFiles)
    elif args.spacedim == 2:
        output_movie_base, output_movie = generate_movie_2D(PlotFiles)

    # convert the avi movie into an animated gif
    convert_avi_to_gif(output_movie_base, output_movie)