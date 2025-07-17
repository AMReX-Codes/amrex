# state file generated using paraview version 5.13.0
import paraview

import glob

paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView("RenderView")
renderView1.ViewSize = [1093, 931]
renderView1.AxesGrid = "Grid Axes 3D Actor"
renderView1.StereoType = "Crystal Eyes"
renderView1.CameraPosition = [4.0871501253875975, 2.2332036879247603, 8.620196826886515]
renderView1.CameraViewUp = [
    0.28405525969665424,
    0.886848733903073,
    -0.36443371497870924,
]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.7320508075688772
renderView1.LegendGrid = "Legend Grid Actor"
renderView1.PolarGrid = "Polar Grid Actor"
renderView1.BackEnd = "OSPRay raycaster"
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name="Layout #1")
layout1.AssignView(0, renderView1)
layout1.SetSize(1093, 931)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# User edit: create a glob which returns a list of matching path names
# in cwd
PlotFiles = sorted(glob.glob("plt" + "[0-9]" * 5))

# create a new 'AMReX/BoxLib Grid Reader'
plt00000 = AMReXBoxLibGridReader(registrationName="plt00000*", FileNames=PlotFiles)

plt00000.CellArrayStatus = ["phi"]

# create a new 'Slice'
slice1 = Slice(registrationName="Slice1", Input=plt00000)
slice1.SliceType = "Plane"
slice1.HyperTreeGridSlicer = "Plane"
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = "Uniform Binning"

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [
    0.25835322071340727,
    0.11812878588415326,
    0.13969875959104724,
]
slice1.SliceType.Normal = [0.924710841412284, -0.027237427802698565, 0.3796945908243607]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, "GeometryRepresentation")

# get 2D transfer function for 'phi'
phiTF2D = GetTransferFunction2D("phi")

# get color transfer function/color map for 'phi'
phiLUT = GetColorTransferFunction("phi")
phiLUT.TransferFunction2D = phiTF2D
phiLUT.RGBPoints = [
    1.0,
    0.231373,
    0.298039,
    0.752941,
    1.4909280363367723,
    0.865003,
    0.865003,
    0.865003,
    1.9818560726735446,
    0.705882,
    0.0156863,
    0.14902,
]
phiLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = "Surface"
slice1Display.ColorArrayName = ["CELLS", "phi"]
slice1Display.LookupTable = phiLUT
slice1Display.SelectNormalArray = "None"
slice1Display.SelectTangentArray = "None"
slice1Display.SelectTCoordArray = "None"
slice1Display.TextureTransform = "Transform2"
slice1Display.OSPRayScaleFunction = "Piecewise Function"
slice1Display.Assembly = "Hierarchy"
slice1Display.SelectedBlockSelectors = [""]
slice1Display.SelectOrientationVectors = "None"
slice1Display.ScaleFactor = 0.2
slice1Display.SelectScaleArray = "None"
slice1Display.GlyphType = "Arrow"
slice1Display.GlyphTableIndexArray = "None"
slice1Display.GaussianRadius = 0.01
slice1Display.SetScaleArray = [None, ""]
slice1Display.ScaleTransferFunction = "Piecewise Function"
slice1Display.OpacityArray = [None, ""]
slice1Display.OpacityTransferFunction = "Piecewise Function"
slice1Display.DataAxesGrid = "Grid Axes Representation"
slice1Display.PolarAxes = "Polar Axes Representation"
slice1Display.SelectInputVectors = [None, ""]
slice1Display.WriteLog = ""

# setup the color legend parameters for each legend in this view

# get color legend/bar for phiLUT in view renderView1
phiLUTColorBar = GetScalarBar(phiLUT, renderView1)
phiLUTColorBar.Title = "phi"
phiLUTColorBar.ComponentTitle = ""

# set color bar visibility
phiLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'phi'
phiPWF = GetOpacityTransferFunction("phi")
phiPWF.Points = [1.0, 0.0, 0.5, 0.0, 1.9818560726735446, 1.0, 0.5, 0.0]
phiPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 0.0
animationScene1.EndTime = 10.0
animationScene1.PlayMode = "Snap To TimeSteps"

# ----------------------------------------------------------------
# restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, useful when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple split view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
