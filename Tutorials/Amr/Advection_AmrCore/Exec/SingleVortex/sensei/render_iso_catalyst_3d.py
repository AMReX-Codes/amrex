
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.4.1 64 bits

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=5
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.4.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1000, 700]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.5, 0.5, 3.2557533687070332]
      renderView1.CameraFocalPoint = [0.5, 0.5, -0.09031184624419736]
      renderView1.CameraParallelScale = 0.8660254037844386
      renderView1.Background = [0.0, 0.0, 0.0]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='pv_image_3d_%t.png', freq=1, fittoscreen=0, magnification=1, width=1000, height=700, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML UniformGrid AMR Reader'
      # create a producer from a simulation input
      mesh_000 = coprocessor.CreateProducer(datadescription, 'mesh')

      # create a new 'Cell Data to Point Data'
      cellDatatoPointData1 = CellDatatoPointData(Input=mesh_000)

      # create a new 'Contour'
      contour1 = Contour(Input=cellDatatoPointData1)
      contour1.ContourBy = ['POINTS', 'phi']
      contour1.ComputeScalars = 1
      contour1.Isosurfaces = [0.99429, 1.1043655555555556, 1.214441111111111, 1.3245166666666668, 1.4345922222222223, 1.5446677777777778, 1.6547433333333332, 1.764818888888889, 1.8748944444444444, 1.98497]
      contour1.PointMergeMethod = 'Uniform Binning'

      # create a new 'Annotate Time'
      annotateTime1 = AnnotateTime()
      annotateTime1.Format = 't = %0.2f'

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'phi'
      phiLUT = GetColorTransferFunction('phi')
      phiLUT.RGBPoints = [0.99429, 0.278431372549, 0.278431372549, 0.858823529412, 1.13595724, 0.0, 0.0, 0.360784313725, 1.2766338000000002, 0.0, 1.0, 1.0, 1.41929172, 0.0, 0.501960784314, 0.0, 1.55996828, 1.0, 1.0, 0.0, 1.70163552, 1.0, 0.380392156863, 0.0, 1.84330276, 0.419607843137, 0.0, 0.0, 1.9849700000000001, 0.878431372549, 0.301960784314, 0.301960784314]
      phiLUT.ColorSpace = 'RGB'
      phiLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'phi'
      phiPWF = GetOpacityTransferFunction('phi')
      phiPWF.Points = [0.99429, 0.0, 0.5, 0.0, 1.9849700000000001, 1.0, 0.5, 0.0]
      phiPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mesh_000
      mesh_000Display = Show(mesh_000, renderView1)
      # trace defaults for the display properties.
      mesh_000Display.Representation = 'AMR Blocks'
      mesh_000Display.ColorArrayName = [None, '']
      mesh_000Display.DiffuseColor = [0.0, 0.0, 0.0]
      mesh_000Display.OSPRayScaleArray = 'GhostType'
      mesh_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
      mesh_000Display.SelectOrientationVectors = 'None'
      mesh_000Display.ScaleFactor = 0.1
      mesh_000Display.SelectScaleArray = 'None'
      mesh_000Display.GlyphType = 'Arrow'
      mesh_000Display.GlyphTableIndexArray = 'None'
      mesh_000Display.DataAxesGrid = 'GridAxesRepresentation'
      mesh_000Display.PolarAxes = 'PolarAxesRepresentation'
      mesh_000Display.ScalarOpacityUnitDistance = 0.0174438098693218

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      mesh_000Display.DataAxesGrid.XTitle = 'X'
      mesh_000Display.DataAxesGrid.YTitle = 'Y'
      mesh_000Display.DataAxesGrid.ZTitle = 'Z'
      mesh_000Display.DataAxesGrid.XTitleBold = 1
      mesh_000Display.DataAxesGrid.XTitleFontSize = 14
      mesh_000Display.DataAxesGrid.YTitleBold = 1
      mesh_000Display.DataAxesGrid.YTitleFontSize = 14
      mesh_000Display.DataAxesGrid.ZTitleBold = 1
      mesh_000Display.DataAxesGrid.ZTitleFontSize = 14
      mesh_000Display.DataAxesGrid.XLabelBold = 1
      mesh_000Display.DataAxesGrid.XLabelFontSize = 14
      mesh_000Display.DataAxesGrid.YLabelBold = 1
      mesh_000Display.DataAxesGrid.YLabelFontSize = 14
      mesh_000Display.DataAxesGrid.ZLabelBold = 1
      mesh_000Display.DataAxesGrid.ZLabelFontSize = 14

      # show data from contour1
      contour1Display = Show(contour1, renderView1)
      # trace defaults for the display properties.
      contour1Display.Representation = 'Surface'
      contour1Display.ColorArrayName = ['POINTS', 'phi']
      contour1Display.LookupTable = phiLUT
      contour1Display.OSPRayScaleArray = 'GhostType'
      contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      contour1Display.SelectOrientationVectors = 'GhostType'
      contour1Display.ScaleFactor = 0.0572519063949585
      contour1Display.SelectScaleArray = 'GhostType'
      contour1Display.GlyphType = 'Arrow'
      contour1Display.GlyphTableIndexArray = 'GhostType'
      contour1Display.DataAxesGrid = 'GridAxesRepresentation'
      contour1Display.PolarAxes = 'PolarAxesRepresentation'
      contour1Display.GaussianRadius = 0.02862595319747925
      contour1Display.SetScaleArray = ['POINTS', 'GhostType']
      contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
      contour1Display.OpacityArray = ['POINTS', 'GhostType']
      contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      contour1Display.SetScalarBarVisibility(renderView1, True)

      # show data from annotateTime1
      annotateTime1Display = Show(annotateTime1, renderView1)
      # trace defaults for the display properties.
      annotateTime1Display.Bold = 1
      annotateTime1Display.FontSize = 12
      annotateTime1Display.WindowLocation = 'LowerLeftCorner'

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for phiLUT in view renderView1
      phiLUTColorBar = GetScalarBar(phiLUT, renderView1)
      phiLUTColorBar.WindowLocation = 'AnyLocation'
      phiLUTColorBar.Position = [0.852, 0.07857142857142851]
      phiLUTColorBar.Title = 'phi'
      phiLUTColorBar.ComponentTitle = ''
      phiLUTColorBar.TitleBold = 1
      phiLUTColorBar.TitleFontSize = 24
      phiLUTColorBar.LabelBold = 1
      phiLUTColorBar.LabelFontSize = 18
      phiLUTColorBar.ScalarBarThickness = 24
      phiLUTColorBar.ScalarBarLength = 0.8357142857142857

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(mesh_000)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'mesh': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
