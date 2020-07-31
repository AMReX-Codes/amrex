
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
      renderView1.ViewSize = [1000, 800]
      renderView1.InteractionMode = '2D'
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.OrientationAxesOutlineColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.6499329008559785, 0.4879518204669307, 10000.0]
      renderView1.CameraFocalPoint = [0.6499329008559785, 0.4879518204669307, 0.0]
      renderView1.CameraParallelScale = 0.5952751765669945
      renderView1.Background = [1.0, 1.0, 1.0]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.Visibility = 1
      renderView1.AxesGrid.XTitle = 'X'
      renderView1.AxesGrid.YTitle = 'Y '
      renderView1.AxesGrid.ZTitle = ''
      renderView1.AxesGrid.XTitleColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.AxesGrid.XTitleBold = 1
      renderView1.AxesGrid.XTitleFontSize = 16
      renderView1.AxesGrid.YTitleColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.AxesGrid.YTitleBold = 1
      renderView1.AxesGrid.YTitleFontSize = 16
      renderView1.AxesGrid.GridColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.AxesGrid.XLabelColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.AxesGrid.XLabelBold = 1
      renderView1.AxesGrid.YLabelColor = [0.0, 0.0, 0.4980392156862745]
      renderView1.AxesGrid.YLabelBold = 1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='pv_image_2d_%t.png', freq=1, fittoscreen=0, magnification=1, width=1000, height=800, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML UniformGrid AMR Reader'
      # create a producer from a simulation input
      mesh_0000 = coprocessor.CreateProducer(datadescription, 'mesh')

      # create a new 'Cell Data to Point Data'
      cellDatatoPointData1 = CellDatatoPointData(Input=mesh_0000)

      # create a new 'Contour'
      contour1 = Contour(Input=cellDatatoPointData1)
      contour1.ContourBy = ['POINTS', 'phi']
      contour1.ComputeScalars = 1
      contour1.Isosurfaces = [0.99517, 1.1054688888888888, 1.2157677777777778, 1.3260666666666667, 1.4363655555555555, 1.5466644444444444, 1.6569633333333333, 1.767262222222222, 1.877561111111111, 1.98786]
      contour1.PointMergeMethod = 'Uniform Binning'

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'phi'
      phiLUT = GetColorTransferFunction('phi')
      phiLUT.RGBPoints = [1.0, 0.278431372549, 0.278431372549, 0.858823529412, 1.1428909412660986, 0.0, 0.0, 0.360784313725, 1.2847826451806859, 0.0, 1.0, 1.0, 1.4286728237982957, 0.0, 0.501960784314, 0.0, 1.570564527712883, 1.0, 1.0, 0.0, 1.7134554689789816, 1.0, 0.380392156863, 0.0, 1.8563464102450802, 0.419607843137, 0.0, 0.0, 1.9992373515111788, 0.878431372549, 0.301960784314, 0.301960784314]
      phiLUT.ColorSpace = 'RGB'
      phiLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'phi'
      phiPWF = GetOpacityTransferFunction('phi')
      phiPWF.Points = [1.0, 0.0, 0.5, 0.0, 1.9992373515111788, 1.0, 0.5, 0.0]
      phiPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from mesh_0000
      mesh_0000Display = Show(mesh_0000, renderView1)
      # trace defaults for the display properties.
      mesh_0000Display.Representation = 'Wireframe'
      mesh_0000Display.AmbientColor = [0.0, 0.0, 0.0]
      mesh_0000Display.ColorArrayName = ['POINTS', '']
      mesh_0000Display.OSPRayScaleArray = 'GhostType'
      mesh_0000Display.OSPRayScaleFunction = 'PiecewiseFunction'
      mesh_0000Display.SelectOrientationVectors = 'None'
      mesh_0000Display.ScaleFactor = 0.1
      mesh_0000Display.SelectScaleArray = 'None'
      mesh_0000Display.GlyphType = 'Arrow'
      mesh_0000Display.GlyphTableIndexArray = 'None'
      mesh_0000Display.DataAxesGrid = 'GridAxesRepresentation'
      mesh_0000Display.PolarAxes = 'PolarAxesRepresentation'
      mesh_0000Display.ScalarOpacityUnitDistance = 0.057873097067582834

      # show data from contour1
      contour1Display = Show(contour1, renderView1)
      # trace defaults for the display properties.
      contour1Display.Representation = 'Surface'
      contour1Display.ColorArrayName = ['POINTS', 'phi']
      contour1Display.LookupTable = phiLUT
      contour1Display.LineWidth = 2.0
      contour1Display.OSPRayScaleArray = 'phi'
      contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      contour1Display.SelectOrientationVectors = 'None'
      contour1Display.ScaleFactor = 0.05746527910232544
      contour1Display.SelectScaleArray = 'None'
      contour1Display.GlyphType = 'Arrow'
      contour1Display.GlyphTableIndexArray = 'None'
      contour1Display.DataAxesGrid = 'GridAxesRepresentation'
      contour1Display.PolarAxes = 'PolarAxesRepresentation'
      contour1Display.GaussianRadius = 0.02873263955116272
      contour1Display.SetScaleArray = ['POINTS', 'phi']
      contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
      contour1Display.OpacityArray = ['POINTS', 'phi']
      contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # show color legend
      contour1Display.SetScalarBarVisibility(renderView1, True)

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for phiLUT in view renderView1
      phiLUTColorBar = GetScalarBar(phiLUT, renderView1)
      phiLUTColorBar.WindowLocation = 'AnyLocation'
      phiLUTColorBar.Position = [0.7951562900333077, 0.0840151515151517]
      phiLUTColorBar.Title = 'phi'
      phiLUTColorBar.ComponentTitle = ''
      phiLUTColorBar.TitleColor = [0.0, 0.0, 0.0]
      phiLUTColorBar.TitleBold = 1
      phiLUTColorBar.TitleFontSize = 32
      phiLUTColorBar.LabelColor = [0.0, 0.0, 0.0]
      phiLUTColorBar.LabelBold = 1
      phiLUTColorBar.LabelFontSize = 24
      phiLUTColorBar.ScalarBarThickness = 28
      phiLUTColorBar.ScalarBarLength = 0.848928571428571

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(contour1)
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
