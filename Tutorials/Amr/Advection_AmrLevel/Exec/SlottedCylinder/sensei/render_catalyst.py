
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# ParaView 5.4.1 64 bits

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=4
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
      renderView1.ViewSize = [480, 480]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.5, 0.5, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.5, 0.5, 1.9306242333430546]
      renderView1.CameraFocalPoint = [0.5, 0.5, -0.8014265742258199]
      renderView1.CameraParallelScale = 0.7071067811865476
      renderView1.Background = [0.0, 0.0, 0.0]

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=480, height=480, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML UniformGrid AMR Reader'
      # create a producer from a simulation input
      amr_mesh_ = coprocessor.CreateProducer(datadescription, 'mesh')

      # create a new 'Outline'
      outline1 = Outline(Input=amr_mesh_)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get color transfer function/color map for 'phi'
      phiLUT = GetColorTransferFunction('phi')
      phiLUT.RGBPoints = [0.9990617484291469, 0.278431372549, 0.278431372549, 0.858823529412, 1.1420481577953918, 0.0, 0.0, 0.360784313725, 1.2840346622010337, 0.0, 1.0, 1.0, 1.4280209765278817, 0.0, 0.501960784314, 0.0, 1.5700074809335236, 1.0, 1.0, 0.0, 1.7129938902997686, 1.0, 0.380392156863, 0.0, 1.8559802996660135, 0.419607843137, 0.0, 0.0, 1.9989667090322585, 0.878431372549, 0.301960784314, 0.301960784314]
      phiLUT.ColorSpace = 'RGB'
      phiLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'phi'
      phiPWF = GetOpacityTransferFunction('phi')
      phiPWF.Points = [0.9990617484291469, 0.0, 0.5, 0.0, 1.9989667090322585, 1.0, 0.5, 0.0]
      phiPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from amr_mesh_
      amr_mesh_Display = Show(amr_mesh_, renderView1)
      # trace defaults for the display properties.
      amr_mesh_Display.Representation = 'Wireframe'
      amr_mesh_Display.ColorArrayName = ['CELLS', 'phi']
      amr_mesh_Display.LookupTable = phiLUT
      amr_mesh_Display.OSPRayScaleArray = 'phi'
      amr_mesh_Display.OSPRayScaleFunction = 'PiecewiseFunction'
      amr_mesh_Display.SelectOrientationVectors = 'None'
      amr_mesh_Display.ScaleFactor = 0.1
      amr_mesh_Display.SelectScaleArray = 'None'
      amr_mesh_Display.GlyphType = 'Arrow'
      amr_mesh_Display.GlyphTableIndexArray = 'None'
      amr_mesh_Display.DataAxesGrid = 'GridAxesRepresentation'
      amr_mesh_Display.PolarAxes = 'PolarAxesRepresentation'
      amr_mesh_Display.ScalarOpacityUnitDistance = 0.08838834764831846
      amr_mesh_Display.ScalarOpacityFunction = phiPWF

      # show data from outline1
      outline1Display = Show(outline1, renderView1)
      # trace defaults for the display properties.
      outline1Display.Representation = 'Surface'
      outline1Display.ColorArrayName = [None, '']
      outline1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      outline1Display.SelectOrientationVectors = 'None'
      outline1Display.ScaleFactor = 0.1
      outline1Display.SelectScaleArray = 'None'
      outline1Display.GlyphType = 'Arrow'
      outline1Display.GlyphTableIndexArray = 'None'
      outline1Display.DataAxesGrid = 'GridAxesRepresentation'
      outline1Display.PolarAxes = 'PolarAxesRepresentation'
      outline1Display.GaussianRadius = 0.05
      outline1Display.SetScaleArray = [None, '']
      outline1Display.ScaleTransferFunction = 'PiecewiseFunction'
      outline1Display.OpacityArray = [None, '']
      outline1Display.OpacityTransferFunction = 'PiecewiseFunction'

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(outline1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'mesh': [1, 1]}
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
