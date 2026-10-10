"""Render the EB surface written by WriteEBSurface (eb.pvtp) with ParaView.

usage: pvbatch plot_eb.py <run_dir> <output.png> [title]
where run_dir has the output of a run with case=N write_surface=1.
"""
import sys

from paraview.simple import *  # noqa: F401,F403

run_dir, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else ""

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [2000, 1260]
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 1
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass

eb = OpenDataFile(run_dir + "/eb.pvtp")
d = Show(eb, view)
d.SetRepresentationType("Surface With Edges")
d.AmbientColor = d.DiffuseColor = [0.2, 0.2, 0.75]
d.EdgeColor = [0.0, 0.0, 0.35]
d.ColorArrayName = [None, ""]

if title:
    t = Text(Text=title)
    td = Show(t, view)
    td.Color = [0, 0, 0]
    td.FontSize = 28
    td.WindowLocation = "Upper Center"

cam = view.GetActiveCamera()
cam.SetViewUp(0.0, 0.0, 1.0)
cam.SetPosition(-0.6, -1.0, 0.8)
cam.SetFocalPoint(0.0, 0.0, 0.0)
view.ResetCamera(False)
Render(view)
SaveScreenshot(out, view, ImageResolution=[2000, 1260])
