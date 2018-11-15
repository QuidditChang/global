import sys
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

timestep=int(sys.argv[1])
age=int(sys.argv[2])

# create a new 'Legacy VTK Reader'
global009580vtk = LegacyVTKReader(FileNames=['global.00.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global019580vtk = LegacyVTKReader(FileNames=['global.01.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global029580vtk = LegacyVTKReader(FileNames=['global.02.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global039580vtk = LegacyVTKReader(FileNames=['global.03.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global049580vtk = LegacyVTKReader(FileNames=['global.04.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global059580vtk = LegacyVTKReader(FileNames=['global.05.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global069580vtk = LegacyVTKReader(FileNames=['global.06.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global079580vtk = LegacyVTKReader(FileNames=['global.07.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global089580vtk = LegacyVTKReader(FileNames=['global.08.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global099580vtk = LegacyVTKReader(FileNames=['global.09.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global109580vtk = LegacyVTKReader(FileNames=['global.10.%d.vtk' % timestep])

# create a new 'Legacy VTK Reader'
global119580vtk = LegacyVTKReader(FileNames=['global.11.%d.vtk' % timestep])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1067, 582]

# get color transfer function/color map for 'Viscosity'
viscosityLUT = GetColorTransferFunction('Viscosity')

# show data in view
global009580vtkDisplay = Show(global009580vtk, renderView1)
# trace defaults for the display properties.
global009580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global009580vtkDisplay.LookupTable = viscosityLUT
#global009580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global009580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global009580vtkDisplay.GlyphType = 'Arrow'
#global009580vtkDisplay.ScalarOpacityUnitDistance = 137563.48109841836
#global009580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global009580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global009580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global009580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
global009580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global039580vtkDisplay = Show(global039580vtk, renderView1)
# trace defaults for the display properties.
global039580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global039580vtkDisplay.LookupTable = viscosityLUT
#global039580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global039580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global039580vtkDisplay.GlyphType = 'Arrow'
#global039580vtkDisplay.ScalarOpacityUnitDistance = 137506.69732141742
#global039580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global039580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global039580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global039580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global039580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global049580vtkDisplay = Show(global049580vtk, renderView1)
# trace defaults for the display properties.
global049580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global049580vtkDisplay.LookupTable = viscosityLUT
#global049580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global049580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global049580vtkDisplay.GlyphType = 'Arrow'
#global049580vtkDisplay.ScalarOpacityUnitDistance = 167117.74680552137
#global049580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global049580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global049580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global049580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global049580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global079580vtkDisplay = Show(global079580vtk, renderView1)
# trace defaults for the display properties.
global079580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global079580vtkDisplay.LookupTable = viscosityLUT
#global079580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global079580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global079580vtkDisplay.GlyphType = 'Arrow'
#global079580vtkDisplay.ScalarOpacityUnitDistance = 167117.54590736775
#global079580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global079580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global079580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global079580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global079580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global089580vtkDisplay = Show(global089580vtk, renderView1)
# trace defaults for the display properties.
global089580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global089580vtkDisplay.LookupTable = viscosityLUT
#global089580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global089580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global089580vtkDisplay.GlyphType = 'Arrow'
#global089580vtkDisplay.ScalarOpacityUnitDistance = 137506.74415500654
#global089580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global089580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global089580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global089580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global089580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global099580vtkDisplay = Show(global099580vtk, renderView1)
# trace defaults for the display properties.
global099580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global099580vtkDisplay.LookupTable = viscosityLUT
#global099580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global099580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global099580vtkDisplay.GlyphType = 'Arrow'
#global099580vtkDisplay.ScalarOpacityUnitDistance = 137563.46585171248
#global099580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global099580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global099580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global099580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global099580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global119580vtkDisplay = Show(global119580vtk, renderView1)
# trace defaults for the display properties.
global119580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global119580vtkDisplay.LookupTable = viscosityLUT
#global119580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global119580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global119580vtkDisplay.GlyphType = 'Arrow'
#global119580vtkDisplay.ScalarOpacityUnitDistance = 137506.74396808128
#global119580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global119580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global119580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global119580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global119580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global019580vtkDisplay = Show(global019580vtk, renderView1)
# trace defaults for the display properties.
global019580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global019580vtkDisplay.LookupTable = viscosityLUT
#global019580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global019580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global019580vtkDisplay.GlyphType = 'Arrow'
#global019580vtkDisplay.ScalarOpacityUnitDistance = 167117.54590736775
#global019580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global019580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global019580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global019580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global019580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global109580vtkDisplay = Show(global109580vtk, renderView1)
# trace defaults for the display properties.
global109580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global109580vtkDisplay.LookupTable = viscosityLUT
#global109580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global109580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global109580vtkDisplay.GlyphType = 'Arrow'
#global109580vtkDisplay.ScalarOpacityUnitDistance = 167117.74680552137
#global109580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global109580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global109580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global109580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global109580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global029580vtkDisplay = Show(global029580vtk, renderView1)
# trace defaults for the display properties.
global029580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global029580vtkDisplay.LookupTable = viscosityLUT
#global029580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global029580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global029580vtkDisplay.GlyphType = 'Arrow'
#global029580vtkDisplay.ScalarOpacityUnitDistance = 137563.39313921705
#global029580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global029580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global029580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global029580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global029580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global059580vtkDisplay = Show(global059580vtk, renderView1)
# trace defaults for the display properties.
global059580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global059580vtkDisplay.LookupTable = viscosityLUT
#global059580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global059580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global059580vtkDisplay.GlyphType = 'Arrow'
#global059580vtkDisplay.ScalarOpacityUnitDistance = 137563.4475703484
#global059580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global059580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global059580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global059580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global059580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# show data in view
global069580vtkDisplay = Show(global069580vtk, renderView1)
# trace defaults for the display properties.
global069580vtkDisplay.ColorArrayName = ['POINTS', 'Viscosity']
global069580vtkDisplay.LookupTable = viscosityLUT
#global069580vtkDisplay.OSPRayScaleArray = 'Viscosity'
#global069580vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
#global069580vtkDisplay.GlyphType = 'Arrow'
#global069580vtkDisplay.ScalarOpacityUnitDistance = 137506.75862178658
#global069580vtkDisplay.SetScaleArray = ['POINTS', 'Viscosity']
#global069580vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
#global069580vtkDisplay.OpacityArray = ['POINTS', 'Viscosity']
#global069580vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'

# show color bar/color legend
global069580vtkDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Viscosity'
viscosityPWF = GetOpacityTransferFunction('Viscosity')

# set active source
SetActiveSource(global009580vtk)

# set active source
SetActiveSource(global119580vtk)

# set active source
SetActiveSource(global009580vtk)

# set active source
SetActiveSource(global119580vtk)

# create a new 'Append Datasets'
appendDatasets1 = AppendDatasets(Input=[global069580vtk, global059580vtk, global019580vtk, global109580vtk, global029580vtk, global009580vtk, global049580vtk, global039580vtk, global079580vtk, global089580vtk, global099580vtk, global119580vtk])

# show data in view
appendDatasets1Display = Show(appendDatasets1, renderView1)
# trace defaults for the display properties.
appendDatasets1Display.ColorArrayName = ['POINTS', 'Viscosity']
appendDatasets1Display.LookupTable = viscosityLUT
#appendDatasets1Display.OSPRayScaleArray = 'Viscosity'
#appendDatasets1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#appendDatasets1Display.GlyphType = 'Arrow'
#appendDatasets1Display.ScalarOpacityUnitDistance = 131581.13592878426
#appendDatasets1Display.SetScaleArray = ['POINTS', 'Viscosity']
#appendDatasets1Display.ScaleTransferFunction = 'PiecewiseFunction'
#appendDatasets1Display.OpacityArray = ['POINTS', 'Viscosity']
#appendDatasets1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(global109580vtk, renderView1)

# hide data in view
Hide(global079580vtk, renderView1)

# hide data in view
Hide(global099580vtk, renderView1)

# hide data in view
Hide(global029580vtk, renderView1)

# hide data in view
Hide(global059580vtk, renderView1)

# hide data in view
Hide(global119580vtk, renderView1)

# hide data in view
Hide(global009580vtk, renderView1)

# hide data in view
Hide(global069580vtk, renderView1)

# hide data in view
Hide(global049580vtk, renderView1)

# hide data in view
Hide(global019580vtk, renderView1)

# hide data in view
Hide(global089580vtk, renderView1)

# hide data in view
Hide(global039580vtk, renderView1)

# show color bar/color legend
appendDatasets1Display.SetScalarBarVisibility(renderView1, True)

# invert the transfer function
viscosityLUT.InvertTransferFunction()

# convert to log space
viscosityLUT.MapControlPointsToLogSpace()

# Properties modified on viscosityLUT
viscosityLUT.UseLogScale = 1

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(Input=appendDatasets1)
isoVolume1.InputScalars = ['POINTS', 'Viscosity']
isoVolume1.ThresholdRange = [0.009999999776482582, 100.0]

# Properties modified on isoVolume1
isoVolume1.InputScalars = ['POINTS', 'Temperature']
isoVolume1.ThresholdRange = [-0.012071309611201286, 0.5]

# show data in view
isoVolume1Display = Show(isoVolume1, renderView1)
# trace defaults for the display properties.
isoVolume1Display.ColorArrayName = ['POINTS', 'Viscosity']
isoVolume1Display.LookupTable = viscosityLUT
#isoVolume1Display.OSPRayScaleArray = 'Viscosity'
#isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#isoVolume1Display.GlyphType = 'Arrow'
#isoVolume1Display.ScalarOpacityUnitDistance = 172290.11255508967
#isoVolume1Display.SetScaleArray = ['POINTS', 'Viscosity']
#isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
#isoVolume1Display.OpacityArray = ['POINTS', 'Viscosity']
#isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(appendDatasets1, renderView1)

# show color bar/color legend
isoVolume1Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)

# create a new 'Threshold'
threshold1 = Threshold(Input=isoVolume1)
threshold1.Scalars = ['POINTS', 'Viscosity']
threshold1.ThresholdRange = [0.009999999776482582, 100.0]

# Properties modified on threshold1
threshold1.Scalars = ['POINTS', 'Depth']
threshold1.ThresholdRange = [200.0, 2765.010009765625]

# show data in view
threshold1Display = Show(threshold1, renderView1)
# trace defaults for the display properties.
threshold1Display.ColorArrayName = ['POINTS', 'Viscosity']
threshold1Display.LookupTable = viscosityLUT
#threshold1Display.OSPRayScaleArray = 'Viscosity'
#threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#threshold1Display.GlyphType = 'Arrow'
#threshold1Display.ScalarOpacityUnitDistance = 189824.67269906844
#threshold1Display.SetScaleArray = ['POINTS', 'Viscosity']
#threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
#threshold1Display.OpacityArray = ['POINTS', 'Viscosity']
#threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(isoVolume1, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)

# set scalar coloring
ColorBy(threshold1Display, ('POINTS', 'Depth'))

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Depth'
depthLUT = GetColorTransferFunction('Depth')

# get opacity transfer function/opacity map for 'Depth'
depthPWF = GetOpacityTransferFunction('Depth')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
depthLUT.ApplyPreset('Warm to Cool (Extended)', True)

# invert the transfer function
depthLUT.InvertTransferFunction()

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
depthLUT.ApplyPreset('Cool to Warm (Extended)', True)

# set active source
SetActiveSource(appendDatasets1)

# create a new 'Iso Volume'
isoVolume2 = IsoVolume(Input=appendDatasets1)
isoVolume2.InputScalars = ['POINTS', 'Viscosity']
isoVolume2.ThresholdRange = [0.009999999776482582, 100.0]

# Properties modified on isoVolume2
isoVolume2.InputScalars = ['POINTS', 'Temperature']
isoVolume2.ThresholdRange = [0.8, 1.2109349966049194]

# show data in view
isoVolume2Display = Show(isoVolume2, renderView1)
# trace defaults for the display properties.
isoVolume2Display.ColorArrayName = ['POINTS', 'Viscosity']
isoVolume2Display.LookupTable = viscosityLUT
#isoVolume2Display.OSPRayScaleArray = 'Viscosity'
#isoVolume2Display.OSPRayScaleFunction = 'PiecewiseFunction'
#isoVolume2Display.GlyphType = 'Arrow'
#isoVolume2Display.ScalarOpacityUnitDistance = 183220.8825024971
#isoVolume2Display.SetScaleArray = ['POINTS', 'Viscosity']
#isoVolume2Display.ScaleTransferFunction = 'PiecewiseFunction'
#isoVolume2Display.OpacityArray = ['POINTS', 'Viscosity']
#isoVolume2Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(appendDatasets1, renderView1)

# show color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, True)

# set scalar coloring
ColorBy(isoVolume2Display, ('POINTS', 'Temperature'))

# rescale color and/or opacity maps used to include current data range
isoVolume2Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Temperature'
temperatureLUT = GetColorTransferFunction('Temperature')

# get opacity transfer function/opacity map for 'Temperature'
temperaturePWF = GetOpacityTransferFunction('Temperature')

# rescale color and/or opacity maps used to exactly fit the current data range
isoVolume2Display.RescaleTransferFunctionToDataRange(False)

# Rescale transfer function
temperatureLUT.RescaleTransferFunction(0.0, 0.7)

# Rescale transfer function
temperaturePWF.RescaleTransferFunction(0.0, 0.7)

# hide color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, False)

# show color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(appendDatasets1)

# create a new 'Contour'
contour1 = Contour(Input=appendDatasets1)
contour1.ContourBy = ['POINTS', 'Viscosity']
contour1.Isosurfaces = [50.00499999988824]
contour1.PointMergeMethod = 'Uniform Binning'

# Properties modified on contour1
contour1.ContourBy = ['POINTS', 'Depth']
contour1.Isosurfaces = [50.0]

# show data in view
contour1Display = Show(contour1, renderView1)
# trace defaults for the display properties.
contour1Display.ColorArrayName = ['POINTS', 'Viscosity']
contour1Display.LookupTable = viscosityLUT
#contour1Display.OSPRayScaleArray = 'Colat'
#contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
#contour1Display.GlyphType = 'Arrow'
#contour1Display.SetScaleArray = ['POINTS', 'Colat']
#contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
#contour1Display.OpacityArray = ['POINTS', 'Colat']
#contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

# hide data in view
Hide(appendDatasets1, renderView1)

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)

# Properties modified on contour1Display
contour1Display.Opacity = 0.18

# Properties modified on renderView1
renderView1.LightSwitch = 1

# Properties modified on renderView1
renderView1.LightIntensity = 0.5000000000000001

# set active source
SetActiveSource(None)

# set active source
SetActiveSource(appendDatasets1)

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=appendDatasets1)

# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Time: %f Ma'
annotateTimeFilter1.Shift = age

# show data in view
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)

# Properties modified on annotateTimeFilter1
annotateTimeFilter1.Format = 'Time: %.0f Ma'

# Rescale transfer function
viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)

# set active source
SetActiveSource(None)

if timestep==9580:
    # create a new 'XML PolyData Reader'
    coastlines_Los_Alamosvtp = XMLPolyDataReader(FileName=['Coastlines_Los_Alamos.vtp'])
    
    # show data in view
    coastlines_Los_AlamosvtpDisplay = Show(coastlines_Los_Alamosvtp, renderView1)
    # trace defaults for the display properties.
    #coastlines_Los_AlamosvtpDisplay.ColorArrayName = [None, '']
    #coastlines_Los_AlamosvtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    #coastlines_Los_AlamosvtpDisplay.GlyphType = 'Arrow'
    #coastlines_Los_AlamosvtpDisplay.SetScaleArray = [None, '']
    #coastlines_Los_AlamosvtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    #coastlines_Los_AlamosvtpDisplay.OpacityArray = [None, '']
    #coastlines_Los_AlamosvtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    
    # Rescale transfer function
    viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)
    
    # create a new 'Calculator'
    calculator1 = Calculator(Input=coastlines_Los_Alamosvtp)
    calculator1.Function = ''
    
    # Properties modified on calculator1
    calculator1.CoordinateResults = 1
    calculator1.Function = 'coordsX*6371000.0*iHat+coordsY*6371000.0*jHat+coordsZ*6371000.0*kHat'
    
    # show data in view
    calculator1Display = Show(calculator1, renderView1)
    # trace defaults for the display properties.
    #calculator1Display.ColorArrayName = [None, '']
    #calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    #calculator1Display.GlyphType = 'Arrow'
    #calculator1Display.SetScaleArray = [None, '']
    #calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
    #calculator1Display.OpacityArray = [None, '']
    #calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'
    
    # hide data in view
    Hide(coastlines_Los_Alamosvtp, renderView1)
    
    # Rescale transfer function
    viscosityLUT.RescaleTransferFunction(0.00999999977648, 100.0)
    
    # set active source
    SetActiveSource(coastlines_Los_Alamosvtp)
    
    # set active source
    SetActiveSource(calculator1)
    
    # Properties modified on calculator1Display
    calculator1Display.Opacity = 0.58

# export view
ExportView('global.%dMa.pdf' % age, view=renderView1)

# current camera placement for renderView1
renderView1.CameraPosition = [12206953.947447909, -25659427.67156836, 6371443.83028981]
renderView1.CameraViewUp = [-0.19884720890143004, 0.14603511793722318, 0.9690890216286491]
renderView1.CameraParallelScale = 11034884.148018954

# save screenshot
SaveScreenshot('global.%dMa.png' % age, magnification=2, quality=100, view=renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [12206953.947447909, -25659427.67156836, 6371443.83028981]
renderView1.CameraViewUp = [-0.19884720890143004, 0.14603511793722318, 0.9690890216286491]
renderView1.CameraParallelScale = 11034884.148018954

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
