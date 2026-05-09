#!/usr/bin/env python
import vtk
import chigger

camera = vtk.vtkCamera()
camera.SetViewUp(0.0000000000, 1.0000000000, 0.0000000000)
camera.SetPosition(0.6986057726, 0.4433031623, 2.7543779941)
camera.SetFocalPoint(0.6986057726, 0.4433031623, 0.0000000000)

reader = chigger.exodus.ExodusReader('inclined.e')
reader.setOptions(displacement_magnitude=20.0, block=['interLayer', 'metalLayer'])

result = chigger.exodus.ExodusResult(reader)
result.setOptions(edges=True, edge_color=[0, 0, 0], variable='conc', block=['interLayer', 'metalLayer'], local_range=True, camera=camera)

cbar = chigger.exodus.ExodusColorBar(result)
cbar.setOptions(colorbar_origin=(0.8, 0.25, 0.0))
cbar.setOptions('primary', lim=[0.0, 1.6107148046729415e-14], font_color=[0,0,0], font_size=16)

time = chigger.annotations.TextAnnotation(position=[0.5,0.9], font_size=24, text_color=[0,0,0], justification='center')

window = chigger.RenderWindow(result,cbar,time)
window.setOptions(size=None, style=None, background=[1, 1, 1])
#window.setOptions(size=[1920, 1080],style=None, background=[0.7058823529411765, 0.7058823529411765, 0.7058823529411765], background2=[0.43529411764705883, 0.43529411764705883, 0.43529411764705883], gradient_background=True)

window.start()
window.resetCameraClippingRange()
for i, t in enumerate(reader.getTimes()):
    reader.setOptions(timestep=i)
    time.setOptions(text='Time = {:5.4f} sec.'.format(t))
    filename = 'output/aa{:05d}.jpg'.format(i)
    window.write(filename)

chigger.utils.img2mov('output/aa*.jpg', 'result.webm', duration=10, num_threads=4)
