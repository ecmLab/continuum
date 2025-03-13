This folder contains scripts used to generate the animation demonstrating the 3D simulation results.

To generate the animation, follow these steps:
1. Download/generate the results (exodus and csv) here.
2. Run the `visualize_*.py` scripts to generate images for each of the field variable. Adjust parameters in these scripts as needed.
3. Run the `curves.py` script to generate C-V, C-t, V-t curves.
4. Run the `collage.py` script to collage all the results onto one big canvas. This will write a bunch of images in the `animation` subdirectory.
5. Use your favorite tool to convert the images in `animation` into a video. The command I used was

```bash
ffmpeg -framerate 2 -pattern_type glob -i 'animation/*.png' -c:v libx264 -pix_fmt yuv420p 3D_SSB.mov
```
