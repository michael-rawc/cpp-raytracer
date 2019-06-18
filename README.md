# cpp-raytracer
Raytrace based image renderer in C++

Basic, non-professionally optimised image renderer.
Hardcoded values in most places.
Will produce a png image with multiple differently coloured spheres and one loaded .obj mesh.

## To compile:
First, navigate to cpp-raytracer main file.
Then, execute: g++ src\renderImage.cpp src\lodepng.cpp -o renderImage
Then, run: renderImage
  
Output image will be out.png

## Changing parameters

To change which object to use in the scene, set the constant objectChoice to either "CubeTop.obj", "lowpolybunny.obj", or "bunny.obj"

To change the Anti-Aliasing value, set the constant SAMPLES to any square number (9 will give 3x3 AA), 1 gives no AA.

To change the reflection/freaction depth value, set the constant DEPTH to any number, default is 2, lower is faster.

To change the resolution of the image (image window scales with resolution), set the consants PIXELHEIGHT and PIXELWEIGHT appropriately.

-> All constants can be found in lines 17-24 of renderImage.cpp.

## Credits to:

Robert. S for OBJ-Loader, found at: https://github.com/Bly7/OBJ-Loader
Used to get information from .obj files.

Lode Vandevenne for lodepng, found at: https://github.com/lvandeve/lodepng
Used to form a .png image for output.
