# cpp-raytracer
Raytrace based image renderer in C++

Basic, non-professionally optimised image renderer.
Hardcoded values in most places.

To compile:
  First, navigate to cpp-raytracer main file.
  Then, execute: g++ src\renderImage.cpp src\lodepng.cpp -o renderImage
  Then, run: renderImage
  
Output image will be in out.png
  
To change which object to use, set the constant objectChoice to either "CubeTop.obj", "lowpolybunny.obj", or "bunny.obj"\n
To change the AA value, set the constant SAMPLES to any square number, 1 gives no AA.\n
To change the depth value, set the constant DEPTH to any number, default is 2, 1 is fastest.\n
To change the resolution of the image (image window scales with resolution), set the consants PIXELHEIGHT and PIXELWEIGHT appropriately.\n
-> All constants can be found in lines 17-24 of renderImage.cpp.

Credits to:
  Robert. S for OBJ-Loader, found at: https://github.com/Bly7/OBJ-Loader
    Used to get information from .obj files.
  Lode Vandevenne for lodepng, found at: https://github.com/lvandeve/lodepng
    Used to form a .png image for output.
