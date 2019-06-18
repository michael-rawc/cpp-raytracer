# cpp-raytracer
Raytrace based image renderer in C++

Basic, non-professionally optimised image renderer.
Hardcoded values in most places.

To compile:
  First, navigate to cpp-raytracer main file
  Then, execute: g++ src\renderImage.cpp src\lodepng.cpp -o renderImage
  Then, run: renderImage
  
Output image will be in out.png
  
To change which object to use, set the constant objectChoice to either "CubeTop.obj", "lowpolybunny.obj", or "bunny.obj"
To change the AA value, set the constant SAMPLES to any square number, 1 gives no AA.
To change the depth value, set the constant DEPTH to any number, default is 2, 1 is fastest.
To change the resolution of the image (image window scales with resolution), set the consants PIXELHEIGHT and PIXELWEIGHT appropriately.
-> All constants can be found in lines 17-24 of renderImage.cpp
