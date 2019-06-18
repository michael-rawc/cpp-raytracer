#pragma once

#include <iostream>
#include <array>

class Camera {
public:
    std::array<double, 3> focalPoint;
    std::array<double, 3> up;
    std::array<double, 3> right;
    std::array<double, 3> forward;
    double focalLength;
    double imageHeight;
    double imageWidth;

    Camera (std::array<double, 3> focalP, std::array<double, 3> upD, std::array<double, 3> rightD,
         std::array<double, 3> forwardD, double focalL, double imageH, double imageW) {
        focalPoint = focalP;
        up = upD;
        right = rightD;
        forward = forwardD;
        focalLength = focalL;
        imageHeight = imageH;
        imageWidth = imageW;

        std::cout << "Camera intialised at: " << focalPoint[0] << ", " <<
            focalPoint[1] << ", " << focalPoint[2] << std::endl;
    }
};
