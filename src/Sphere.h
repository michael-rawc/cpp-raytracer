#pragma once

#include <iostream>
#include <array>

class Sphere {
public:
    std::array<double, 3> centre;
    double radius;
    std::array<double, 3> ambient;
    std::array<double, 3> diffuse;
    std::array<double, 3> specular;
    double opticalDensity;
    double dissolve;

    Sphere(std::array<double, 3> middle, double poofiness, std::array<double, 3> kA, std::array<double, 3> kD,
            std::array<double, 3> kS, double opt, double diss) {
        centre = middle;
        radius = poofiness;
        ambient = kA;
        diffuse = kD;
        specular = kS;
        opticalDensity = opt;
        dissolve = diss;

        std::cout << "Sphere intialised at: " << centre[0] << ", " <<
            centre[1] << ", " << centre[2] << " with radius: " << radius << std::endl;
    }
};
