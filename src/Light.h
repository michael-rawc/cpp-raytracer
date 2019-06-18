#pragma once

#include <iostream>
#include <array>

class Light {
public:
    double intensity;
    std::array<double, 3> direction;
    std::array<double, 3> point;

    Light (double intense) {
        intensity = intense;

        std::cout << "Ambient light intialised with intensity: "<< intensity << std::endl;
    }

    Light (double intense, std::array<double, 3> d) {
        intensity = intense;
        direction = d;

        std::cout << "Directional light initialised with intensity: " << intensity << ", direction: "
            << direction[0] << ", " <<
                direction[1] << ", " << direction[2] << std::endl;
    }

    Light (double intense, std::array<double, 3> d, std::array<double, 3> v) {
        intensity = intense;
        direction = d;
        point = v;

        std::cout << "Directional light initialised with intensity: " << intensity << ", direction: "
            << direction[0] << ", " <<
                direction[1] << ", " << direction[2] <<
                ", at: " << point[0] << ", " <<
                    point[1] << ", " << point[2] << std::endl;
    }
};
