#include <iostream>
#include <windows.h>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <random>
#include <time.h>
#include "OBJ_Loader.h" // This is an external library to load a .obj file, can be found here: https://github.com/Bly7/OBJ-Loader
#include "lodepng.h" // This is an external library to form a .png image, can be found here: https://github.com/lvandeve/lodepng
#include "Camera.h"
#include "Sphere.h"
#include "Light.h"

using namespace std;

const int SAMPLES = 2; // Value used to set how many subpixels to use for AA (SAMPLESxSAMPLES AA)
const int DEPTH = 2; // Value to define how deeply to calculate reflections and refractions
const double REFRAC = 1.0; // Refractive index of environment 1 for ~air
const int PIXELHEIGHT = 180, PIXELWIDTH = 240; // Ouput image size parameters
const int N = 100; // Falloff parameter for specular light
const double PUNY = 0.00000001; // Value to prevent maths failures
const string objectChoice = "/resources/lowpolybunny.obj";
const string outputFilename = "out.png";


vector<double> colourForRay(Camera, Light, vector<Light>, vector<Sphere>, vector<objl::Mesh>, vector<vector<array<objl::Vertex, 3>>>, vector<Sphere>, array<double, 3>, array<double, 3>, int, double);
bool isShadowed(array<double, 3>, vector<Sphere>, vector<objl::Mesh>, vector<vector<array<objl::Vertex, 3>>>, vector<Sphere>, array<double, 3>);
void outputImage(vector<vector<vector<int>>>, vector<unsigned char>);
vector<array<double, 3>> get3DPixLocs(Camera, int, int, int, int, int);
double dotProd(array<double, 3>, array<double, 3>);
array<double, 3> crossProd(array<double, 3>, array<double, 3>);
array<double, 3> vecMinus(array<double, 3>, array<double, 3>);
array<double, 3> vecAdd(array<double, 3>, array<double, 3>);
array<double, 3> vecMultiply(double, array<double, 3>);
double vecMagnitude(array<double, 3>);
array<double, 3> vecNormalise(array<double, 3>);

int main(int argc, char * argv[]) {
    // Could potentially get user to input all const values required instead of hardcode them

    vector<vector<vector<int>>> pixels(PIXELHEIGHT, vector<vector<int>>(PIXELWIDTH, vector<int>(4)));
    // Create image for pixels
    vector<unsigned char> image;
    vector<objl::Mesh> meshes;
    vector<Sphere> boundingSpheres;
    vector<vector<array<objl::Vertex, 3>>> meshTriangles;
    vector<Sphere> spheres;
    vector<Light> directionals;

    // Create camera object
    Camera camera({0,0,0}, {0,0,1}, {1,0,0}, {0,-1,0}, 1, 1, ((double) PIXELWIDTH / (double) PIXELHEIGHT));

    // Load .obj file(s) into the scene
    objl::Loader objLoader;

    if (objLoader.LoadFile("./" + objectChoice)) {
        cout << objectChoice << " loaded" << endl;
        for (int i = 0; i < objLoader.LoadedMeshes.size(); i++) {
            meshes.push_back(objLoader.LoadedMeshes[i]);
        }
    }
    else {
        cout << "File could not be found/loaded" << endl;
        return(1);
    }

    // Fill meshTriangles and boundingSpheres vectors
    for (int i = 0; i < meshes.size(); i++) {
                              //X, Y, Z
        array<double, 3> max = {0, 0, 0};
        array<double, 3> min = {0, 0, 0};
        // Fill meshTriangles vector with the triangles from the mesh
        vector<array<objl::Vertex, 3>> triangles;
        for (int j = 0; j < meshes[i].Vertices.size(); j+=3) {
            array<objl::Vertex, 3> tri;
            tri = {meshes[i].Vertices[j],
                    meshes[i].Vertices[j+1],
                    meshes[i].Vertices[j+2]};
            triangles.push_back(tri);
            if (j == 0) {
                max = {tri[0].Position.X,
                    tri[0].Position.Y,
                    tri[0].Position.Z};
                min = {tri[0].Position.X,
                    tri[0].Position.Y,
                    tri[0].Position.Z};
            }
            for (int k = 0; k < 3; k++) {
                if (tri[k].Position.X > max[0]) {
                    max[0] = tri[k].Position.X;
                }
                else if (tri[k].Position.X < min[0]) {
                    min[0] = tri[k].Position.X;
                }
                if (tri[k].Position.Y > max[1]) {
                    max[1] = tri[k].Position.Y;
                }
                else if (tri[k].Position.Y < min[1]) {
                    min[1] = tri[k].Position.Y;
                }
                if (tri[k].Position.Z > max[2]) {
                    max[2] = tri[k].Position.Z;
                }
                else if (tri[k].Position.Z < min[2]) {
                    min[2] = tri[k].Position.Z;
                }
            }
        }
        // Push all triangles for the mesh
        meshTriangles.push_back(triangles);

        // Finding centre of mesh/bounding sphere
        array<double, 3> meshCentre = {(max[0] + min[0])/2,
                (max[1] + min[1])/2,
                (max[2] + min[2])/2};
        // Finding radius of bounding sphere
        double radius = 0;
        // Finding largest distance from centre to each point of the bounding box
        if (vecMagnitude(vecMinus({max[0], max[1], max[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({max[0], max[1], max[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({max[0], min[1], max[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({max[0], min[1], max[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({max[0], max[1], min[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({max[0], max[1], min[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({max[0], min[1], min[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({max[0], min[1], min[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({min[0], max[1], max[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({min[0], max[1], max[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({min[0], min[1], max[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({min[0], min[1], max[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({max[0], min[1], min[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({max[0], min[1], min[2]}, meshCentre));
        }
        else if (vecMagnitude(vecMinus({min[0], min[1], min[2]}, meshCentre)) > radius) {
            radius = vecMagnitude(vecMinus({min[0], min[1], min[2]}, meshCentre));
        }
        // Adding bounding sphere
        cout << "Bounding ";
        boundingSpheres.push_back(Sphere(meshCentre, radius, {0,0,0}, {0,0,0}, {0,0,0}, 1, 1));
    }

    cout << "Number of meshes: " << meshes.size() << endl;

    // Create sphere object(s)
    spheres.push_back(Sphere({4.2,6,0}, 3, {0,0,255}, {0,0,255}, {255,255,255}, 1, 1));
    spheres.push_back(Sphere({0,5,-1.7}, 0.8, {100,100,100}, {100,100,100}, {255,255,255}, 1.5, 0.05));
    spheres.push_back(Sphere({0,5,-1.7}, 0.3, {200,200,255}, {200,200,255}, {255,255,255}, 1.33, 0.1));
    spheres.push_back(Sphere({0,16,2}, 5, {255,0,0}, {255,0,0}, {255,255,255}, 1, 1));

    //spheres.push_back(Sphere({0,20,-49}, 50, {150,150,150}, {150,150,150}, {255,255,255})); // Floor sphere
    spheres.push_back(Sphere({0,20,-53}, 50, {150,150,150}, {200,150,150}, {255,255,255}, 1, 1)); // Floor sphere other

    // Creating light object(s)
    Light ambient = Light(0.2); // Ambient
    directionals.push_back(Light(0.6, {0, 1, -1})); // Directional
    directionals.push_back(Light(0.3, {3, 3, 1}));  // Directional
    directionals.push_back(Light(0.8, {3, -3, -1}, {-15, 15, 5})); // Directional point light

    cout << "reflection/refraction depth: " << DEPTH <<", samples per pixel: " << SAMPLES << ", image resolution: " << PIXELHEIGHT << "x" << PIXELWIDTH << endl;
    cout << "outputting to " << outputFilename << endl;
    cout << "Performing ray tests.. " << endl;
    // Iterate through each pixel and test for object
    for (int i = 0; i < PIXELHEIGHT; i ++) {
        for (int j = 0; j < PIXELWIDTH; j++) {
            vector<array<double, 3>> qs = get3DPixLocs(camera, i, j, PIXELHEIGHT, PIXELWIDTH, SAMPLES);
            vector<double> totalSamples = {0, 0, 0, 0};
            array<double, 3> d = vecMinus(get3DPixLocs(camera, i, j, PIXELHEIGHT, PIXELWIDTH, 1)[0], camera.focalPoint);
            // Anti-Aliasing / Super sampling
            for (int k = 0; k < qs.size(); k++) {
                // Calculate 3D Pixel Location (q) and direction vector (d)
                array<double, 3> q = qs[k];
                // Find pixel RGBa value for the ray from q with direction d
                vector<double> colourRay = colourForRay(camera, ambient, directionals, spheres, meshes, meshTriangles, boundingSpheres, d, q, DEPTH, REFRAC);
                // Adding to running total
                totalSamples = {totalSamples[0] + colourRay[0],
                    totalSamples[1] + colourRay[1],
                    totalSamples[2] + colourRay[2],
                    totalSamples[3] + colourRay[3]};
            }
            // Adding averages RGBa values to pixels vector
            pixels[i][j] = {(int) (totalSamples[0]/qs.size()), (int) (totalSamples[1]/qs.size()), (int) (totalSamples[2]/qs.size()), (int) (totalSamples[3]/qs.size())};
        }
        cout << "\rrow: " << i+1 << " completed";
    }
    cout << endl;

    cout << "Done" << endl;

    outputImage(pixels, image);
    Beep(523, 300);
    return(0);
}

// Helper Functions

vector<double> colourForRay(Camera camera, Light ambient, vector<Light> directionals, vector<Sphere> spheres, vector<objl::Mesh> meshes, vector<vector<array<objl::Vertex, 3>>> triangles, vector<Sphere> bSpheres, array<double, 3> d, array<double, 3> q, int depth, double refracIndex) {
    vector<double> pixel;
    pixel = {0, 0, 0, 255};

    // Boolean to inform system when no t has been found yet
    bool firstT = true;
    double bestT;

    // Normalising d vector
    d = vecNormalise(d);
    // Now have ray direction

    // Sphere interception
    // For all spheres in the scene
    for (int k = 0; k < spheres.size(); k++) {
        array<double, 3> p = vecMinus(q, spheres[k].centre);

        // Using discriminant to see if ray intercepts sphere
        double b = 2 * dotProd(d, p);
        double c = dotProd(p, p) - pow(spheres[k].radius,2);
        double discriminant = pow(b,2) - 4*c;

        // Check if ray intercepts Sphere
        if (discriminant > 0) {
            // Find t values and pick one
            double t;
            // It will intercept twice, use closest
            if (-b < sqrt(discriminant)) {
                t = (-b + sqrt(discriminant))/2;
            }
            else {
                t = (-b - sqrt(discriminant))/2;
            }
            // Use t value
            // Check if this is the firstT value found (Used to update t initially)
            if (firstT && t > 0) {
                bestT = t;
                firstT = false;
            }
            if (t <= bestT && t > 0) {
                // Update best t value if t is better than previous t
                bestT = t;

                // Finding point on sphere that intercepts ray
                array<double, 3> point = vecAdd(q, vecMultiply(t, d));
                // Find normal to sphere at point
                array<double, 3> normalS = vecMinus(point, spheres[k].centre);
                normalS = vecNormalise(normalS);


                // Calculate lighting
                // TODO make function to calculate lighting, lightAtSphere()

                // Adding lighting values
                array<double, 3> totalLight = {0, 0, 0};
                // For the ambient light
                array<double, 3> amb = vecMultiply(0.2, vecMultiply(ambient.intensity, spheres[k].ambient));
                // Add to total light
                totalLight = vecAdd(totalLight, amb);

                // For each directionals light source
                for (int l = 0; l < directionals.size(); l++) {
                    array<double, 3> lightD = vecMultiply(-1, directionals[l].direction);
                    lightD = vecNormalise(lightD);

                    // Testing for shadows
                    bool shadow = isShadowed(lightD, spheres, meshes, triangles, bSpheres, point);

                    // Begin finding diffusion and specular lighting
                    double nDotL = dotProd(normalS, lightD);

                    // Test to ensure no self-occlusion and not a shadowed pixel
                    if (nDotL >= 0 && !shadow) {
                        // Set RGB values from angle of incidence (nDotL) and diffuse values (colour from sphere)
                        array<double, 3> diff = vecMultiply(nDotL, spheres[k].diffuse);
                        // Add changes to RGB values from specular lighting
                        array<double, 3> h = vecMultiply(1/vecMagnitude(vecAdd(vecMultiply(-1, d), lightD)), vecAdd(vecMultiply(-1, d), lightD));
                        double nDotH = dotProd(normalS, h);
                        // Add specular lighting to diffuse lighting
                        array<double, 3> diffSpec = vecAdd(diff, vecMultiply(pow(nDotH,N), spheres[k].specular));
                        // Multiply RGB values by intensity of the light
                        array<double, 3> diffSpecL = vecMultiply(directionals[l].intensity, diffSpec);
                        // Add to total light
                        totalLight = vecAdd(totalLight, vecMultiply(spheres[k].dissolve, vecMultiply(0.5, diffSpecL)));
                        // Do depth of reflections
                        if (depth > 0) {
                            double dDotN = dotProd(d, normalS);
                            if (spheres[k].dissolve > 0) {
                                // Reflections
                                array<double, 3> reflection = vecMinus(d, vecMultiply(2*dDotN, normalS));
                                reflection = vecNormalise(reflection);
                                vector<double> reflectedColour = colourForRay(camera, ambient, directionals, spheres, meshes, triangles, bSpheres, reflection, vecAdd(point, vecMultiply(PUNY, reflection)), depth-1, refracIndex);
                                array<double, 3> vecRefColour = {reflectedColour[0], reflectedColour[1], reflectedColour[2]};
                                array<double, 3> ref = vecMultiply(spheres[k].dissolve, vecMultiply(0.1, vecRefColour));
                                totalLight = vecAdd(totalLight, ref);
                            }
                            // Refractions
                            if (1-spheres[k].dissolve > 0) {
                                double refractiveRatio = refracIndex/spheres[k].opticalDensity;
                                array<double, 3> refraction;
                                if (refractiveRatio == 1) {
                                    // When the refractive indexes are the same
                                    refraction = d;
                                }
                                else {
                                    double cosPhi = sqrt(1 - (pow(refractiveRatio, 2)*(1 - pow(dDotN, 2))));
                                    refraction = vecAdd(vecMultiply(refractiveRatio, d), vecMultiply((refractiveRatio*dDotN - cosPhi), normalS));
                                    refraction = vecNormalise(refraction);
                                }
                                vector<Sphere> newSpheres = spheres;
                                newSpheres.erase(newSpheres.begin() + k); // remove sphere hit from function call
                                vector<double> refractedColour = colourForRay(camera, ambient, directionals, newSpheres, meshes, triangles, bSpheres, refraction, vecAdd(point, vecMultiply(PUNY, refraction)), depth-1, spheres[k].opticalDensity);
                                array<double, 3> vecRefrColour = {refractedColour[0], refractedColour[1], refractedColour[2]};
                                array<double, 3> refr = vecMultiply(1-spheres[k].dissolve, vecRefrColour);
                                totalLight = vecAdd(totalLight, refr);
                            }
                        }
                    }
                }
                // If any of totalLight[i] are > 255 then set 255
                for (int l = 0; l < totalLight.size(); l++) {
                    if (totalLight[l] > 255) {
                        totalLight[l] = 255;
                    }
                }
                // Use totalLight RGB values and set pixels
                pixel = {totalLight[0], totalLight[1], totalLight[2], 255};

            }
        }

    }
    // end checking spheres

    // Mesh interception
    // For all meshes in the scene
    for (int k = 0; k < meshes.size(); k++) {
        // If ray intercepts bounding sphere
        array<double, 3> p = vecMinus(q, bSpheres[k].centre);

        // Using discriminant to see if ray intercepts sphere
        double b = 2 * dotProd(d, p);
        double c = dotProd(p, p) - pow(bSpheres[k].radius,2);
        double discriminant = pow(b,2) - 4*c;

        // Check if ray intercepts bounding sphere
        if (discriminant > 0) {
            for (int l = 0; l < triangles[k].size(); l++) {
                // Making normal vector for first vertex (and subsequently all vertices) in triangle
                array<double, 3> p2minusp1 = {triangles[k][l][1].Position.X - triangles[k][l][0].Position.X,
                    triangles[k][l][1].Position.Y - triangles[k][l][0].Position.Y,
                    triangles[k][l][1].Position.Z - triangles[k][l][0].Position.Z};
                array<double, 3> p3minusp1 = {triangles[k][l][2].Position.X - triangles[k][l][0].Position.X,
                    triangles[k][l][2].Position.Y - triangles[k][l][0].Position.Y,
                    triangles[k][l][2].Position.Z - triangles[k][l][0].Position.Z};
                array<double, 3> normal1 = crossProd(p2minusp1, p3minusp1);
                normal1 = vecNormalise(normal1);

                // Checking if the ray intercepts the plane (if n.d == 0, no interception)
                if (dotProd(normal1, d) != 0) {
                    // Ray does intercept plane

                    // Find t at interception
                    array<double, 3> pminuse = {triangles[k][l][0].Position.X - q[0],
                        triangles[k][l][0].Position.Y - q[1],
                        triangles[k][l][0].Position.Z - q[2]};
                    double t = (dotProd(pminuse, normal1)) / dotProd(normal1, d);

                    // Use t value
                    // Check if this is the firstT value found (Used to update t initially)
                    if (firstT && t > 0) {
                        bestT = t;
                        firstT = false;
                    }
                    if (t <= bestT && t > PUNY) {
                        // Check if bestT value is better than previous bestT

                        // Find if interception of plane is in the triangle
                        // Find vector, a at interception
                        array<double, 3> a = {q[0] + t * d[0],
                            q[1] + t * d[1],
                            q[2] + t * d[2]};

                        // Make vectors for triangle (vis)
                        array<double, 3> v1 = {triangles[k][l][1].Position.X - triangles[k][l][0].Position.X,
                            triangles[k][l][1].Position.Y - triangles[k][l][0].Position.Y,
                            triangles[k][l][1].Position.Z - triangles[k][l][0].Position.Z};
                        array<double, 3> v2 = {triangles[k][l][2].Position.X - triangles[k][l][1].Position.X,
                            triangles[k][l][2].Position.Y - triangles[k][l][1].Position.Y,
                            triangles[k][l][2].Position.Z - triangles[k][l][1].Position.Z};
                        array<double, 3> v3 = {triangles[k][l][0].Position.X - triangles[k][l][2].Position.X,
                            triangles[k][l][0].Position.Y - triangles[k][l][2].Position.Y,
                            triangles[k][l][0].Position.Z - triangles[k][l][2].Position.Z};

                        // Using cross product to check if ray intercepts triangle
                        array<double, 3> aminusp1 = {a[0] - triangles[k][l][0].Position.X,
                            a[1] - triangles[k][l][0].Position.Y,
                            a[2] - triangles[k][l][0].Position.Z};
                        array<double, 3> aminusp2 = {a[0] - triangles[k][l][1].Position.X,
                            a[1] - triangles[k][l][1].Position.Y,
                            a[2] - triangles[k][l][1].Position.Z};
                        array<double, 3> aminusp3 = {a[0] - triangles[k][l][2].Position.X,
                            a[1] - triangles[k][l][2].Position.Y,
                            a[2] - triangles[k][l][2].Position.Z};
                        array<double, 3> cross1 = crossProd(aminusp1, v1);
                        array<double, 3> cross2 = crossProd(aminusp2, v2);
                        array<double, 3> cross3 = crossProd(aminusp3, v3);

                        // If always in the same direction, then it intercepts.
                        if (dotProd(cross1, cross2) > 0 && dotProd(cross2, cross3) > 0
                            && dotProd(cross3, cross1) > 0) {
                            // Ray does intercept the triangle!
                            // Update bestT value
                            bestT = t;

                            // Already found point on triangle where ray intercepts (a)
                            // Already found normal to triangle at a (normal1)

                            // Calculate lighting
                            // TODO make function to calculate lighting, lightAtMesh()

                            // Adding lighting values
                            array<double, 3> totalLight = {0, 0, 0};
                            // For the ambient light
                            array<double, 3> meshAmbient = {meshes[k].MeshMaterial.Ka.X*255,
                                meshes[k].MeshMaterial.Ka.Y*255,
                                meshes[k].MeshMaterial.Ka.Z*255};
                            array<double, 3> amb = vecMultiply(0.2, vecMultiply(ambient.intensity, meshAmbient));
                            // Add to total light
                            totalLight = vecAdd(totalLight, amb);

                            // For each directional light source
                            for (int l = 0; l < directionals.size(); l++) {
                                array<double, 3> lightD = vecMultiply(-1, directionals[l].direction);
                                lightD = vecNormalise(lightD);

                                // Testing for shadows
                                bool shadow = isShadowed(lightD, spheres, meshes, triangles, bSpheres, a);

                                // Begin finding diffusion and specular lighting
                                double nDotL = dotProd(normal1, lightD);

                                // Test to ensure no self-occlusion and not a shadowed pixel
                                if (nDotL >= 0 && !shadow) {
                                    // Set RGB values from angle of incidence and diffuse values (colour)
                                    array<double, 3> meshDiffuse = {meshes[k].MeshMaterial.Kd.X*255,
                                        meshes[k].MeshMaterial.Kd.Y*255,
                                        meshes[k].MeshMaterial.Kd.Z*255};
                                    array<double, 3> meshSpecular = {meshes[k].MeshMaterial.Ks.X*255,
                                        meshes[k].MeshMaterial.Ks.Y*255,
                                        meshes[k].MeshMaterial.Ks.Z*255};
                                    array<double, 3> diff = vecMultiply(nDotL, meshDiffuse);
                                    // Add changes to RGB values from specular lighting
                                    array<double, 3> h = vecMultiply(1/vecMagnitude(vecAdd(vecMultiply(-1, d), lightD)), vecAdd(vecMultiply(-1, d), lightD));
                                    double nDotH = dotProd(normal1, h);
                                    array<double, 3> diffSpec = vecAdd(diff, vecMultiply(pow(nDotH, N), meshSpecular));
                                    // Multiply RGB values by intensity of the light
                                    array<double, 3> diffSpecL = vecMultiply(directionals[l].intensity, diffSpec);
                                    // Add to total light
                                    totalLight = vecAdd(totalLight, vecMultiply(meshes[k].MeshMaterial.d, vecMultiply(0.5, diffSpecL)));
                                    // Do depth of reflections
                                    if (depth > 0) {
                                        double dDotN = dotProd(d, normal1);
                                        if (meshes[k].MeshMaterial.d > 0) {
                                            // Reflections
                                            array<double, 3> reflection = vecMinus(d, vecMultiply(2*dDotN, normal1));
                                            reflection = vecNormalise(reflection);
                                            vector<double> reflectedColour = colourForRay(camera, ambient, directionals, spheres, meshes, triangles, bSpheres, reflection, vecAdd(a, vecMultiply(PUNY, reflection)), depth-1, refracIndex);
                                            array<double, 3> vecRefColour = {reflectedColour[0], reflectedColour[1], reflectedColour[2]};
                                            array<double, 3> ref = vecMultiply(meshes[k].MeshMaterial.d, vecMultiply(0.1, vecRefColour));
                                            totalLight = vecAdd(totalLight, ref);
                                        }
                                        if (1-meshes[k].MeshMaterial.d > 0) {
                                            // Refractions
                                            double refractiveRatio = refracIndex/meshes[k].MeshMaterial.Ni;
                                            array<double, 3> refraction;
                                            if (refractiveRatio == 1) {
                                                refraction = d;
                                            }
                                            else {
                                                double cosPhi = sqrt(1 - (pow(refractiveRatio, 2)*(1 - pow(dDotN, 2))));
                                                refraction = vecAdd(vecMultiply(refractiveRatio, d), vecMultiply((refractiveRatio*dDotN - cosPhi), normal1));
                                                refraction = vecNormalise(refraction);
                                            }
                                            vector<objl::Mesh> newMeshes = meshes;
                                            newMeshes.erase(newMeshes.begin() + k); // Remove mesh hit from function call
                                            vector<vector<array<objl::Vertex, 3>>> newTriangles = triangles;
                                            newTriangles.erase(newTriangles.begin() + k); // Remove triangles form mesh hit
                                            vector<Sphere> newBSpheres = bSpheres;
                                            newBSpheres.erase(newBSpheres.begin() + k); // Remove bounding sphere from mesh hit

                                            vector<double> refractedColour = colourForRay(camera, ambient, directionals, spheres, newMeshes, newTriangles, newBSpheres, refraction, vecAdd(a, vecMultiply(PUNY, refraction)), depth-1, meshes[k].MeshMaterial.Ni);
                                            array<double, 3> vecRefrColour = {refractedColour[0], refractedColour[1], refractedColour[2]};
                                            array<double, 3> refr = vecMultiply(1-meshes[k].MeshMaterial.d, vecRefrColour);
                                            totalLight = vecAdd(totalLight, refr);
                                        }
                                    }
                                }
                            }
                            // If any of totalLight[i] are > 255 then set 255
                            for (int l = 0; l < totalLight.size(); l++) {
                                if (totalLight[l] > 255) {
                                    totalLight[l] = 255;
                                }
                            }
                            // Add colour of triangle to the array of pixels
                            pixel = {totalLight[0], totalLight[1], totalLight[2], 255};
                        }
                    }
                }
            }
        }
    }
    // end checking mesh

    return(pixel);
}

bool isShadowed(array<double, 3> lightD, vector<Sphere> spheres, vector<objl::Mesh> meshes, vector<vector<array<objl::Vertex, 3>>> triangles, vector<Sphere> bSpheres, array<double, 3> point) {
    // Testing for shadows
    // Shadow testing for spheres
    for (int m = 0; m < spheres.size(); m++) {
        array<double, 3> p = vecMinus(point, spheres[m].centre);

        // Using discriminant to see if ray intercepts sphere
        double b = 2 * dotProd(lightD, p);
        double c = dotProd(p, p) - pow(spheres[m].radius,2);
        double discriminant = pow(b,2) - 4*c;
        if (discriminant > 0) {
            // Ray does intercept a sphere in two places
            // Therefore is in a shadow
            if ((-b - sqrt(discriminant))/2 > PUNY && (-b + sqrt(discriminant))/2 > PUNY) {
                return(true);
            }
        }
    }

    // Shadow testing for meshes
    for (int m = 0; m < meshes.size(); m++) {
        // If ray intercepts bounding sphere
        array<double, 3> p = vecMinus(point, bSpheres[m].centre);

        // Using discriminant to see if ray intercepts sphere
        double b = 2 * dotProd(lightD, p);
        double c = dotProd(p, p) - pow(bSpheres[m].radius,2);
        double discriminant = pow(b,2) - 4*c;

        // Check if ray intercepts bounding sphere
        if (discriminant > 0) {
            // Shadow testing for each triangle
            for (int n = 0; n < triangles[m].size(); n++) {
                // Making normal vector for first vertex (and subsequently all vertices) in triangle
                array<double, 3> p2minusp1 = {triangles[m][n][1].Position.X - triangles[m][n][0].Position.X,
                    triangles[m][n][1].Position.Y - triangles[m][n][0].Position.Y,
                    triangles[m][n][1].Position.Z - triangles[m][n][0].Position.Z};
                array<double, 3> p3minusp1 = {triangles[m][n][2].Position.X - triangles[m][n][0].Position.X,
                    triangles[m][n][2].Position.Y - triangles[m][n][0].Position.Y,
                    triangles[m][n][2].Position.Z - triangles[m][n][0].Position.Z};
                array<double, 3> normal = crossProd(p2minusp1, p3minusp1);
                normal = vecNormalise(normal);

                // Checking if the ray intercepts the plane (if n.d == 0, no interception)
                if (dotProd(normal, lightD) != 0) {
                    // Ray does intercept plane

                    // Find t at interception
                    array<double, 3> pminuse = {triangles[m][n][0].Position.X - point[0],
                        triangles[m][n][0].Position.Y - point[1],
                        triangles[m][n][0].Position.Z - point[2]};
                    double t = (dotProd(pminuse, normal)) / dotProd(normal, lightD);

                    if (t > PUNY) {
                        // Find if interception of plane is in the triangle
                        // Find vector, aa at interception
                        array<double, 3> a = {point[0] + t * lightD[0],
                            point[1] + t * lightD[1],
                            point[2] + t * lightD[2]};

                        // Make vectors for triangle (vis)
                        array<double, 3> v1 = {triangles[m][n][1].Position.X - triangles[m][n][0].Position.X,
                            triangles[m][n][1].Position.Y - triangles[m][n][0].Position.Y,
                            triangles[m][n][1].Position.Z - triangles[m][n][0].Position.Z};
                        array<double, 3> v2 = {triangles[m][n][2].Position.X - triangles[m][n][1].Position.X,
                            triangles[m][n][2].Position.Y - triangles[m][n][1].Position.Y,
                            triangles[m][n][2].Position.Z - triangles[m][n][1].Position.Z};
                        array<double, 3> v3 = {triangles[m][n][0].Position.X - triangles[m][n][2].Position.X,
                            triangles[m][n][0].Position.Y - triangles[m][n][2].Position.Y,
                            triangles[m][n][0].Position.Z - triangles[m][n][2].Position.Z};

                        // Using cross product to check if ray intercepts triangle
                        array<double, 3> aminusp1 = {a[0] - triangles[m][n][0].Position.X,
                            a[1] - triangles[m][n][0].Position.Y,
                            a[2] - triangles[m][n][0].Position.Z};
                        array<double, 3> aminusp2 = {a[0] - triangles[m][n][1].Position.X,
                            a[1] - triangles[m][n][1].Position.Y,
                            a[2] - triangles[m][n][1].Position.Z};
                        array<double, 3> aminusp3 = {a[0] - triangles[m][n][2].Position.X,
                            a[1] - triangles[m][n][2].Position.Y,
                            a[2] - triangles[m][n][2].Position.Z};
                        array<double, 3> cross1 = crossProd(aminusp1, v1);
                        array<double, 3> cross2 = crossProd(aminusp2, v2);
                        array<double, 3> cross3 = crossProd(aminusp3, v3);

                        // If always in the same direction, then it intercepts.
                        if (dotProd(cross1, cross2) > 0 && dotProd(cross2, cross3) > 0
                            && dotProd(cross3, cross1) > 0) {
                            // Ray does intercept the triangle!
                            // Therefore is in a shadow
                            return(true);
                        }
                    }
                }
            }
        }
    }
    // If is not in shadow, return false
    return(false);
}

void outputImage(vector<vector<vector<int>>> pixels, vector<unsigned char> image) {
    // Adding pixels to image
    for (int i = 0; i < PIXELHEIGHT; i++) {
        for (int j = 0; j < PIXELWIDTH; j++) {
            image.push_back(pixels[i][j][0]);
            image.push_back(pixels[i][j][1]);
            image.push_back(pixels[i][j][2]);
            image.push_back(pixels[i][j][3]);
        }
    }

    // Output pixels to a png file
    lodepng::encode(outputFilename, image, PIXELWIDTH, PIXELHEIGHT);
    image.clear();
    cout << "Render saved to .png" << endl;
}

vector<array<double, 3>> get3DPixLocs(Camera cam, int x, int y, int pixHeight, int pixWidth, int samples) {
    vector<array<double, 3>> qs;
    srand(time(NULL)); // Seeding random variable for jittering within grid
    uniform_real_distribution<> dist(0,1);
    for (int i = 1; i < 2*sqrt(samples); i+=2) {
        for (int j = 1; j < 2*sqrt(samples); j+=2) {
            double rand01;
            if (samples == 1) {
                rand01 = 0.5; // When sampling with 1 (no AA) just use 0.5
            }
            else {
                rand01 = ((double)rand()/(double)RAND_MAX); // Random variable
            }
            // Calculating parameters needed for q values using random variable
            double b = cam.imageHeight * ((x - ((rand01*2 + i-1)/2*sqrt(samples)))/pixHeight - (0.5));
            double r = cam.imageWidth * ((y - ((rand01*2 + j-1)/2*sqrt(samples)))/pixWidth - (0.5));

            array<double, 3> dist = {cam.focalPoint[0] + cam.focalLength * (-cam.forward[0]) + r * cam.right[0]
                - b * cam.up[0],
                cam.focalPoint[1] + cam.focalLength * (-cam.forward[1]) + r * cam.right[1]
                    - b * cam.up[1],
                cam.focalPoint[2] + cam.focalLength * (-cam.forward[2]) + r * cam.right[2]
                    - b * cam.up[2]};
            qs.push_back(dist);
        }
    }
    return(qs);
}

double dotProd(array<double, 3> vector1, array<double, 3> vector2) {
    return(vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2]);
}

array<double, 3> crossProd(array<double, 3> vector1, array<double, 3> vector2) {
    array<double, 3> result = {vector1[1] * vector2[2] - vector1[2] * vector2[1],
        vector1[2] * vector2[0] - vector1[0] * vector2[2],
        vector1[0] * vector2[1] - vector1[1] * vector2[0]};
    return(result);
}

array<double, 3> vecMinus(array<double, 3> first, array<double, 3> second) {
    array<double, 3> result = {first[0] - second[0], first[1] - second[1], first[2] - second[2]};
    return(result);
}

array<double, 3> vecAdd(array<double, 3> first, array<double, 3> second) {
    array<double, 3> result = {first[0] + second[0], first[1] + second[1], first[2] + second[2]};
    return(result);
}

array<double, 3> vecMultiply(double first, array<double, 3> second) {
    array<double, 3> result = {first * second[0], first * second[1], first * second[2]};
    return(result);
}

double vecMagnitude(array<double, 3> vector) {
    return(sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2)));
}

array<double, 3> vecNormalise(array<double, 3> vector) {
    array<double, 3> result = {vector[0] / sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2)),
        vector[1] / sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2)),
        vector[2] / sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2))};
    return(result);
}
