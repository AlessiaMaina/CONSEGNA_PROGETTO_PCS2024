#pragma once
#include <iostream>
#include <Eigen/Eigen>
#include <vector>
#include <array>

using namespace std;
using namespace Eigen;


namespace DiscreteAndFractureNetworkLibrary
{

struct strFractures
{
    unsigned int NumberOfFractures = 0;                      // Indica il numero delle fratture identificate
    vector<unsigned int> FractureId = {};                    // E' un vettore che indica l'identificativo di ogni frattura
    vector<unsigned int> NumberOfVertices = {};              // E' un vettore che indica il numero dei vertici per ogni frattura
    vector<vector<Vector3d>> VerticesOfFractures = {};       // E' un vettore che per ogni frattura memorizza un vettore di coordinate 3d
};

struct BoundingSphere
{
    Vector3d center;          // Indica le coordinate 3D del centro della Bounding Sphere
    double radius;            // Indica con double il valore che assume il raggio della Bounding Sphere
};

struct strTraces
{
    unsigned int NumberOfTraces = 0;                         // Indica il numero delle tracce identificate
    vector<unsigned int> TraceId = {};                       // Indica l'identificativo di ogni traccia
    vector<Vector2i> FractureIds = {};                       // Indica la coppia degli identificativi delle fratture che definiscono la traccia
    vector<array<Vector3d, 2>> VerticesOfTraces = {};        // E' un vettore che per ogni traccia memorizza un array contenente le coordinate 3d dei 2 estremi

    vector<array<bool,2>> Tips = {};                         // E' un vettore che memorizza per ogni traccia, i valori booleani che essa assume rispettivamente per la prima e per la seconda frattura
    vector<double> TraceLength = {};                         // E' un vettore che indica la lunghezza per ogni traccia
    vector<pair<double, unsigned int>> passing = {};         // E' un vettore che memorizza coppie lunghezza-id associate alle tracce passanti
    vector<pair<double, unsigned int>> notPassing = {};      // E' un vettore che memorizza coppie lunghezza-id associate alle tracce NON passanti
};


}

