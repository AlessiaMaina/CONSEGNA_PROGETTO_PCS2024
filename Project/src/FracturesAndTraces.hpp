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
    unsigned int NumberOfFractures = 0;                                 // Indica il numero delle fratture identificate
    vector<unsigned int> FractureId = {};                               // Indica l'identificativo di ogni frattura
    vector<unsigned int> NumberOfVertices = {};                         // Indica il numero dei vertici di una frattura
    map <unsigned int, vector<Vector3d>> VerticesOfFractures = {};      // Indica una struttura con una chiave (numero intero) e con dei valori, che sono le coordinate 3D dei vertici della frattura
};

struct BoundingSphere
{
    Vector3d center;          // Indica le coordinate 3D del centro della Bounding Sphere
    double radius;            // Indica con double il valore che assume il raggio della Bounding Sphere
};

struct strTraces
{
    unsigned int NumberOfTraces = 0;                                    // Indica il numero delle tracce identificate
    vector<unsigned int> TraceId = {};                                  // Indica l'identificativo di ogni traccia
    vector<Vector2i> FractureIds = {};                                  // Indica la coppia degli identificativi delle fratture che definiscono la traccia
    map<unsigned int, array<Vector3d, 2>> VerticesOfTraces = {};        // Indica le coordinate 3D dei 2 vertici che identificano la traccia (punto di inizio e punto di fine)

    vector<unsigned int> FracturesWithTraceId = {};                     // Indica gli identificativi delle fratture che contengono tracce
    vector<bool> Tips = {};                                             // Indica un vettore che contiene valori booleani per le Tips
    vector<double> TraceLength = {};                                    // Indica le lunghezze delle tracce, di tipo "double"
};


}

