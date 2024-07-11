#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary
{

struct strPolygonalMesh
{

    // Definizione delle Celle 0D
    unsigned int NumberOfCell0Ds = 0;                           // Indica il numero di Celle 0D
    vector<unsigned int> IdOfCell0D = {};                       // Indica l'identificativo delle Celle 0D
    vector<Vector3d> CoordinatesOfCell0Ds = {};                 // Indica un vettore contenente coordinate 3D relativo alle Celle 0D

    // Definizione delle Celle 1D
    unsigned int NumberOfCell1Ds = 0;                           // Indica il numero di Celle 1D
    vector<unsigned int> IdOfCell1D = {};                       // Indica l'identificativo delle Celle 1D
    vector<Vector2i> VerticesOfCell1DIds = {};                  // Indica un vettore che contiene gli indici dei vertici che compongono i segmenti delle Celle 1D

    // Definizione delle Celle 2D
    unsigned int NumberOfCell2Ds = 0;                           // Indica il numero di Celle 2D
    vector<unsigned int> IdOfCell2D = {};                       // Indica l'identificativo delle Celle 2D
    vector<vector<unsigned int>> VerticesOfCell2DIds = {};      // Indica gli identificativi dei vertici che definiscono la Cella 2D
    vector<vector<unsigned int>> EdgesOfCell2DIds = {};         // Indica gli identificativi dei lati che definiscono la Cella 2D
    vector<unsigned int> NumberOfVertices = {};                 // Indica il numero dei vertici che definiscono la Cella 2D
    vector<unsigned int> NumberOfEdges = {};                    // Indica il numero dei lati che definiscono la Cella 2D

};

}

