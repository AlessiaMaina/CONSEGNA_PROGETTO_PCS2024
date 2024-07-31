#pragma once

#include "FracturesAndTraces.hpp"

using namespace std;
using namespace DiscreteAndFractureNetworkLibrary;

#include "UCDUtilities.hpp"

using namespace Gedim;

void exportParaview(strFractures& frac, strTraces& trac, const string& fileFractures, const string& fileTraces)
{
    // ESPORTAZIONE IN PARAVIEW DELLE FRATTURE
    unsigned int N = 0;

    // Ciclo FOR che serve per calcolare il numero TOTALE di vertici, prendendo in considerazione ogni frattura
    for (size_t i = 0; i < frac.VerticesOfFractures.size(); ++i)
    {
        // Indica il numero di vertici per la frattura corrente
        N += frac.VerticesOfFractures[i].size();
    }

    // Creazione di una matrice di dimensione 3 x N, in cui ogni colonna corrisponde alle coordinate
    // di un vertice e ogni riga corrisponde rispettivamente alle coordinate x y z
    MatrixXd frac_points(3,N);

    unsigned int col = 0;

    // Inserimento delle coordinate di tutti i vertici all'interno della matrice

    for (const auto& vertex : frac.VerticesOfFractures)
    {
        // Per ogni vertice, vengono memorizzate le coordinate 3d nelle relative righe
        for (const auto& p : vertex)
        {
            frac_points(0, col) = p.x();
            frac_points(1, col) = p.y();
            frac_points(2, col) = p.z();

            ++col;
        }
    }

    // Creazione di un vettore che al suo interno contiene un vettore con gli id di ogni vertice
    vector<vector<unsigned int>> idVerticesAllFrac;
    unsigned int start = 0;

    // Ciclo FOR per ogni numero di vertici relativo a ogni singola frattura
    for (unsigned int n : frac.NumberOfVertices)
    {
        // Creazione di un vettore per memorizzare gli id dei vertici della frattura corrente
        vector<unsigned int> idVerticesSingleFrac;

        // Per ogni vertice, viene inserito il relativo id all'interno di idVerticesSingleFrac
        for (unsigned int i = 0; i < n; ++i)
        {
            idVerticesSingleFrac.push_back(start + i);
        }

        // Per ogni n, il relativo idVerticesSingleFrac viene inserito in idVerticesAllFrac
        idVerticesAllFrac.push_back(idVerticesSingleFrac);

        start += n;
    }

    vector<UCDProperty<double>> points_properties;
    vector<UCDProperty<double>> frac_properties;

    VectorXi frac_materials(idVerticesAllFrac.size());

    // Ciclo FOR che ad ogni poligono associa un indice sequenziale (materiale)
    for (int i = 0; i < frac_materials.size(); ++i)
    {
        frac_materials(i) = i;
    }

    Gedim::UCDUtilities UCD;
    UCD.ExportPolygons(fileFractures, frac_points, idVerticesAllFrac, points_properties, frac_properties, frac_materials);


    // ESPORTAZIONE IN PARAVIEW DELLE TRACCE

    size_t n = trac.VerticesOfTraces.size();
    size_t numDouble = n * 2;

    // Creazione di una matrice di dimensione 3 x numero totale di vertici che definiscono ogni traccia,
    // dove ogni coppia di colonne indica rispettivamente il primo e il secondo vertice della traccia
    // dove ogni riga corrisponde rispettivamente alle coordinate x y z
    MatrixXd trac_points(3, numDouble);

    for (size_t i = 0; i < n; ++i)
    {
        // Il primo vertice di ogni traccia viene inserito nelle colonne pari
        trac_points.col(i * 2) = trac.VerticesOfTraces[i][0];
        // Il secondo vertice di ogni traccia viene inserito nelle colonne dispari
        trac_points.col(i * 2 + 1) = trac.VerticesOfTraces[i][1];
    }

    // Creazione di una matrice che al suo interno memorizza gli indici rappresentanti gli estremi di ogni traccia
    MatrixXi indexOfEdges(2, n);
    for (size_t i = 0; i < n; ++i)
    {
        // L'id del primo vertice di ogni traccia viene memorizzato nella prima riga
        indexOfEdges(0, i) = i * 2;
        // L'id del secondo vertice di ogni traccia viene memorizzato nella seconda riga
        indexOfEdges(1, i) = i * 2 + 1;
    }

    VectorXi trac_materials(n);

    // Ciclo FOR che ad ogni traccia associa un indice sequenziale (materiale)
    for (unsigned int i = 0; i < n; ++i)
    {
        trac_materials(i) = i;
    }

    UCD.ExportSegments(fileTraces, trac_points, indexOfEdges, points_properties, frac_properties, trac_materials);
}
