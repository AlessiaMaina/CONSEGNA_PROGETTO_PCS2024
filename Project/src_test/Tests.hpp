#pragma once

#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "FracturesAndTraces.hpp"
#include "Utils.hpp"
#include "PolygonalMesh.hpp"
#include "SortingAlgorithm_MERGESORT.hpp"


using namespace std;
using namespace Eigen;
using namespace DiscreteAndFractureNetworkLibrary;


// TEST sull'importo delle fratture
TEST(TESTImportFractures, FileInputValido)
{
    strFractures fractures;
    string inputFilePath = "./DFN/FR3_data.txt";
    bool result = importListFractures(fractures, inputFilePath);
    EXPECT_TRUE(result);
    EXPECT_EQ(fractures.NumberOfFractures, 3);
    EXPECT_EQ(fractures.FractureId[0], 0);
    EXPECT_EQ(fractures.NumberOfVertices[0], 4);

    EXPECT_TRUE(areVectorsEqual(fractures.VerticesOfFractures[fractures.FractureId[0]][0], Vector3d(0, 0, 0)));
    EXPECT_TRUE(areVectorsEqual(fractures.VerticesOfFractures[fractures.FractureId[0]][1], Vector3d(1, 0, 0)));
    EXPECT_TRUE(areVectorsEqual(fractures.VerticesOfFractures[fractures.FractureId[0]][2], Vector3d(1, 1, 0)));
    EXPECT_TRUE(areVectorsEqual(fractures.VerticesOfFractures[fractures.FractureId[0]][3], Vector3d(0, 1, 0)));
}

TEST(TESTImportFractures, FileInputNonValido)
{
    strFractures fractures;
    string inputFilePath = "./DFN/FR20_data.txt";      // Esempio di un file inesistente
    bool result = importListFractures(fractures, inputFilePath);
    EXPECT_FALSE(result);
}

// TEST sulla funzione che calcola le Bounding Spheres+
TEST(TESTComputeBoundingSphere, SetVerticiValido)
{
    // Definizione di un set di vertici noto
    vector<Vector3d> vertices =
    {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };

    // Calcolo della Bounding Sphere
    BoundingSphere result;
    bool exceptionThrown = false;
    try
    {
        result = computeBoundingSphere(vertices);
    } catch (...)
    {
        exceptionThrown = true;
    }

    // Verifica che nessuna eccezione sia stata lanciata
    EXPECT_FALSE(exceptionThrown);

    // Verifica che il centro e il raggio della sfera siano quelli attesi
    Vector3d expectedCenter(0.5, 0.5, 0.0);
    double expectedRadius = 0.7071;  // Calcolo del raggio approssimato

    EXPECT_TRUE(result.center.isApprox(expectedCenter, 0.001));
    EXPECT_NEAR(result.radius, expectedRadius, 0.001);      // Utilizzo tolleranza pari a 0.001
}


TEST(TESTComputeBoundingSphere, SetVerticiSenzaValori)
{
    vector<Vector3d> vertices;
    EXPECT_THROW(computeBoundingSphere(vertices), invalid_argument);
}


TEST(TESTComputeBoundingSphere, SetCon2Vertici)
{
    vector<Vector3d> vertices = {Vector3d(0, 0, 0), Vector3d(1, 0, 0)};
    BoundingSphere result = computeBoundingSphere(vertices);
    // Calcolo del centro atteso
    Vector3d expectedCenter = (vertices[0] + vertices[1]) / 2.0;
    // Calcolo del raggio atteso
    double expectedRadius = (vertices[0] - vertices[1]).norm() / 2.0;
    EXPECT_TRUE(result.center.isApprox(expectedCenter, 1e-10));         // Utilizzo tolleranza pari a 10^(-10)
    EXPECT_NEAR(result.radius, expectedRadius, 1e-10);
}


// TEST sulla funzione che verifica se un punto si trova nella frattura
TEST(TESTIsPointInFracture, PuntoInterno)
{
    vector<Vector3d> fractureVertices =
    {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };

    Vector3d pointInside(0.5, 0.5, 0.0);
    EXPECT_TRUE(isPointInFracture(pointInside, fractureVertices));
}

TEST(TESTIsPointInFracture, PuntoEsterno)
{
    vector<Vector3d> fractureVertices =
    {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };

    Vector3d pointOutside(2.0, 2.0, 0.0);
    EXPECT_FALSE(isPointInFracture(pointOutside, fractureVertices));
}

// TEST per verificare il vettore normale
TEST(TESTnormal, Vertici)
{
    vector<Vector3d> vertices =
    {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(1.0, 0.0, 0.0),
        Vector3d(1.0, 1.0, 0.0),
        Vector3d(0.0, 1.0, 0.0)
    };

    Vector3d expectedNormal(0.0, 0.0, 1.0);   // Vettore normale atteso per questi vertici

    Vector3d computedNormal = normal(vertices);     // Calcolo del vettore normale

    // Verifica che il vettore calcolato e quello atteso siano uguali
    EXPECT_TRUE(areVectorsEqual(computedNormal, expectedNormal));
}

// TEST per verificare il calcolo effettivo delle tracce
TEST(TESTDefinitionOfTraces, Generico)
{
    // Test eseguito per la traccia tra la frattura ID 0 e la frattura ID 1
    strFractures frac;
    strTraces trac;
    frac.NumberOfFractures = 2;
    frac.FractureId = {0, 1};
    frac.NumberOfVertices = {4, 4};
    frac.VerticesOfFractures[0] = {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)};
    frac.VerticesOfFractures[1] = {Vector3d(0.8, 0.0, -0.1), Vector3d(0.8, 0.0, 2.9999999999999999e-01), Vector3d(0.8, 1, 2.9999999999999999e-01), Vector3d(0.8, 1, -1)};
    trac.NumberOfTraces = 0;

    definitionOfTraces(frac, trac);

    array<Vector3d, 2> VerticesOfTraces = {Vector3d(0.8, 0.0, 1), Vector3d(0.0, 0.8, 0.0)};
    EXPECT_EQ(trac.NumberOfTraces, 1);
    EXPECT_EQ(trac.TraceId.size(), 1);
    EXPECT_EQ(trac.TraceId[0], 0);
    EXPECT_EQ(trac.FractureIds.size(), 1);
    EXPECT_EQ(trac.FractureIds[0], Vector2i(0, 1));
    EXPECT_EQ(trac.VerticesOfTraces[0], VerticesOfTraces);
}


// TEST per verificare se un punto si trova o meno su un segmento
TEST(TESTIsPointInSegment, PointOnSegment)
{
    Vector3d point(0.5, 0.0, 0.0); // Punto che si trova effettivamente sul segmento
    Vector3d p1(0.0, 0.0, 0.0);
    Vector3d p2(1.0, 0.0, 0.0);

    ASSERT_TRUE(isPointInSegment(point, p1, p2));
}

TEST(TESTIsPointInSegment, PointNotOnSegment)
{
    Vector3d point(0.5, 0.1, 0.0); // Punto che non si trova sul segmento
    Vector3d p1(0.0, 0.0, 0.0);
    Vector3d p2(1.0, 0.0, 0.0);

    ASSERT_FALSE(isPointInSegment(point, p1, p2));
}

// Test case per verificare il comportamento della funzione computeTypeTrace
TEST(TESTComputeTypeTrace, NonPassingTrace)
{
    strFractures frac;
    strTraces trac;

    // Configura un esempio di frattura e tracce dove la traccia non è passante
    frac.NumberOfFractures = 1;
    frac.NumberOfVertices.push_back(4);
    frac.VerticesOfFractures.insert({
        1,
        {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)}
    });

    trac.NumberOfTraces = 1;
    trac.VerticesOfTraces.insert({
        1,
        {{Vector3d(0.2, 0.2, 0), Vector3d(0.8, 0.8, 0)}}
    });

    computeTypeTrace(frac, trac);

    // Verifica che il risultato sia conforme alle aspettative
    ASSERT_EQ(trac.Tips.size(), 1);
    EXPECT_FALSE(trac.Tips[0]);     // Il valore atteso è che la traccia NON sia passante
}


TEST(TESTComputeTypeTrace, PassingTrace)
{
    strFractures frac;
    strTraces trac;

    // Configura un esempio di frattura e tracce dove la traccia è passante
    frac.NumberOfFractures = 1;
    frac.NumberOfVertices.push_back(4);
    frac.VerticesOfFractures.insert({
        1,
        {Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(1, 1, 0), Vector3d(0, 1, 0)}
    });

    trac.NumberOfTraces = 1;
    trac.VerticesOfTraces.insert({
        1,
        {{Vector3d(0.2, 0.2, 0), Vector3d(0.8, 0.2, 0)}}
    });

    // Esegue la funzione computeTypeTrace sul dato di test
    computeTypeTrace(frac, trac);

    // Verifica che il risultato sia conforme alle aspettative
    ASSERT_EQ(trac.Tips.size(), 1);
    EXPECT_TRUE(trac.Tips[0]); // Ci si aspetta che la traccia sia passante
}


// TEST per verificare il calcolo della lunghezza delle tracce
TEST(TESTOrderTraces, CasoBase)
{
    strTraces trac;

    // Configura un esempio di tracce, ovvero le tracce che sono presenti sulla frattura ID 0
    trac.NumberOfTraces = 2;
    trac.VerticesOfTraces =
    {
        {0, {Vector3d(0, 0, 0), Vector3d(1, 0, 0)}},
        {1, {Vector3d(0, 0, 0), Vector3d(0.3, 0.3, 0)}}
    };

    trac.FractureIds = {Vector2i(0, 0), Vector2i(0, 0)};
    trac.Tips = {false, true};
    trac.TraceLength = {1.0, 0.3161837};
    trac.TraceId = {0, 1};

    orderTraces(trac);

    // Verifica che il risultato sia conforme alle aspettative
    ASSERT_EQ(trac.Tips.size(), 2);
    // Le tracce risultano essere ordinate per lunghezza decrescente
    EXPECT_EQ(trac.TraceId[0], 0);
    EXPECT_EQ(trac.TraceId[1], 1);
    EXPECT_EQ(trac.Tips[0], false);
    EXPECT_EQ(trac.Tips[1], true);
}









