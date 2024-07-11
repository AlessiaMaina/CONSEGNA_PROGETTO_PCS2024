#pragma once

#include "FracturesAndTraces.hpp"
#include "importExport.hpp"
#include "PolygonalMesh.hpp"



using namespace std;
using namespace DiscreteAndFractureNetworkLibrary;
using namespace PolygonalLibrary;


bool testEdgesOfFracture(const strFractures& fractures);

bool testVerticesOfFracture(const strFractures& fractures);

bool areVectorsEqual(const Vector3d& v1,
                     const Vector3d& v2);

BoundingSphere computeBoundingSphere(const vector<Vector3d>& vertices);

bool proximityOfFractures(const BoundingSphere& sphere1,
                          const BoundingSphere& sphere2,
                          const double radiusTolerance);

Vector3d normal(const vector<Vector3d>& fracture);

double calculateTriangleArea(const Vector3d& A,
                             const Vector3d& B,
                             const Vector3d& C);

bool isPointInTriangle(const Vector3d& P,
                       const Vector3d& A,
                       const Vector3d& B,
                       const Vector3d& C);

bool isPointInFracture(const Vector3d& P,
                       vector<Vector3d>& polyVert);

VectorXd solveWithDecompositionPALU(const MatrixXd& matrixA,
                                    const VectorXd& columnVector);

VectorXd solveWithDecompositionQR(const MatrixXd& matrixA,
                                  const VectorXd& columnVector);

void definitionOfTraces(strFractures& frac,
                        strTraces& trac);

bool isPointInSegment(const Vector3d& point,
                      const Vector3d& p1,
                      const Vector3d& p2);

void computeTypeTrace(const strFractures &fractures,
                      strTraces& traces);

void lengthTraces(strTraces& traces);

void orderTraces(strTraces& traces);





