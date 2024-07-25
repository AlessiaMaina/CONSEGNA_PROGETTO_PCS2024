#include "Utils.hpp"
#include "importExport.hpp"
#include "SortingAlgorithm_MERGESORT.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <Eigen/Eigen>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>



using namespace std;
using namespace Eigen;
using namespace DiscreteAndFractureNetworkLibrary;
using namespace SortLibrary;


// *************************************************************************************************************************************
// PRIMA PARTE DEL PROGETTO ************************************************************************************************************
// *************************************************************************************************************************************



const double machineEpsilon = numeric_limits<double>::epsilon();          // Definizione della tolleranza
const double defaultTolerance = pow(10,-10);

// TEST 1: Verificare nella frattura la presenza di lati tutti diversi da zero
bool testEdgesOfFracture(const strFractures& fractures)
{
    for (const auto& fracture : fractures.VerticesOfFractures)
    {
        const auto& vertices = fracture.second;
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            Vector3d v_Start = vertices[i];
            Vector3d v_End = vertices[i + 1];
            if ((v_End - v_Start).norm() <= machineEpsilon)         // La norma della differenza tra i 2 vertici NON deve essere minore o uguale alla tolleranza
            {
                cerr << "ATTENTION: Fracture has sides of zero lenght. This is impossible! " << endl;
                return false;
            }
        }
    }

    return true;
}


// TEST 2: Verificare nella frattura che sono presenti almeno 3 vertici, altrimenti la frattura NON è definita
bool testVerticesOfFracture(const strFractures& fractures)
{
    for (const auto& fracture : fractures.VerticesOfFractures)
    {
        size_t numberOfVertices = fracture.second.size();
        if (numberOfVertices < 3)
        {
            cerr << "ATTENTION: Fracture has less than 3 vertices and is undefined. This is impossible!" << endl;
        }
    }

    return true;
}

// Definizione di una funzione serve a verificare se due vettori 3d sono "quasi" uguali
bool areVectorsEqual(const Vector3d& v1, const Vector3d& v2)
{
    if((v1-v2).squaredNorm()<= machineEpsilon)
    {
        return true;
    }
    return false;
}

double distance(const Vector3d& a, const Vector3d& b)
{
    double distAB = (a - b).norm();
    return distAB;
}

BoundingSphere computeBoundingSphere(const vector<Vector3d>& vertices)
{
    if (vertices.empty())
    {
        // Se il vettore di coordinate 3D "vertices" è vuoto, viene generata un'eccezione
        throw invalid_argument("ERROR: vector<Vector3d> vertices is empty");
    }

    // Dato un vertice, viene trovato il punto più distante da esso che presenta una distanza massima
    Vector3d x0 = vertices[0];
    Vector3d farVert = vertices[0];
    double maxDist = 0.0;

    for (const auto& vertex : vertices)
    {
        double dist = distance(vertex, x0);

        if (dist > maxDist)
        {
            // Aggiornamento della distanza massima e del punto più lontano
            maxDist = dist;
            farVert = vertex;
        }
    }

    // Dato farthestPoint, trova il punto più lontano da esso per che presenta una distanza massima
    Vector3d secFarVert = farVert;
    double maxDistance = 0.0;

    for (const auto& vertex : vertices)
    {
        // Calcolo della distanza tra un vertice e quello più lontano
        double dist = distance(vertex, farVert);
        if (dist > maxDistance)
        {
            // Aggiornamento delle variabili
            maxDistance = dist;
            secFarVert = vertex;
        }
    }

    // Definizione della sfera iniziale, avente centro pari alla media dei 2 punti più lontani
    // e avente come raggio la metà della distanza tra farthestPoint e secondFarthestPoint
    Vector3d center = (farVert + secFarVert) / 2.0;
    double radius = distance(farVert, secFarVert) / 2.0;

    // Espansione della sfera per includere tutti i vertici
    for (const auto& vertex : vertices)
    {
        double dist = distance(vertex, center);
        if (dist > radius)
        {
            double newRadius = (radius + dist) / 2.0;
            // Fattore di scala: è un fattore che serve per determinare l'espansione uniforme della sfera il cui centro
            // deve spostarsi in relazione al vertice più distante. Viene utilizzato per aggiornare iterativamente
            // la posizione del centroide della Bounding Sphere
            double scaleFactor = (newRadius - radius) / dist;

            center = (1 - scaleFactor) * center + scaleFactor * vertex;
            radius = newRadius;
        }
    }

    return {center, radius};
}


bool proximityOfFractures(BoundingSphere& sphere1, BoundingSphere& sphere2, const double radiusTolerance)
{
    double distBetweenCenters = distance(sphere1.center, sphere2.center);
    double radiusSum = sphere1.radius + sphere2.radius + radiusTolerance;

    // Se la distanza tra i centri delle 2 Bounding Spheres è minore o uguale della somma dei raggi,
    // inclusa la tolleranza, allora le Bounding Spheres si intersecano
    if (distBetweenCenters <= radiusSum)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Definizione del vettore normale, ovvero un vettore perpendicolare alla superficie della frattura
Vector3d normal(const vector<Vector3d>& fracture)
{
    Vector3d uVert = fracture[0];
    Vector3d commonVert = fracture[1];
    Vector3d vVert = fracture[2];

    Vector3d u = commonVert - uVert;        // Indica il vettore u
    Vector3d v = commonVert - vVert;        // Indica il vettore v

    Vector3d crossProd = (u.cross(v));      // Calcolo del prodotto vettoriale tra il vettore u e v
    Vector3d n = crossProd.normalized();
    // Il prodotto vettoriale deve essere normalizzato, al fine di ottenere il vettore normale al poligono

    return n;
}



// Definizione della funzione per calcolare l'area di un triangolo dati i vertici
double calculateTriangleArea(const Vector3d& A, const Vector3d& B, const Vector3d& C)
{
    // Calcolo dei vettori AB e AC
    Vector3d AB = B - A;
    Vector3d AC = C - A;

    // Calcolo del loro prodotto vettoriale, ovvero AB x AC
    Vector3d crossProduct = AB.cross(AC);

    // L'area del triangolo dati i vertici è pari alla metà della norma del prodotto vettoriale
    double area = crossProduct.norm() / 2.0;

    return area;
}

// Definizione di una funzione che verifica se un punto è interno o meno ad un triangolo
bool isPointInTriangle(const Vector3d& P, const Vector3d& A, const Vector3d& B, const Vector3d& C)
{
    double areaABC = calculateTriangleArea(A, B, C);

    double areaABP = calculateTriangleArea(A, B, P);
    double areaACP = calculateTriangleArea(A, C, P);
    double areaBCP = calculateTriangleArea(B, C, P);

    // Definizione della somma totale delle aree dei triangoli che si formano tra il punto P e i vertici A, B e C
    double sumAreas = areaABP + areaACP + areaBCP;

    // Il punto è considerato interno al triangolo se la differenza tra la somma totale delle aree e l'area di ABC
    // è inferiore alla tolleranza scelta
    if(abs(sumAreas - areaABC) < machineEpsilon)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// Definizione di una funzione che verifica se un punto è interno o meno a una frattura (ovvero un poligono)
// utilizzando la triangolazione, nel dettaglio la scomposizione in triangoli della frattura
// Si ricorda che un poligono con n vertici può essere suddiviso in "n-2" triangoli
bool isPointInFracture(const Vector3d& P, vector<Vector3d>& fractureVertices)
{
    for (unsigned int i = 1; i < fractureVertices.size()-1; ++i)
    {
    // Il primo vertice viene utilizzato come vertice fisso per originare i triangoli
	// Gli argomenti sono: Punto, Primo Vertice, Secondo Vertice, Terzo Vertice
	if (isPointInTriangle(P, fractureVertices[0], fractureVertices[i], fractureVertices[i+1]))
        {
            return true;
        }
    }
    return false;
}



VectorXd solveWithDecompositionPALU(const MatrixXd& matrixA, const VectorXd& columnVector)
{
    VectorXd solutionWithPALU = matrixA.partialPivLu().solve(columnVector);
    return solutionWithPALU;

}

VectorXd solveWithDecompositionQR(const MatrixXd& matrixA, const VectorXd& columnVector)
{
    VectorXd solutionWithQR = matrixA.colPivHouseholderQr().solve(columnVector);
    return solutionWithQR;

}

// Definizione di una funzione per definire tutti i parametri di una traccia
void definitionOfTraces(strFractures& frac, strTraces& trac)
{
    unsigned int count = 0;		// Inizializzazione di un contatore a zero

    // Utilizzando la struct BoundingSphere, si crea un vettore che viene inizializzato con dimensione pari
    // al numero delle fratture
    vector<BoundingSphere> boundingSpheres(frac.NumberOfFractures);

    // Si effettua il calcolo delle Bounding Spheres per tutte le fratture
    for (unsigned int i = 0; i < frac.NumberOfFractures; ++i)
    {
        boundingSpheres[i] = computeBoundingSphere(frac.VerticesOfFractures[i]);
    }

    for (unsigned int i = 0; i < frac.NumberOfFractures - 1; i++)
    {
        for (unsigned int j = i + 1; j < frac.NumberOfFractures; j++)
        {
            // Condizione che verifica se le Bounding Spheres sono sufficientemente vicine da creare un'intersezione
            if (proximityOfFractures(boundingSpheres[i], boundingSpheres[j], machineEpsilon))
            {
                // Definizione della matrice A e del vettore colonna b, per la risoluzione del sistema lineare
                // al fine di trovare l'intersezione delle 2 fratture, dati 2 piani passanti per esse
                
                // Definizione delle mappe di Gauss (versori normali) delle 2 fratture
                Vector3d n1 = normal(frac.VerticesOfFractures[i]);
                Vector3d n2 = normal(frac.VerticesOfFractures[j]);
                Vector3d t = n1.cross(n2);
                
                Matrix3d A;	// Matrice dove n1 è contenuto 1a riga, n2 nella 2a riga e infine t nella 3a riga
                A.row(0) = n1.transpose();
                A.row(1) = n2.transpose();
                A.row(2) = t.transpose();

                double b1 = n1.transpose() * boundingSpheres[i].center;
                double b2 = n2.transpose() * boundingSpheres[j].center;
                Vector3d b = {b1, b2, 0};

                // Se si verifica questa condizione, ovvero che il determinante è pari a 0, dove zero viene
                // approssimato alla tolleranza di default, e ciò implica che t,n1 e n2 siano complanari, il sistema
                // non ammette soluzione e di conseguenza non avviene l'intersezione
                // Complanarità: si verifica quando il prodotto misto è pari a 0, ovvero t*(n1xn2). Non viene indicata nella seguente condizione IF,
                // in qanto è sufficiente la condizione sul determinante

                if (A.determinant() < machineEpsilon)
                {
                    continue;
                }

                Vector3d extrTrace = solveWithDecompositionPALU(A,b);

                unsigned int NumOfExtremes = 0;	   // Inizializzazione di un contatore per il numero di estremi trovati
                array<Vector3d, 2> extremes;	   // Definizione di unArray per memorizzare gli estremi della traccia

                // Ciclo FOR che serve per eseguire i controlli sulla prima frattura. E' un ciclo che itera sui vertici
                // della prima frattura al fine di trovare gli estremi che definiscono la traccia di intersezione
                for (unsigned int k = 0; k < frac.NumberOfVertices[i]; k++)
                {
                    // Con extr1 ed extr2 vengono indicati 2 estremi consecutivi della frattura
                    Vector3d extr1 = frac.VerticesOfFractures[i][k];
                    Vector3d extr2 = frac.VerticesOfFractures[i][(k + 1) % frac.NumberOfVertices[i]];

                    MatrixXd Ainters(3, 2);
                    Ainters.col(0) = (extr2 - extr1);
                    Ainters.col(1) = t;

                    // Se si verifica questa condizione, ovvero che il prodotto vettoriale tra il vettore direzionale della frattura 
                    // e il vettore tangente t è prossimo a 0 (quindi il vettore direzionale e t sono paralleli),
                    // vuol dire che non viene creata un'intersezione
                    if ((extr2 - extr1).cross(t).squaredNorm() < machineEpsilon)
                    {
                        continue;
                    }

                    // Rappresenta la differenza tra il punto di intersezione extrTrace e uno degli estremi extr1
                    Vector3d bInters = (extrTrace - extr1);

                    Vector2d xSolution = solveWithDecompositionQR(Ainters, bInters);
                    // Il valore di alfa deve appartenere all'intervallo [0 - defTol, 1 + defTol]
                    if (xSolution[0] > (0 - machineEpsilon) && xSolution[0] < (1 + machineEpsilon) && NumOfExtremes < 2)
                    {
                        // Vettore 3d che calcola il punto di intersezione
                        Vector3d pointIntersec = extr1 + (xSolution[0] * (extr2 - extr1));
                        // Verifica se il punto pointIntersec, di intersezione, si trova all'interno della frattura j
                        if (NumOfExtremes == 0 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[j]))
                        {
                            // Si verifica l'aggiunta del punto di intersezione a extremes
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                        }
                        // Verifica se il punto di intersezione è diverso dal primo e dentro la frattura j
                        // NumOfExtremes == 1: verifica se è già stato trovato il primo estremo
                        // extremes[0] != pointIntersec: verifica se il punto di intersezione non è lo stesso del primo estremo
                        else if (NumOfExtremes == 1 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[j]) && extremes[0] != pointIntersec)
                        {
                            // In tal modo, viene registrata una nuova traccia che collega le fratture coinvolte i e j
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                            trac.NumberOfTraces += 1;
                            trac.TraceId.push_back(count);
                            Vector2i idFracturesTrace = {frac.FractureId[i], frac.FractureId[j]};
                            trac.FractureIds.push_back(idFracturesTrace);
                            trac.VerticesOfTraces[count] = extremes;
                            count += 1;
                        }
                    }
                }

                // Ciclo FOR che serve per eseguire i controlli sulla seconda frattura. E' un ciclo che itera sui vertici
                // della seconda frattura al fine di trovare gli estremi che definiscono la traccia di intersezione
                for (unsigned int k = 0; k < frac.NumberOfVertices[j]; k++)
                {
                    // Con extr1 ed extr2 vengono indicati 2 estremi consecutivi della frattura
                    Vector3d extr1 = frac.VerticesOfFractures[j][k];
                    Vector3d extr2 = frac.VerticesOfFractures[j][(k + 1) % frac.NumberOfVertices[j]];

                    MatrixXd Ainters(3, 2);
                    Ainters.col(0) = (extr2 - extr1);
                    Ainters.col(1) = t;

                    // Se si verifica questa condizione, ovvero che il prodotto vettoriale tra il vettore direzionale della frattura 
                    // e il vettore tangente t è prossimo a 0 (quindi il vettore direzionale e t sono paralleli),
                    // vuol dire che non viene creata un'intersezione
                    if ((extr2 - extr1).cross(t).squaredNorm() < machineEpsilon)
                    {
                        continue;
                    }

                    // Rappresenta la differenza tra il punto di intersezione extrTrace e uno degli estremi extr1
                    Vector3d bInters = (extrTrace - extr1);

                    Vector2d xSolution = solveWithDecompositionQR(Ainters, bInters);
                    // Il valore di alfa deve appartenere all'intervallo [0 - defTol, 1 + defTol]
                    if (xSolution[0] > (0 - machineEpsilon) && xSolution[0] < (1 + machineEpsilon) && NumOfExtremes < 2)
                    {
                        // Vettore 3d che calcola il punto di intersezione
                        Vector3d pointIntersec = extr1 + (xSolution[0] * (extr2 - extr1));
                        // Verifica se il punto pointIntersec, di intersezione, si trova all'interno della frattura j
                        if (NumOfExtremes == 0 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[i]))
                        {
                            // Si verifica l'aggiunta del punto di intersezione a extremes
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                        }
                        // Verifica se il punto di intersezione è diverso dal primo e dentro la frattura j
                        // NumOfExtremes == 1: verifica se è già stato trovato il primo estremo
                        // extremes[0] != pointIntersec: verifica se il punto di intersezione non è lo stesso del primo estremo
                        else if (NumOfExtremes == 1 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[i]) && extremes[0] != pointIntersec)
                        {
                            // In tal modo, viene registrata una nuova traccia che collega le fratture coinvolte i e j
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                            trac.NumberOfTraces += 1;
                            trac.TraceId.push_back(count);
                            Vector2i idFracturesTrace = {frac.FractureId[i], frac.FractureId[j]};
                            trac.FractureIds.push_back(idFracturesTrace);
                            trac.VerticesOfTraces[count] = extremes;
                            count += 1;
                        }
                    }
                }
            }
        }
    }
}



// Definizione di una funzione che verifica se un punto si trova su un segmento
bool isPointInSegment(const Vector3d& point, const Vector3d& p1, const Vector3d& p2)
{
    Vector3d v1 = point - p1;
    Vector3d v2 = p2 - p1;
    double crossProduct = v1.cross(v2).norm();

    // Se si verifica tale condizione, significa che il punto non è esattamente sul segmento
    if (crossProduct > defaultTolerance)
    {
        return false;
    }

    double dotProduct = v1.dot(v2);

    // Se si verifica almeno una delle 2 condizioni, allora il punto non appartiene al segmento
    // Prima condizione: il prodotto scalare è negativo, quindi il punto è prima di p1
    // Seconda condizione: il prodotto scalare è maggiore della distanza tra p1 e p2, quindi il punto è dopo p2
    if (dotProduct < 0 || dotProduct > v2.squaredNorm())
    {
        return false;
    }

    return true;
}

// Definizione di una funzione che calcola la tipologia della traccia, passante o non passante
void computeTypeTrace(const strFractures& fractures, strTraces& traces)
{
    // Ciclo FOR principale, che itera su tutte le tracce presenti in traces
    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        // Estrae le estremità iniziale e finale della traccia corrente
        const Vector3d& traceStart = traces.VerticesOfTraces.at(i)[0];
        const Vector3d& traceEnd = traces.VerticesOfTraces.at(i)[1];

        // Variabile booleana che indica se la traccia è passante (false) o non passante (true)
        bool Tips = true;  // Si assume inizialmente che la traccia sia non passante

        // Ciclo FOR che itera su tutte le frattura
        for (unsigned int k = 0; k < fractures.NumberOfFractures; ++k)
        {
            // Variabile per verificare se gli estremi sono su lati diversi della frattura
            bool foundDifferentSides = false;
            bool endsOnDifferentSides = false;

            // Itera su tutti i lati della frattura corrente, tenendo in considerazione i vertici
            for (unsigned int l = 0; l < fractures.NumberOfVertices.at(k); ++l)
            {
                // Estrae i vertici di un lato della frattura
                Vector3d p1F = fractures.VerticesOfFractures.at(k).at(l);
                Vector3d p2F = fractures.VerticesOfFractures.at(k).at((l + 1) % fractures.NumberOfVertices.at(k));

                // Verifica se traceStart si trova su questo lato della frattura
                if (isPointInSegment(traceStart, p1F, p2F))
                {
                    // Indicatore che serve per iterare sui lati successivi fino alla fine dei vertici della frattura
                    unsigned int m = (l + 1) % fractures.NumberOfVertices.at(k);

                    while (m != l)
                    {
                        // Estrae i vertici del lato successivo della frattura
                        Vector3d next_p1F = fractures.VerticesOfFractures.at(k).at(m);
                        Vector3d next_p2F = fractures.VerticesOfFractures.at(k).at((m + 1) % fractures.NumberOfVertices.at(k));

                        // Verifica se traceEnd si trova sul lato successivo e non sul lato di traceStart
                        if (isPointInSegment(traceEnd, next_p1F, next_p2F) && !isPointInSegment(traceEnd, p1F, p2F))
                        {
                            endsOnDifferentSides = true;
                            break;
                        }

                        m = (m + 1) % fractures.NumberOfVertices.at(k);  // Viene eseguito il controllo sul lato successivo
                    }

                    // Se entrambi gli estremi si trovano su lati diversi della frattura, la traccia è passante
                    if (endsOnDifferentSides)
                    {
                        foundDifferentSides = true;
                        break;
                    }
                }
            }

            // Se i due estremi della traccia sono su lati diversi della frattura, la traccia è passante
            if (foundDifferentSides)
            {
                Tips = false;
                break;
            }
        }

        // Il valore di Tips ottenuto viene aggiunto alla lista
        traces.Tips.push_back(Tips);
    }
}



// Definizione di una funzione per calcolare la lunghezza della traccia
void lengthTraces(strTraces& traces)
{
    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        const Vector3d& traceStart = traces.VerticesOfTraces[i][0];
        const Vector3d& traceEnd = traces.VerticesOfTraces[i][1];

        // Dati gli estremi della traccia, viene calcolata la loro lunghezza
        double length = distance(traceEnd, traceStart);
        traces.TraceLength.push_back(length);
    }
}



// Definizione di una funzione per ordinare le tracce secondo l'algoritmo Merge Sort
void orderTraces(strTraces& traces)
{
    vector<pair<double, unsigned int>> passing;
    vector<pair<double, unsigned int>> notPassing;

    // Popola vettori passante e nonPassante
    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        // Se il valore i di Tips + False, la traccia è passante. Altrimenti la traccia è non passante
        if (traces.Tips[i] == false)
        {
            passing.push_back({traces.TraceLength[i], i });
        }
        else
        {
            notPassing.push_back({traces.TraceLength[i], i });
        }
    }

    // Ordina vettori per lunghezza decrescente
    SortLibrary::MergeSort(passing);
    SortLibrary::MergeSort(notPassing);

    // Strutture che servono per ricostruire le tracce ordinate secondo Merge Sort
    vector<unsigned int> orderedTraceId;
    vector<Vector2i> orderedFractureIds;
    map<unsigned int, array<Vector3d, 2>> orderedVerticesOfTraces;
    vector<bool> orderedTips;
    vector<double> orderedLengths;


    // Inserimento delle strutture nella tipologia Passante
    for (const auto& trace : passing)
    {
        unsigned int idx = trace.second;
        orderedFractureIds.push_back(traces.FractureIds[idx]);
        orderedVerticesOfTraces[idx] = traces.VerticesOfTraces[idx];
        orderedTips.push_back(traces.Tips[idx]);
        orderedLengths.push_back(traces.TraceLength[idx]);
        orderedTraceId.push_back(traces.TraceId[idx]);
    }

    // Inserimento delle strutture nella tipologia Non Passante
    for (const auto& trace : notPassing)
    {
        unsigned int idx = trace.second;
        orderedFractureIds.push_back(traces.FractureIds[idx]);
        orderedVerticesOfTraces[idx] = traces.VerticesOfTraces[idx];
        orderedTips.push_back(traces.Tips[idx]);
        orderedLengths.push_back(traces.TraceLength[idx]);
        orderedTraceId.push_back(traces.TraceId[idx]);
    }

    // Aggiornamento della struttura con i vettori ordinati
    traces.FractureIds = orderedFractureIds;

    // Pulizia e aggiornamento di VerticesOfTraces
    traces.VerticesOfTraces.clear();
    for (const auto& pair : orderedVerticesOfTraces)
    {
        traces.VerticesOfTraces[pair.first] = pair.second;
    }

    traces.TraceId = orderedTraceId;
    traces.Tips = orderedTips;
    traces.TraceLength = orderedLengths;
}



