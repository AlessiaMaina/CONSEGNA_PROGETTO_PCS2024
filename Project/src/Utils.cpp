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


// Definizione delle tolleranze
const double machineEpsilon = numeric_limits<double>::epsilon();
const double defaultTolerance = pow(10,-10);


double distance(const Vector3d& a, const Vector3d& b)
{
    double distAB = (a - b).norm();
    return distAB;
}



// TEST 1: Verificare nella frattura la presenza di lati tutti diversi da zero
bool testEdgesOfFracture(const strFractures& fractures)
{
    // Itera attraverso ogni frattura nel vettore di fratture
    for (const vector<Vector3d>& fracture : fractures.VerticesOfFractures)
    {
        const size_t vertexCount = fracture.size();

        // Itera attraverso i vertici della frattura
        for (size_t i = 0; i < vertexCount - 1; ++i)
        {
            const Vector3d& v_Start = fracture[i];
            const Vector3d& v_End = fracture[i + 1];

            // Controlla la distanza tra i due vertici è minore o uguale della tolleranza
            if (distance(v_End, v_Start) <= machineEpsilon)
            {
                cerr << "ATTENTION: Fracture has sides of zero length. This is impossible!" << endl;
                return false;
            }
        }
    }

    return true;
}


// TEST 2: Verificare nella frattura che sono presenti almeno 3 vertici, altrimenti la frattura NON è definita
bool testVerticesOfFracture(const strFractures& fractures)
{
    for (const vector<Vector3d>& fracture : fractures.VerticesOfFractures)
    {
        size_t numberOfVertices = fracture.size();

        // Verifica se il numero di vertici è inferiore a 3
        if (numberOfVertices < 3)
        {
            cerr << "ATTENTION: Fracture has less than 3 vertices and is undefined. This is impossible!" << endl;
            return false;
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


// Definizione di una funzione che calcola la Bounding Sphere relativa ad una frattura (un poligono)
// ovvero una sfera circoscritta per un insieme di vertici nello spazio tridimensionale
BoundingSphere computeBoundingSphere(const vector<Vector3d>& vertices)
{
    if (vertices.empty())
    {
        // Se il vettore di coordinate 3D "vertices" è vuoto, viene generata un'eccezione
        throw invalid_argument("ERROR: vector<Vector3d> vertices is empty");
    }

    // Inizializzazione del vettore 3d contenente le coordinate del centro
    Vector3d center = Vector3d::Zero();
    // Inizializzazione del valore double per trovare il raggio massimo della Bounding Sphere
    double radius = 0.0;

    // Viene eseguito il calcolo del centro della Bounding Sphere pari alla media dei vertici
    for (const auto& vertex : vertices)
    {
        center += vertex;
    }

    center /= vertices.size();

    // Iterazione su ogni elemento di vertices
    for (const auto& vertex : vertices)
    {
        // Calcolo della distanza tra il centro e il vertice corrente
        double dist = distance(vertex, center);

        // Se la distanza calcolata è maggiore del raggio attuale, allora il raggio attuale
        // viene aggiornato con questa nuova distanza (fino ad ottenere il raggio massimo)
        if (dist > radius)
        {
            radius = dist;
        }
    }

    return {center, radius};
}


// Definizione di una funzione che verifica se 2 tracce effettivamente si intersecano
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
    // Il prodotto vettoriale deve essere normalizzato, al fine di ottenere il vettore unitario
    // perpendicolare (normale) al poligono

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
    if(abs(sumAreas - areaABC) < defaultTolerance)
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
    // Gli elementi di isPointInTriangle sono: Punto, Primo Vertice, Secondo Vertice, Terzo Vertice
	if (isPointInTriangle(P, fractureVertices[0], fractureVertices[i], fractureVertices[i+1]))
        {
            return true;
        }
    }
    return false;
}


// Definizione della tecnica di decomposizione PA LU
VectorXd solveWithDecompositionPALU(const MatrixXd& matrixA, const VectorXd& columnVector)
{
    VectorXd solutionWithPALU = matrixA.partialPivLu().solve(columnVector);
    return solutionWithPALU;

}

// Definizione della tecnica di decomposizione QR
VectorXd solveWithDecompositionQR(const MatrixXd& matrixA, const VectorXd& columnVector)
{
    VectorXd solutionWithQR = matrixA.colPivHouseholderQr().solve(columnVector);
    return solutionWithQR;

}


// Definizione di una funzione per definire tutti i parametri di una traccia
void definitionOfTraces(strFractures& frac, strTraces& trac)
{
    // Inizializzazione di un contatore a zero, per tenere in considerazione le tracce create
    unsigned int count = 0;

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
                
                // Definizione della normale alla prima frattura
                Vector3d n1 = normal(frac.VerticesOfFractures[i]);
                // Definizione della normale alla seconda frattura
                Vector3d n2 = normal(frac.VerticesOfFractures[j]);
                // Definizione del prodotto vettoriale tra la normale della prima frattura e la normale della seconda
                Vector3d t = n1.cross(n2);
                
                // Matrice QUADRATA dove n1 è contenuto nella 1a riga, n2 nella 2a riga e infine t nella 3a riga
                Matrix3d A;
                A.row(0) = n1.transpose();
                A.row(1) = n2.transpose();
                A.row(2) = t.transpose();

                double b1 = n1.transpose() * boundingSpheres[i].center;
                double b2 = n2.transpose() * boundingSpheres[j].center;
                Vector3d b = {b1, b2, 0};

                // Se il determinante è prossimo a 0, xiò indica che i piani delle fratture
                // sono paralleli (o quasi paralleli) e quindi NON si verifica intersezione

                if (A.determinant() < machineEpsilon)
                {
                    // Tale coppia di fratture NON viene quindi presa in considerazione
                    // e si passa alla coppia di fratture successiva
                    continue;
                }

                // Risoluzione del sistema lineare con la decomposizione PA LU,
                // in quanto A è una matrice quadrata NON singolare
                // DEFINIZIONE: Una matrice quadrata A è NON singolare se det(A) è diverso da 0
                Vector3d extrTrace = solveWithDecompositionPALU(A,b);

                unsigned int NumOfExtremes = 0;	   // Inizializzazione di un contatore per il numero di estremi trovati
                array<Vector3d, 2> extremes;	   // Definizione di un Array per memorizzare gli estremi della traccia

                // Ciclo FOR per la prima frattura. E' un ciclo che itera sui vertici
                // della prima frattura al fine di trovare gli estremi che definiscono la traccia di intersezione
                for (unsigned int k = 0; k < frac.NumberOfVertices[i]; k++)
                {
                    // Con extr1 ed extr2 vengono indicati 2 estremi consecutivi della frattura i
                    Vector3d extr1 = frac.VerticesOfFractures[i][k];
                    Vector3d extr2 = frac.VerticesOfFractures[i][(k + 1) % frac.NumberOfVertices[i]];

                    // Definizione di una matrice RETTANGOLARE di 3 righe e 2 colonne
                    MatrixXd Ainters(3, 2);

                    // Indica il vettore direzionale del lato corrente della frattura i
                    Ainters.col(0) = (extr2 - extr1);
                    // Indica il vettore tangente, necessario per determinare l'intersezione
                    Ainters.col(1) = t;

                    // Se il prodotto vettoriale tra il vettore direzionale e t è prossimo a zero,
                    // i vettori sono paralleli e non c'è intersezione
                    if ((extr2 - extr1).cross(t).squaredNorm() < machineEpsilon)
                    {
                        // Tale lato NON viene preso in considerazione e si passa al lato successivo
                        continue;
                    }

                    // Indica la differenza tra il punto di intersezione extrTrace ed extr1
                    Vector3d bInters = (extrTrace - extr1);

                    // Risoluzione del sistema lineare con la decomposizione QR,
                    // in quanto Ainters è una matrice RETTANGOLARE SOVRADETERMINATA
                    // Trova i parametri per determinare il punto di intersezione sul lato della frattura i
                    Vector2d xSolution = solveWithDecompositionQR(Ainters, bInters);

                    // Il valore di alfa deve appartenere all'intervallo [0 - defTol, 1 + defTol]
                    if (xSolution[0] > (0 - machineEpsilon) && xSolution[0] < (1 + machineEpsilon) && NumOfExtremes < 2)
                    {
                        // Vettore 3d che determina le coordinate esatte del punto di intersezione
                        Vector3d pointIntersec = extr1 + (xSolution[0] * (extr2 - extr1));

                        // Verifica se non è stato trovato alcun estremo e se il punto di intersezione
                        // si trova all'interno della frattura j
                        // In questo modo è assicurato che l'intersezione sia valida per entrambe le fratture
                        if (NumOfExtremes == 0 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[j]))
                        {
                            // Il punto di intersezione viene aggiunto agli estremi (variabile extremes)
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                        }
                        // Verifica se il secondo estremo è valido e diverso dal primo, e se si
                        // trova all'interno della frattura j
                        else if (NumOfExtremes == 1 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[j]) && extremes[0] != pointIntersec)
                        {
                            // Il secondo punto di intersezione viene aggiunto agli estremi, in extremes[1]
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;

                            // Aggiornamento del numero totale di tracce trovate
                            trac.NumberOfTraces += 1;

                            // Aggiunge l'id della nuova traccia a TraceId
                            trac.TraceId.push_back(count);

                            // Creazione di una coppia di id delle fratture che identificano la traccia
                            Vector2i idFracturesTrace = {frac.FractureId[i], frac.FractureId[j]};

                            // Aggiunge la coppia di id a FractureIds
                            trac.FractureIds.push_back(idFracturesTrace);

                            // Aggiunge gli estremi della traccia alla posizione count di VerticesOfTraces
                            trac.VerticesOfTraces.push_back(extremes);

                            count += 1;
                        }
                    }
                }

                // Ciclo FOR per la seconda frattura. E' un ciclo che itera sui vertici
                // della seconda frattura al fine di trovare gli estremi che definiscono la traccia di intersezione
                for (unsigned int k = 0; k < frac.NumberOfVertices[j]; k++)
                {
                    // Con extr1 ed extr2 vengono indicati 2 estremi consecutivi della frattura j
                    Vector3d extr1 = frac.VerticesOfFractures[j][k];
                    Vector3d extr2 = frac.VerticesOfFractures[j][(k + 1) % frac.NumberOfVertices[j]];

                    // Definizione di una matrice RETTANGOLARE di 3 righe e 2 colonne
                    MatrixXd Ainters(3, 2);

                    // Indica il vettore direzionale del lato corrente della frattura j
                    Ainters.col(0) = (extr2 - extr1);
                    // Indica il vettore tangente, necessario per determinare l'intersezione
                    Ainters.col(1) = t;

                    // Se il prodotto vettoriale tra il vettore direzionale e t è prossimo a zero,
                    // i vettori sono paralleli e non c'è intersezione
                    if ((extr2 - extr1).cross(t).squaredNorm() < machineEpsilon)
                    {
                        // Tale lato NON viene preso in considerazione e si passa al lato successivo
                        continue;
                    }

                    // Indica la differenza tra il punto di intersezione extrTrace ed extr1
                    Vector3d bInters = (extrTrace - extr1);

                    // Risoluzione del sistema lineare con la decomposizione QR,
                    // in quanto Ainters è una matrice RETTANGOLARE SOVRADETERMINATA
                    // Trova i parametri per determinare il punto di intersezione sul lato della frattura j
                    Vector2d xSolution = solveWithDecompositionQR(Ainters, bInters);

                    // Il valore di alfa deve appartenere all'intervallo [0 - defTol, 1 + defTol]
                    // In questo modo, il punto di intersezione si trova effettivamente tra extr1 e extr2
                    if (xSolution[0] > (0 - machineEpsilon) && xSolution[0] < (1 + machineEpsilon) && NumOfExtremes < 2)
                    {
                        // Vettore 3d che determina le coordinate esatte del punto di intersezione
                        Vector3d pointIntersec = extr1 + (xSolution[0] * (extr2 - extr1));     

                        // Verifica se non è stato trovato alcun estremo e se il punto di intersezione
                        // si trova all'interno della frattura i
                        // In questo modo è assicurato che l'intersezione sia valida per entrambe le fratture
                        if (NumOfExtremes == 0 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[i]))
                        {
                            // Il primo punto di intersezione viene aggiunto agli estremi, in extremes[0]
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;
                        }
                        // Verifica se il secondo estremo è valido e diverso dal primo, e se si
                        // trova all'interno della frattura i
                        else if (NumOfExtremes == 1 && isPointInFracture(pointIntersec, frac.VerticesOfFractures[i]) && extremes[0] != pointIntersec)
                        {
                            // Il secondo punto di intersezione viene aggiunto agli estremi, in extremes[1]
                            extremes[NumOfExtremes] = pointIntersec;
                            NumOfExtremes += 1;

                            // Aggiornamento del numero totale di tracce trovate
                            trac.NumberOfTraces += 1;

                            // Aggiunge l'id della nuova traccia a TraceId
                            trac.TraceId.push_back(count);

                            // Creazione di una coppia di id delle fratture che identificano la traccia
                            Vector2i idFracturesTrace = {frac.FractureId[i], frac.FractureId[j]};

                            // Aggiunge la coppia di id a FractureIds
                            trac.FractureIds.push_back(idFracturesTrace);

                            // Aggiunge gli estremi della traccia alla posizione count di VerticesOfTraces
                            trac.VerticesOfTraces.push_back(extremes);

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
    Vector3d v2 = point - p1;
    Vector3d v1 = p2 - p1;
    double crossProduct = v1.cross(v2).norm();

    // Se si verifica tale condizione, significa che il punto non è esattamente sul segmento
    if (crossProduct > defaultTolerance)
    {
        return false;
    }

    double dotProduct = v2.dot(v1);
    double squaredLength = v1.dot(v1);

    // Se si verifica almeno una delle 2 condizioni, allora il punto non appartiene al segmento
    // Prima condizione: il prodotto scalare è negativo, quindi il punto è prima di p1
    // Seconda condizione: il prodotto scalare è maggiore della distanza tra p1 e p2, quindi il punto è dopo p2
    if (dotProduct < 0 || dotProduct > squaredLength)
    {
        return false;
    }
    return true;
}


// Definizione di una funzione che calcola la tipologia della traccia, passante o non passante
void computeTypeTrace(const strFractures& fractures, strTraces& traces)
{
    // La dimensione del vettore Tips viene modificata al fine di contenere
    // traces.NumberOfTraces elementi
    traces.Tips.resize(traces.NumberOfTraces);

    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        // Definizione dei 2 estremi della traccia
        const Vector3d& traceStart = traces.VerticesOfTraces.at(i)[0];
        const Vector3d& traceEnd = traces.VerticesOfTraces.at(i)[1];

        // Identificazione delle 2 fratture che formano la traccia
        const Vector2i& fracturePair = traces.FractureIds.at(i);

        // Inizializzazione l'array di booleani per questa traccia
        array<bool, 2> currentTips = {true, true};

        // Controllo per entrambe le fratture nel fracturePair
        for (unsigned int j = 0; j < 2; ++j)
        {
            unsigned int fractureIndex = fracturePair[j];
            // Variabile booleana per determinare se i 2 estremi della traccia si trovano
            // su lati differenti della frattura
            bool endsOnDifferentSides = false;

            // Controllo per ogni lato della frattura
            for (unsigned int l = 0; l < fractures.NumberOfVertices[fractureIndex]; ++l)
            {
                // Definizione dei vertici per il primo lato
                const Vector3d& p1F = fractures.VerticesOfFractures.at(fractureIndex)[l];
                const Vector3d& p2F = fractures.VerticesOfFractures.at(fractureIndex)[(l + 1) % fractures.NumberOfVertices[fractureIndex]];

                if (isPointInSegment(traceStart, p1F, p2F))
                {
                    unsigned int m = (l + 1) % fractures.NumberOfVertices[fractureIndex];

                    // Il ciclo WHILE viene eseguito finchè m non ritorna al punto iniziale l
                    while (m != l)
                    {
                        // Definizione dei vertici per il secondo lato
                        const Vector3d& next_p1F = fractures.VerticesOfFractures.at(fractureIndex)[m];
                        const Vector3d& next_p2F = fractures.VerticesOfFractures.at(fractureIndex)[(m + 1) % fractures.NumberOfVertices[fractureIndex]];

                        // Se il secondo vertice della traccia si trova su un lato differente,
                        // allora la variabile endsOnDifferentSide assume valore "true"
                        if (isPointInSegment(traceEnd, next_p1F, next_p2F))
                        {
                            endsOnDifferentSides = true;
                            break;
                        }

                        // Aggiornamento dell'indice m per passare al prossimo lato del poligono
                        m = (m + 1) % fractures.NumberOfVertices[fractureIndex];
                    }

                    // Se la condizione è vera (endsOnDifferentSides assume valore "true"),
                    // allora l'elemento corrispondente di currentTips viene impostato a false
                    if (endsOnDifferentSides)
                    {
                        currentTips[j] = false;
                        break;
                    }
                }
            }
        }

        // Il risultato viene aggiunto a vector<array<bool,2>> Tips
        traces.Tips[i] = currentTips;
    }
}


// Definizione di una funzione per calcolare la lunghezza della traccia
void lengthTraces(strTraces& traces)
{
    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        // Definizione dei 2 estremi consecutivi della traccia
        const Vector3d& traceStart = traces.VerticesOfTraces[i][0];
        const Vector3d& traceEnd = traces.VerticesOfTraces[i][1];

        // Forniti gli estremi della traccia, viene calcolata la loro lunghezza tramite
        // la funzione "distance" definita precedentemente
        double length = distance(traceEnd, traceStart);

        traces.TraceLength.push_back(length);
    }

}


// Definizione di una funzione per ordinare le tracce
void orderTraces(strTraces& traces)
{
    // Inserimento delle tracce nella categoria passing e notPassing
    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
    {
        const auto& tips = traces.Tips[i];
        const auto& length = traces.TraceLength[i];

        // CASO 1: la traccia è PASSANTE per entrambe le fratture
        if (tips[0] == false && tips[1] == false)
        {
            traces.passing.emplace_back(length, i);
        }
        // CASO 2: la traccia è NON PASSANTE per entrambe le fratture
        else if (tips[0] == true && tips[1] == true)
        {
            traces.notPassing.emplace_back(length, i);
        }
        else
        {
            // CASO 3: la traccia è PASSANTE o NON PASSANTE per la prima frattura
            if (tips[0] == false)
            {
                traces.passing.emplace_back(length, i);
            } else
            {
                traces.notPassing.emplace_back(length, i);
            }
            // CASO 4: la traccia è PASSANTE o NON PASSANTE per la seconda frattura
            if (tips[1] == false)
            {
                traces.passing.emplace_back(length, i);
            } else
            {
                traces.notPassing.emplace_back(length, i);
            }
        }
    }

    // DEBUG: Stampa dei vettori non ordinati
    cout << "NON ordinato - Passing: \n";
    for (const auto& trace : traces.passing)
    {
        cout << "{Length: " << trace.first << ", Idx: " << trace.second << "} \n";
    }
    cout << endl;

    cout << "NON ordinato - NotPassing: \n";
    for (const auto& trace : traces.notPassing)
    {
        cout << "{Length: " << trace.first << ", Idx: " << trace.second << "} \n";
    }
    cout << endl;



    // Ordina i vettori in base alla lunghezza decrescente, solo se i passing e notPassing contengono alemno 1 elemento
    if(traces.passing.size() > 0)
    {
        SortLibrary::MergeSort(traces.passing);
    }
    if(traces.notPassing.size() > 0)
    {
        SortLibrary::MergeSort(traces.notPassing);
    }


    // DEBUG: Stampa dei vettori ordinati
    cout << "Ordinato - Passing: \n";
    for (const auto& trace : traces.passing)
    {
        cout << "{Length: " << trace.first << ", Idx: " << trace.second << "} \n";
    }
    cout << endl;

    cout << "Ordinato - NotPassing: \n";
    for (const auto& trace : traces.notPassing)
    {
        cout << "{Length: " << trace.first << ", Idx: " << trace.second << "} \n";
    }
    cout << endl;
}


