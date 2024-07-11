#include "FracturesAndTraces.hpp"
#include "importExport.hpp"
#include "Utils.hpp"


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;
using namespace DiscreteAndFractureNetworkLibrary;


bool importListFractures(strFractures& fractures, const string& inputFilePath)
{
    ifstream file;
    file.open(inputFilePath);
    if(file.fail())
    {
        cerr << "ATTENTION: IMPORT of data from file " << inputFilePath << " goes wrong!" << endl;
        return false;
    }

    string line;

    getline(file, line);                // Legge la riga # Number of Fractures
    getline(file, line);

    // Conversione della riga in un elemento "unsigned int"
    istringstream convIntoNumFrac(line);
    unsigned int NumOfFractures;
    convIntoNumFrac >> NumOfFractures;

    // Avviene il ridimensionamento di FracturesId e NumberOfVertices
    fractures.NumberOfFractures = NumOfFractures;
    fractures.FractureId.resize(fractures.NumberOfFractures);
    fractures.NumberOfVertices.resize(fractures.NumberOfFractures);


    for(unsigned int i = 0; i < NumOfFractures; i++)
    {
        getline(file, line);            // Legge la riga # FractureId; NumVertices

        getline(file, line);
        istringstream convert(line);

        char del1;
        unsigned int FracId;
        unsigned int numOfVertices;

        convert >> FracId >> del1 >> numOfVertices;
        getline(file, line);

        fractures.FractureId[i] = FracId;
        fractures.NumberOfVertices[i] = numOfVertices;

        // Debug di output per verificare se gli Id e il numero delle fratture sono stati importati correttamente
        cerr << "Fracture ID: " << FracId << ", Number of vertices: " << numOfVertices << endl;

        // Creazione di un vettore di coordinate 3d per i vertici della frattura
        vector<Vector3d> FracVert(numOfVertices);

        // Estrazione dei valori delle coordinate e inserimento all'interno di un Vector3d
        getline(file, line);
        istringstream firstCoordinate(line);
        getline(file, line);
        istringstream secondCoordinate(line);
        getline(file, line);
        istringstream thirdCoordinate(line);

        double xCoord, yCoord, zCoord;
        char del2;

        // Ciclo FOR che, per ogni vertice, stampa le relative coordinate in 3d
        for (unsigned int k = 0; k < numOfVertices; k++)
        {
            firstCoordinate  >> xCoord >> del2;
            secondCoordinate >> yCoord >> del2;
            thirdCoordinate  >> zCoord >> del2;

            FracVert[k] = Vector3d(xCoord, yCoord, zCoord);
        }

        fractures.VerticesOfFractures[FracId] = FracVert;
    }

    file.close();
    return true;
}



// E' il file indicato dal README.pdf nel seguente modo:
// # Number of Traces
// ...(value)
// # TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2
// ...(values)



bool generateTracesInfo(strTraces& trac, const string& outputFilePath)
{
    ofstream outFile;
    outFile.open(outputFilePath);

    if (!outFile.is_open())
    {
        return false;
    }

    outFile << "# Number of Traces" << endl;

    // Il valore di setprecision è pari a quello utilizzato nei valori in input delle coordinate dei vertici
    outFile << setprecision(16) << scientific;

    outFile << trac.NumberOfTraces << endl;
    outFile << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

    // Ciclo FOR che, per ogni traccia, stampa il relativo identificativo e le coordinate degli estremi che la delimitano
    for (unsigned int i = 0; i < trac.NumberOfTraces; i++)
    {
        outFile << i
             << "; " << trac.FractureIds[i][0] << "; " << trac.FractureIds[i][1] << "; "
             << trac.VerticesOfTraces[i][0][0] << "; " << trac.VerticesOfTraces[i][0][1] << "; " << trac.VerticesOfTraces[i][0][2] << "; "
             << trac.VerticesOfTraces[i][1][0] << "; " << trac.VerticesOfTraces[i][1][1] << "; " << trac.VerticesOfTraces[i][1][2] << endl;
    }

    outFile.close();
    return true;
}



// E' il file indicato dal README.pdf nel seguente modo:
// # FractureId; NumTraces
// ...(values)
// # TraceId; Tips; Length
// ...(values)


bool generateTracesTips(const strFractures& fractures, strTraces& traces, const string& outputFilePath)
{
    ofstream outFile;
    outFile.open(outputFilePath);

    if (!outFile.is_open())
    {
        return false;
    }

    // Il valore di setprecision è pari a quello utilizzato nei valori in input delle coordinate dei vertici
    outFile << setprecision(16) << scientific;

    for (unsigned int k = 0; k < fractures.NumberOfFractures; ++k)
    {
        outFile << "# FractureId; NumTraces" << endl;
        unsigned int numTraces = 0;
        for (const auto& trace : traces.FractureIds)
        {
            // Se si verifica almeno una delle 2 condizioni, il numero delle tracce viene incrementato di 1
            // l'operatore || indica l'operatore OR logico
            if (trace.x() == k || trace.y() == k)
            {
                numTraces++;
            }
        }

        outFile << fractures.FractureId[k] << "; " << numTraces << endl;
        outFile << "# TraceId; Tips; Lenght" << endl;

        for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)        
        {
            // Se si verifica almeno una delle 2 condizioni, vengono stampati i dati relativi alla traccia
            if (traces.FractureIds[i].x() == k || traces.FractureIds[i].y() == k)
            {
                outFile << traces.TraceId[i] << "; " << (traces.Tips[i] ? "true" : "false") << "; " << traces.TraceLength[i] << endl;
            }
        }
    }

    outFile.close();
    return true;
}

