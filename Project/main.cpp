#include "Utils.hpp"
#include "importExport.hpp"
#include "FracturesAndTraces.hpp"
#include "SortingAlgorithm_MERGESORT.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <Eigen/Eigen>
#include <vector>
#include <cmath>
#include <algorithm>


using namespace std;
using namespace Eigen;
using namespace DiscreteAndFractureNetworkLibrary;
using namespace SortLibrary;


int main()
{
    // *************************************************************************************************************************************
    // PRIMA PARTE DEL PROGETTO ************************************************************************************************************
    // *************************************************************************************************************************************


    // Dichiarazione delle strutture utilizzate
    strFractures fractures;
    strTraces traces;


    // IMPORT del file contenente le infomazioni sulle fratture
    string fileInput = "./DFN/FR362_data.txt";

    if(!importListFractures(fractures, fileInput))
    {
        cerr << "ATTENTION: Error importing the file " << fileInput << endl;
        return 1;
    }
    else
    {
        cout << "CHECK: fractures are imported correctly!" << endl;
    }

    // Verifica se i 2 test fondamentali per le fratture sono stati eseguiti
    bool test1Passed = testEdgesOfFracture(fractures);
    bool test2Passed = testVerticesOfFracture(fractures);

    // Se almeno uno dei 2 test viene eseguito con esito NEGATIVO, allora viene restituito errore
    if (!test1Passed || !test2Passed)
    {
        cerr << "ATTENTION: Tests failed. This is an error!" << endl;
        return 1;
    }
    else
    {
        cout << "CHECK: tests on edges and vertices are correct!" << endl;
    }


    definitionOfTraces(fractures, traces);


    // EXPORT del file contenente le informazioni sulle tracce
    string outputTraceInfo = "./Traces_Info_FR362.txt";

    if(!generateTracesInfo(traces, outputTraceInfo))
    {
        return 3;
    }


    // EXPORT del file contenente le tipologie delle tracce
    string outputTraceTips = "./Traces_Tips_FR362.txt";

    computeTypeTrace(fractures, traces);
    lengthTraces(traces);
    orderTraces(traces);


    if (!generateTracesTips(fractures, traces, outputTraceTips))
    {
        return 4;
    }

    return 0;
}




