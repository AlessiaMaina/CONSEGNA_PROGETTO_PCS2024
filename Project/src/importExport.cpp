#include "FracturesAndTraces.hpp"
#include "importExport.hpp"
#include "Utils.hpp"
#include "SortingAlgorithm_MERGESORT.hpp"


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <Eigen/Eigen>

#include <set>


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
    cout << "Il NUMERO DI TRACCE e' pari a: " << trac.NumberOfTraces << endl;

    // Ciclo FOR che, per ogni traccia, stampa il relativo identificativo e le coordinate degli estremi che la delimitano
    for (unsigned int i = 0; i < trac.NumberOfTraces; i++)
    {
        outFile << trac.TraceId[i] << "; "
             << trac.FractureIds[i][0] << "; " << trac.FractureIds[i][1] << "; "
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



bool generateTracesTips(const strFractures& fractures, strTraces& traces, const string& outputFilePath) {
    ofstream outFile(outputFilePath);

    if (!outFile.is_open())
    {
        return false;
    }

    // Il valore di setprecision è pari a quello utilizzato nei valori in input delle coordinate dei vertici
    outFile << setprecision(16) << scientific;

    // Variabile per memorizzare le tracce che sono prima NON passanti e poi passanti
    vector<pair<double, unsigned int>> tracesCaseTF;   // lunghezza -> ID della traccia


    // Iterazione su ogni frattura
    for (unsigned int k = 0; k < fractures.NumberOfFractures; ++k)
    {
        unsigned int numTraces = 0;
        set<unsigned int> countTraces;  // Per evitare il doppio conteggio

        // Conta il numero di tracce associate alla frattura corrente
        for (const auto& trace : traces.FractureIds)
        {
            // Verifica se la frattura corrente (indicata con k) è una delle 2 fratture associate alla traccia
            if (trace[0] == k || trace[1] == k)
            {
                // Se si verifica la 1a condizione (trace[0] == k), allora viene aggiunto trace[1]
                // Se si verifica la 2a condizione (trace[1] == k), allora viene aggiunto trace[0]
                countTraces.insert(trace[0] == k ? trace[1] : trace[0]);
            }
        }

        // In tal modo countTraces contiene le fratture che sono associate alla frattura k
        // Calcolandone la dimensione, si ottiene il numero totale di tracce totali
        // create dalla frattura k con le altre fratture
        numTraces = countTraces.size();

        if (numTraces == 0)
        {
            // Se non ci sono tracce associate alla frattura corrente k, passa alla successiva iterazione del ciclo
            continue;
        }

        outFile << "# FractureId; NumTraces" << endl;
        outFile << fractures.FractureId[k] << "; " << numTraces << endl;
        outFile << "# TraceId; Tips; Length" << endl;

        // Set utilizzato per tenere traccia delle tracce già stampate
        set<unsigned int> printedTraces;

        // Verifica che tracesCaseFT non sia vuoto
        if (!tracesCaseTF.empty())
        {
            // Esecuzione dell'ordinamento di tracesCaseFT, tramite l'algoritmo Merge Sort
            SortLibrary::MergeSort(tracesCaseTF);
            for (const auto& trac : tracesCaseTF)
            {
                // Estrazione dell'id della traccia dalla coppia trac e memorizzazione nella variabile traceId
                unsigned int traceId = trac.second;
                // Estrazione della lunghezza della traccia dalla coppia trac e memorizzazione nella variabile length
                double length = trac.first;

                // Verifica se l'ID della traccia non è già stato stampato
                if (printedTraces.find(traceId) == printedTraces.end())
                {
                    for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
                    {
                        // Verifica se l'id della traccia corrente corrisponde a quello di traceId
                        // e se la traccia è associata alla frattura corrente k
                        if (traces.TraceId[i] == traceId && (traces.FractureIds[i][0] == k || traces.FractureIds[i][1] == k)) {

                            // Se il secondo valore di Tips associato alla traccia è false,
                            // stampa la traccia e inserisce il suo id nelle tracce stampate, ovvero in printedTraces
                            if (traces.Tips[i][1] == false)
                            {
                                outFile << traceId << "; false; " << length << endl;
                                printedTraces.insert(traceId);
                                break;
                            }
                        }
                    }
                }
            }
        }

        // Ciclo FOR che stampa le tracce passanti associate alla frattura corrente

        // In questo caso, nel ciclo FOR viene usato un iteratore per scorrere il vector
        // In tal modo, si possono rimuovere gli elementi dal vector durante l'iterazione senza
        // invalidare l'iteratore. Se si fosse usato "for (const auto& trac : traces.passing)",
        // non si sarebbe potuto modificare direttamente il vector
        for (auto it = traces.passing.begin(); it != traces.passing.end(); )
        {
            // Estrazione dell'id della traccia dalla coppia it e memorizzazione nella variabile traceId
            unsigned int traceId = it->second;
            // Estrazione della lunghezza della traccia dalla coppia it e memorizzazione nella variabile length
            double length = it->first;

            // Variabile booleana che serve a tener conto se la traccia è stata stampata o meno
            bool tracePrinted = false;
            // Variabile booleana (inizializzata a true) per verificare se Tips è {false, false}
            bool caseFF = true;

            for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
            {
                // Verifica se l'id della traccia corrente corrisponde a quello di traceId
                // e se la traccia è associata alla frattura corrente k
                if (traces.TraceId[i] == traceId && (traces.FractureIds[i][0] == k || traces.FractureIds[i][1] == k)) {

                    // Se si verifica se Tips è {false, false}, allora caseFF assume valore true
                    if (traces.Tips[i][0] == false && traces.Tips[i][1] == false)
                    {
                        caseFF = true;
                    }
                    else
                    {
                        caseFF = false;
                    }

                    // Se il primo valore di Tips è true, allora la coppia (length, traceId)
                    // inerente a quella traccia viene inserita all'interno di traceCaseFT
                    if (traces.Tips[i][0] == true)
                    {
                        tracesCaseTF.push_back(make_pair(length, traceId));
                        // Si continua con la prossima traccia senza stampare
                        continue;
                    }

                    outFile << traceId << "; false; " << length << endl;
                    // Stampa l'id della traccia all'interno delle tracce già stampate
                    printedTraces.insert(traceId);
                    // Indica che la traccia è stata stampata
                    tracePrinted = true;
                    break;
                }
            }
            // Verifica se l'id della traccia è stato stampato e Tips della traccia NON è il caso
            // {false, false}
            if (tracePrinted && !caseFF)
            {
                it = traces.passing.erase(it);
                // Rimozione della traccia stampata da passing in modo tale da non prenderla più in considerazione
            }
            else
            {
                ++it;
            }
        }


        // Ciclo FOR che stampa le tracce non passanti associate alla frattura corrente
        for (const auto& trac : traces.notPassing)
        {
            // Estrazione dell'id della traccia dalla coppia trac e memorizzazione nella variabile traceId
            unsigned int traceId = trac.second;
            // Estrazione della lunghezza della traccia dalla coppia trac e memorizzazione nella variabile length
            double length = trac.first;

            // Verifica se l'id della traccia (traceId) è già presente nel set printedTraces
            if (printedTraces.find(traceId) != printedTraces.end())
            {
                continue;
            }

            for (unsigned int i = 0; i < traces.NumberOfTraces; ++i)
            {
                // Verifica se l'id della traccia corrente corrisponde a quello di traceId
                // e se la traccia è associata alla frattura corrente k
                if (traces.TraceId[i] == traceId && (traces.FractureIds[i][0] == k || traces.FractureIds[i][1] == k))
                {
                    outFile << traceId << "; true; " << length << endl;
                    // Stampa l'id della traccia all'interno delle tracce già stampate printedTraces
                    printedTraces.insert(traceId);
                    break;
                }
            }
        }

        outFile << "\n" << endl;
    }

    outFile.close();
    return true;
}







