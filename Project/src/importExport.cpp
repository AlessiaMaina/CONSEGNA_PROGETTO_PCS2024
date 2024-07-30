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
    getline(file, line);        // Legge la riga # Number of Fractures
    getline(file, line);

    // Conversione del contenuto della variabile line (stringa) in un intero senza segno e lo assegna a NumFrac
    unsigned int NumFrac = 0;
    NumFrac = stoi(line);

    unsigned int FracId = 0;
    // Definizione del delimitatore, che in questo caso è ;
    char del = ';';
    unsigned int NumVert = 0;

    // Avviene il ridimensionamento basato su NumFrac per ogni variabile delle fratture
    fractures.NumberOfFractures = NumFrac;
    fractures.FractureId.resize(NumFrac);
    fractures.NumberOfVertices.resize(NumFrac);
    fractures.VerticesOfFractures.resize(NumFrac);

    for(unsigned int i = 0; i < NumFrac; i++)
    {
        getline(file, line);        // Legge le righe # FractureId; NumVertices

        getline(file, line);

        // Consversione dei dati della riga in variabili FracId e NumVert, separati da "del"
        istringstream converter(line);
        converter >> FracId >> del >> NumVert;

        // DEBUG per verificare se gli Id e il numero delle fratture sono stati importati correttamente
        cout << "Fracture ID: " << FracId << ", Number of vertices: " << NumVert << endl;

        // Memorizzazione dell'id e del numero di vertici per la frattura corrente
        fractures.FractureId[i] = FracId;
        fractures.NumberOfVertices[i] = NumVert;
        // Ridimensionamento del vettore di vertici per la frattura corrente, in base a NumVert
        fractures.VerticesOfFractures[i].resize(NumVert);

        getline(file, line);        // Legge le righe # Vertices

        // Ciclo FOR per la lettura delle coordinate 3d dei vertici
        for(unsigned int k = 0; k < 3; k++)
        {
            getline(file, line);
            istringstream converter(line);

            // Ciclo FOR che itera per tutti i vertici
            for(unsigned int j = 0; j < NumVert; j++)
            {
                double coordinates;
                // Estrazione del valore della coordinata
                converter >> coordinates >> del;

                if(k == 0)      // CASO coordinata x
                {
                    fractures.VerticesOfFractures[i][j].x() = coordinates;
                }
                else if (k==1)  // CASO coordinata y
                {
                    fractures.VerticesOfFractures[i][j].y() = coordinates;
                }
                else            // CASO coordinata z
                {
                    fractures.VerticesOfFractures[i][j].z() = coordinates;
                }
            }
        }
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







