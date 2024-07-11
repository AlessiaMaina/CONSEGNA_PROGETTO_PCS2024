#pragma once

#include "FracturesAndTraces.hpp"

using namespace std;
using namespace DiscreteAndFractureNetworkLibrary;


bool importListFractures(strFractures& fractures,
                         const string& inputFilePath);

bool generateTracesInfo(strTraces& trac,
                        const string& outputFilePath);

bool generateTracesTips(const strFractures& fractures,
                        strTraces& traces,
                        const string& outputFilePath);

