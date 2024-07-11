list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp")
list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/importExport.cpp")


list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/FracturesAndTraces.hpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/importExport.hpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/SortingAlgorithm_MERGESORT.hpp")


set(DFN_sources ${DFN_sources} PARENT_SCOPE)
set(DFN_headers ${DFN_headers} PARENT_SCOPE)
set(DFN_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
