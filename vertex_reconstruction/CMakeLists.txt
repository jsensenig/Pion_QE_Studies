cmake_minimum_required( VERSION 3.15 )
project(reco_vertex)

set( CMAKE_CXX_STANDARD 17 )

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities)

add_executable( reco_vtx_study
        CMakeLists.txt
        reco_vtx_study.cpp
        reco_vtx_study.h
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp)


target_link_libraries(reco_vtx_study
        PUBLIC ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
