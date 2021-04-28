cmake_minimum_required( VERSION 3.15 )
project(mc_study)

set( CMAKE_CXX_STANDARD 17 )

find_package( ROOT 6.22 CONFIG REQUIRED )

include_directories(. ../contrib ../utilities)

add_executable( mc_study
        CMakeLists.txt
        mc_study_1.cpp
        mc_study_1.h
        ../utilities/Utilities.hpp
        ../utilities/Histograms.cpp
        ../utilities/Histograms.hpp)


target_link_libraries(mc_study
        PUBLIC ROOT::Core
        ROOT::Physics
        ROOT::Tree
        ROOT::Hist
        ROOT::TreePlayer
        )
