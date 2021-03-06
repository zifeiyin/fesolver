CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

SET(CMAKE_INSTALL_PREFIX  "/Users/zyin/Github/fesolver/platform" )

FIND_PACKAGE(deal.II)
IF(NOT ${deal.II_FOUND})
    MESSAGE(  
        FATAL_ERROR "\n"
        "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    )
ENDIF()

IF(NOT DEAL_II_WITH_PETSC OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
    MESSAGE(FATAL_ERROR "\n"
        "Error! This tutorial requires a deal.II library that was configured with the following options:
            DEAL_II_WITH_PETSC = ON
            DEAL_II_PETSC_WITH_COMPLEX = OFF
        However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
            DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
            DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
        which conflict with the requirements."
    )
ENDIF()

IF(NOT DEAL_II_WITH_UMFPACK) # keep in one line
  MESSAGE(FATAL_ERROR "Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_UMFPACK = ON
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_UMFPACK = ${DEAL_II_WITH_UMFPACK}
which conflict with the requirements."
    )
ENDIF()

# compiler flags
set(CMAKE_CXX_FLAGS_DEBUG " -g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE " -O3 -Wall")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(finiteElementSolver)

ADD_EXECUTABLE( FEsolver 
                src/parameters/parameters.cpp
                src/parameters/read_parameters.cpp
                src/parameters/declare_parameters.cpp
                src/solver/solver.cpp
                src/mesh/mesh.cpp
                src/main/main.cpp )

DEAL_II_SETUP_TARGET(FEsolver)

install ( TARGETS FEsolver DESTINATION bin )
