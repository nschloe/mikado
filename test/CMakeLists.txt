FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(${Mikado_SOURCE_DIR}/src/)
INCLUDE_DIRECTORIES(SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  )

SET(internal_LIBS "mikado")

SET(HELPERSTEST_EXECUTABLE "helpersTest")
ADD_EXECUTABLE(${HELPERSTEST_EXECUTABLE} helpers.cpp main.cpp)
ADD_TEST(helpersTestSerial ${HELPERSTEST_EXECUTABLE})
TARGET_LINK_LIBRARIES(${HELPERSTEST_EXECUTABLE} ${internal_LIBS})

SET(LINEARTEST_EXECUTABLE "linearTest")
ADD_EXECUTABLE (${LINEARTEST_EXECUTABLE} linear.cpp main.cpp)
# add tests
ADD_TEST(linearTestSerial ${LINEARTEST_EXECUTABLE})
# allow root for docker
ADD_TEST(linearTestMpi3
  ${Trilinos_MPI_EXEC} --allow-run-as-root -oversubscribe --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 3 ${LINEARTEST_EXECUTABLE}
  )
TARGET_LINK_LIBRARIES(${LINEARTEST_EXECUTABLE} ${internal_LIBS})

# copy over yaml file
configure_file(nonlinear-config.yml nonlinear-config.yml COPYONLY)

SET(NONLINEARTEST_EXECUTABLE "nonlinearTest")
ADD_EXECUTABLE(${NONLINEARTEST_EXECUTABLE} nonlinear.cpp main.cpp)
ADD_TEST(nonlinearTestSerial ${NONLINEARTEST_EXECUTABLE})
TARGET_LINK_LIBRARIES(${NONLINEARTEST_EXECUTABLE} ${internal_LIBS})
