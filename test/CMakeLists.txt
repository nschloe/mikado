FIND_PACKAGE(MPI REQUIRED)

INCLUDE_DIRECTORIES(${Mikado_SOURCE_DIR}/src/)
INCLUDE_DIRECTORIES(SYSTEM
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  )

SET(internal_LIBS "mikado")

SET(MIKADOTEST_EXECUTABLE "mikadoTest")
ADD_EXECUTABLE (${MIKADOTEST_EXECUTABLE} t.cpp main.cpp)
# add tests
ADD_TEST(mikadoTestSerial ${MIKADOTEST_EXECUTABLE})
# Check out FindMPI.cmake for proper syntax.
# ADD_TEST(ioTestMpi2
#   ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 2 ${IOTEST_EXECUTABLE}
#   )
# ADD_TEST(ioTestMpi7
#   ${Trilinos_MPI_EXEC} --noprefix ${Trilinos_MPI_EXEC_NUMPROCS_FLAG} 7 ${IOTEST_EXECUTABLE}
#   )

# Set executable linking information.
TARGET_LINK_LIBRARIES(
  ${MIKADOTEST_EXECUTABLE}
  ${internal_LIBS}
  )
