add_executable(WaveEquation)
target_sources(WaveEquation PRIVATE main.cpp Solver.cpp Timer.cpp)
target_compile_features(WaveEquation PUBLIC cxx_std_20)
target_link_libraries(WaveEquation PUBLIC ${WaveEquation_LINK_LIBRARIES})

add_executable(WaveEquationOMP)
target_sources(WaveEquationOMP PRIVATE main_omp.cpp SolverOmp.cpp Timer.cpp)
target_compile_features(WaveEquationOMP PUBLIC cxx_std_20)
target_link_libraries(WaveEquationOMP PUBLIC ${WaveEquation_LINK_LIBRARIES})
