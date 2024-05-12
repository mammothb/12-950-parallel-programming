cmake_minimum_required(VERSION 3.25)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Use -std=c++xx instead of -std=gnu++xx
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

project(
  WaveEquation
  VERSION 0.0.1
  LANGUAGES CXX)

include(FetchContent)

add_compile_options(
  -pedantic
  -Wattributes
  -Wextra
  -Wcast-align
  -Wcast-qual
  -Wctor-dtor-privacy
  -Wdisabled-optimization
  -Wformat=2
  -Winit-self
  -Wlogical-op
  -Wmissing-declarations
  -Wmissing-include-dirs
  -Wno-attributes
  -Wnoexcept
  -Wnon-virtual-dtor
  -Wold-style-cast
  -Woverloaded-virtual
  -Wredundant-decls
  -Wshadow
  -Wsign-promo
  -Wstrict-null-sentinel
  -Wstrict-overflow=5
  -Wswitch-default
  -Wundef
  -Wvla
  $<$<CONFIG:Release>:-O3>)

set(WaveEquation_LINK_LIBRARIES "")
# OpenMP
find_package(OpenMP)
list(APPEND WaveEquation_LINK_LIBRARIES OpenMP::OpenMP_CXX)
# fmt
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG 10.2.1
  SYSTEM)
FetchContent_MakeAvailable(fmt)
list(APPEND WaveEquation_LINK_LIBRARIES fmt::fmt)
# Boost
find_package(Boost REQUIRED program_options)
list(APPEND WaveEquation_LINK_LIBRARIES Boost::program_options)

add_subdirectory(src)
