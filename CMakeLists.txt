cmake_minimum_required(VERSION 3.18)
project(g2g VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Build type defaults ────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ── Compiler flags ────────────────────────────────────────────
add_compile_options(-Wall -Wextra -Wpedantic)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-O3 -march=native -DNDEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-g -O0 -fsanitize=address,undefined)
  add_link_options(-fsanitize=address,undefined)
endif()

# ── Dependencies ──────────────────────────────────────────────
# libhandlegraph + libsnarls: install via conda or vcpkg
# conda install -c bioconda -c conda-forge vg
# Until then, we stub the interface in include/g2g/snarl_stub.hpp

find_package(Threads REQUIRED)

# Optional: detect if libhandlegraph is available
find_library(HANDLEGRAPH_LIB handlegraph HINTS ENV CONDA_PREFIX PATH_SUFFIXES lib)
find_path(HANDLEGRAPH_INCLUDE handlegraph/handle.hpp HINTS ENV CONDA_PREFIX PATH_SUFFIXES include)

if(HANDLEGRAPH_LIB AND HANDLEGRAPH_INCLUDE)
  message(STATUS "Found libhandlegraph: ${HANDLEGRAPH_LIB}")
  add_compile_definitions(G2G_HAVE_HANDLEGRAPH)
else()
  message(STATUS "libhandlegraph not found — snarl decomposition will use built-in stub")
endif()

# ── Source files ──────────────────────────────────────────────
set(G2G_SOURCES
  src/graph.cpp
  src/sequence.cpp
  src/decompose.cpp
  src/anchor.cpp
  src/matching.cpp
  src/align.cpp
  src/io.cpp
)

add_library(g2g_lib STATIC ${G2G_SOURCES})
target_include_directories(g2g_lib PUBLIC include)

if(HANDLEGRAPH_LIB)
  target_include_directories(g2g_lib PUBLIC ${HANDLEGRAPH_INCLUDE})
  target_link_libraries(g2g_lib PUBLIC ${HANDLEGRAPH_LIB})
endif()

target_link_libraries(g2g_lib PUBLIC Threads::Threads)

# ── Executables ───────────────────────────────────────────────
add_executable(g2g src/main.cpp)
target_link_libraries(g2g PRIVATE g2g_lib)

# ── Tests ─────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
