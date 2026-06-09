cmake_minimum_required(VERSION 3.18)
include(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

foreach(t sequence matching align)
  add_executable(test_${t} test_${t}.cpp)
  target_link_libraries(test_${t} PRIVATE g2g_lib GTest::gtest_main)
  gtest_discover_tests(test_${t})
endforeach()

add_executable(test_decompose test_decompose.cpp)
target_link_libraries(test_decompose PRIVATE g2g_lib GTest::gtest_main)
gtest_discover_tests(test_decompose)
