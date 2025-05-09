cmake_minimum_required(VERSION 3.17)
project(pmt_simulation VERSION 0.1)

find_package(OpenCV CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(pmt_simulation 
    main.cxx 
    src/pmt_hits.cxx 
    src/signal_generation.cxx
)

target_include_directories(pmt_simulation PUBLIC
    ${OpenCV_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(pmt_simulation
    ${OpenCV_LIBS}
    opencv_core
    opencv_imgcodecs
)