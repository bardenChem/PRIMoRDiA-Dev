cmake_minimum_required(VERSION 3.16)

project(primordia_app LANGUAGES CXX)

# Gera compile_commands.json para ferramentas como IntelliSense.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(
    GLOB SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(PRIMoRDiA_1.50v ${SOURCES})

# O projeto requer C++17 por causa de std::filesystem.
target_compile_features(
    PRIMoRDiA_1.50v
    PRIVATE
    cxx_std_17
)

target_include_directories(
    PRIMoRDiA_1.50v
    PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams)

target_link_libraries(
    PRIMoRDiA_1.50v
    PRIVATE
    Eigen3::Eigen
    OpenMP::OpenMP_CXX
    Boost::iostreams
    stdc++fs
)

# Facilita inspecionar variáveis durante o Debug.
target_compile_options(
    PRIMoRDiA_1.50v
    PRIVATE
    $<$<CONFIG:Debug>:-O0>
)

# GCC antigo precisava dessa biblioteca para std::filesystem.
if(
    CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9
)
    target_link_libraries(
        PRIMoRDiA_1.50v
        PRIVATE
        stdc++fs
    )
endif()