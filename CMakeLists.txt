cmake_minimum_required(VERSION 3.16)

project(primordia_app VERSION 1.50 LANGUAGES CXX)

# For single-configuration generators (e.g. Unix Makefiles),
# use Release by default when the user does not specify a build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(
        CMAKE_BUILD_TYPE
        Release
        CACHE STRING
        "Choose the type of build."
        FORCE
    )
endif()

include(GNUInstallDirs)

# Gera compile_commands.json para ferramentas como IntelliSense.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(
    GLOB SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(PRIMoRDiA_1.50v ${SOURCES})

target_compile_definitions(
    PRIMoRDiA_1.50v
    PRIVATE
    PRIMORDIA_VERSION="${PROJECT_VERSION}"
    PRIMORDIA_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
)

set_target_properties(
    PRIMoRDiA_1.50v
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}"
)

# O projeto requer C++17 por causa de std::filesystem.
target_compile_features(
    PRIMoRDiA_1.50v
    PRIVATE
    cxx_std_17
)

file(
    MAKE_DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/build_info.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/build_info.h"
    @ONLY
)

target_include_directories(
    PRIMoRDiA_1.50v
    PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_BINARY_DIR}/generated"

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

install(
    TARGETS PRIMoRDiA_1.50v
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
