cmake_minimum_required(VERSION 3.0.0)
project(PVHSS VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(./include)
link_directories(./lib)

add_executable(PVHSS main.cpp timetester.cpp PVHSS.cpp CAPSNARK.cpp HSS.cpp helper.cpp)
target_link_libraries(PVHSS ntl m gmp pthread relic)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)