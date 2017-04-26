set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(SOURCE_FILES
    main.cpp)

add_executable(qm ${SOURCE_FILES})

add_executable(grid_test grid_test.cpp grid.hpp)

add_executable(common_test common_test.cpp)

add_executable(test test.cpp)
