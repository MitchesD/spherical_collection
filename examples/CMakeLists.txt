set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)

add_executable(basic basic.cpp)
target_link_libraries(basic PRIVATE ${PROJECT_NAME})
target_include_directories(basic PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
