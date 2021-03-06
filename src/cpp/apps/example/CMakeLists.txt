project(example)

set(src
    DiffeoFunctionMatching.cpp
    DiffeoFunctionMatching.h
    Diffeo_functions.cpp
    Diffeo_functions.h
    impl.cpp
    main.cpp
)

add_executable(${PROJECT_NAME} ${src})

add_dependencies(${PROJECT_NAME}
    ddmatchLib
)

target_link_libraries(${PROJECT_NAME}
    ddmatchLib
)
