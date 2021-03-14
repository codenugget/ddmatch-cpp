project(example)

set(src
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
