project(fft_test)

set(src
    main.cpp
    impl.cpp
)

add_executable(${PROJECT_NAME} ${src})

#add_dependencies(${PROJECT_NAME}
#    ddmatchLib
#)

#target_link_libraries(${PROJECT_NAME}
#    ddmatchLib
#)
