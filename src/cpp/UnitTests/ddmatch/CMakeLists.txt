project(UnitTesting_ddmatch)
set(src
  main.cpp
)

add_executable(${PROJECT_NAME} ${src})

add_dependencies(${PROJECT_NAME}
  ddmatchLib
)

if (UNIX)
target_link_libraries(${PROJECT_NAME}
  ddmatchLib
  gtest
  pthread
)
endif (UNIX)

if (MSVC)
target_link_libraries(${PROJECT_NAME}
  ddmatchLib
  gtest
)
endif(MSVC)
