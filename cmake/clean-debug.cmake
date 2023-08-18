set(debug ${CMAKE_SOURCE_DIR}/bin/flow-simulator.x.dSYM
          ${CMAKE_SOURCE_DIR}/bin/flow-simulator-test.x.dSYM
          ${CMAKE_SOURCE_DIR}/bin/flow-simulator-stress-test.x.dSYM
)

foreach(folder ${debug})
  if (EXISTS ${folder})
     file(REMOVE_RECURSE ${folder})
  endif()
endforeach(folder)
