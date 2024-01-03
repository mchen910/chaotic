set(CHAOTIC_SOURCE_FILES    
    solver/dataframe.cpp           
    solver/ODE.cpp                 
    solver/solver.cpp              
)

set(CHAOTIC_HEADER_FILES    
    ../include/dataframe.h             
    ../include/ODE.h                   
    ../include/solver.h           
)

add_library(chaotic STATIC ${CHAOTIC_SOURCE_FILES} ${CHAOTIC_HEADER_FILES})
target_include_directories(chaotic PUBLIC ../include)