
set(SOURCE_FILES
    Bessel.c
    Bessel.h
    Geometry.c
    Geometry.h
    Interpolation.c
    Interpolation.h
    MathConstants.h
    Misc.c
    Misc.h
    RootFinding.c
    RootFinding.h
    SIE_Math.h
    Statistics.c
    Statistics.h
    Trig.c
    Trig.h)

add_library(sie_math STATIC ${SOURCE_FILES})