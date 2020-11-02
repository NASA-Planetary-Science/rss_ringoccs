rm -f *.so
rm -f *.o

CUSTFLAGS=\
"-DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes"\
" -I../include/ -I../special_functions/ -c"

gcc ${CUSTFLAGS} _diffraction_correction.c
gcc ${CUSTFLAGS} _fraunhofer_diffraction.c
gcc ${CUSTFLAGS} _fresnel_diffraction_gap.c
gcc ${CUSTFLAGS} _fresnel_diffraction_left_straightedge.c
gcc ${CUSTFLAGS} _fresnel_diffraction_right_straightedge.c
gcc ${CUSTFLAGS} _fresnel_diffraction_ringlet.c
gcc ${CUSTFLAGS} _fresnel_diffraction_ringlet_phase.c
gcc ${CUSTFLAGS} _fresnel_diffraction_square_wave.c
gcc ${CUSTFLAGS} _fresnel_kernel.c
gcc ${CUSTFLAGS} _fresnel_scale.c
gcc ${CUSTFLAGS} _fresnel_transform.c
gcc ${CUSTFLAGS} _fresnel_transform_ellipse.c
gcc ${CUSTFLAGS} _fresnel_transform_interpolate.c
gcc ${CUSTFLAGS} _fresnel_transform_legendre_even.c
gcc ${CUSTFLAGS} _fresnel_transform_legendre_odd.c
gcc ${CUSTFLAGS} _fresnel_transform_newton.c
gcc ${CUSTFLAGS} _fresnel_transform_perturbednewton.c

gcc -shared \
_diffraction_correction.o \
_fraunhofer_diffraction.o \
_fresnel_diffraction_gap.o \
_fresnel_diffraction_left_straightedge.o \
_fresnel_diffraction_right_straightedge.o \
_fresnel_diffraction_ringlet.o \
_fresnel_diffraction_ringlet_phase.o \
_fresnel_diffraction_square_wave.o \
_fresnel_kernel.o \
_fresnel_scale.o \
_fresnel_transform.o \
_fresnel_transform_ellipse.o \
_fresnel_transform_interpolate.o \
_fresnel_transform_legendre_even.o \
_fresnel_transform_legendre_odd.o \
_fresnel_transform_newton.o \
_fresnel_transform_perturbednewton.o \
-L../special_functions/ -lrssringoccsspecialfunctions -lfftw3 \
-o librssringoccsdiffractioncorrection.so

mv -f librssringoccsdiffractioncorrection.so /usr/local/lib/librssringoccsdiffractioncorrection.so

rm -f *.o