OMP=-openmp
REAL=-fdefault-real-8
DOUBLE=-fdefault-double-8

sfaqdip.x: const.f interp.f sfaqdip.f
	gfortran ${OMP} ${REAL} ${DOUBLE} const.f interp.f sfaqdip.f -o sfaqdip.x

sfaqwave.x: const.f interp.f sfaqwave.f
	gfortran ${OMP} ${REAL} ${DOUBLE} const.f interp.f sfaqwave.f -o sfaqwave.x

spectr.x: const.f spectr.f
	gfortran ${OMP} ${REAL} ${DOUBLE} const.f spectr.f -o spectr.x

timefreq.x: const.f timefreq.f
	gfortran ${OMP} ${REAL} ${DOUBLE} const.f timefreq.f -o timefreq.x

all: sfaqdip.x sfaqwave.x spectr.x timefreq.x
