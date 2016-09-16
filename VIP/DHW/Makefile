
all: static_instab_correction.so

static_instab_correction.so:
	f2py -c -m static_instab_correction static_instab_correction.f90

clean:
	@rm static_instab_correction.so
