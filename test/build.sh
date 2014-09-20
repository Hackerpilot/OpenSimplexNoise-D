dmd main.d\
	../src/open_simplex_noise.d\
	../../Derelict3/import/derelict/sdl2/*.d\
	../../Derelict3/import/derelict/util/*.d\
	-I../src/\
	-I../../Derelict3/import/\
	-L-lSDL2\
	-L-ldl
