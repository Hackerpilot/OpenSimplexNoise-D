import std.stdio;
import open_simplex_noise;
import derelict.sdl2.sdl;
import derelict.sdl2.types;
import std.conv:to;

void main()
{
	DerelictSDL2.load();
	SDL_Init(SDL_INIT_VIDEO);
	scope(exit) SDL_Quit();
	SDL_VideoInit(null);
	scope(exit) SDL_VideoQuit();

	auto osn = new OpenSimplexNoise!float;

	enum WIDTH = 512;
	enum HEIGHT = WIDTH;
	enum DEPTH = 3;
	enum BYTES_PER_ROW = DEPTH * WIDTH;

	ubyte[BYTES_PER_ROW * HEIGHT] outputBuffer;

	enum float step = 1.0f / 16.0f;
	for (size_t ix = 0; ix < WIDTH; ix++)
	{
		size_t xx = ix * BYTES_PER_ROW;
		for (size_t iy = 0; iy < HEIGHT; iy++)
		{
			size_t yy = xx + (iy * DEPTH);
			ubyte v = cast(ubyte) ((osn.eval((cast(float) ix) * step, (cast(float) iy) * step) + 1) * (ubyte.max / 2));
			outputBuffer[yy    ] = v;
			outputBuffer[yy + 1] = v;
			outputBuffer[yy + 2] = v;
			static if (DEPTH > 3)
				outputBuffer[yy + 3] = 255;
		}
	}
	SDL_Surface* surface = SDL_CreateRGBSurfaceFrom(outputBuffer.ptr,
		WIDTH, HEIGHT, DEPTH * 8, BYTES_PER_ROW, 0, 0, 0, 0);
//	writeln(to!string(SDL_GetError()));
	assert (surface);
	scope(exit) SDL_FreeSurface(surface);
	SDL_SaveBMP(surface, "out.bmp");
//	writeln(to!string(SDL_GetError()));
}
