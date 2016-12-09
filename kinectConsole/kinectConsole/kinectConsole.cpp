#include "stdafx.h"
#include <iostream>
#include <chrono>

#include <SDL.h>

#include "app.h"

using namespace std::chrono;
typedef steady_clock Clock;

void DrawPixelbuffer(SDL_Texture* texture, SDL_Renderer* renderer,
	const uint32* pixelBuffer)
{
	//upload the pixel buffer to a texture
	void* pixels;
	int pitch;
	SDL_LockTexture(texture, nullptr, &pixels, &pitch);
	if (pitch == SCRWIDTH * 4)
		memcpy(pixels, pixelBuffer, SCRWIDTH * SCRHEIGHT * 4);
	else
	{
		const uint32* src = pixelBuffer;
		char* dst = (char*)pixels;
		for (int y = 0; y < SCRHEIGHT; ++y)
		{
			memcpy(dst, src, SCRWIDTH * 4);
			src += SCRWIDTH;
			dst += pitch;
		}
	}
	SDL_UnlockTexture(texture);

	//copy the texture to the frame buffer
	SDL_RenderCopy(renderer, texture, nullptr, nullptr);

	//present the frame buffer on the screen
	SDL_RenderPresent(renderer);
}

#undef main
int main(int, char**){
	//initialize SDL
	SDL_Init(SDL_INIT_VIDEO);

	//create a window
	SDL_Window* window = SDL_CreateWindow(
		"title", 100, 100, SCRWIDTH, SCRHEIGHT, 0);
	if (window == nullptr)
		return 1;

	//create a renderer
	SDL_Renderer* renderer = SDL_CreateRenderer(
		window, -1, SDL_RENDERER_ACCELERATED);
	if (renderer == nullptr)
		return 2;

	//create a texture
	SDL_Texture* texture = SDL_CreateTexture(
		renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING,
		SCRWIDTH, SCRHEIGHT);
	if (texture == nullptr)
		return 3;

	//allocate a pixel buffer
	uint32* pixelBuffer = new uint32[SCRWIDTH * SCRHEIGHT];
	if (pixelBuffer == nullptr)
		return 4;

	//clear the pixel buffer
	memset(pixelBuffer, 0, SCRWIDTH * SCRHEIGHT * 4);

	//draw pixel buffer to the screen
	DrawPixelbuffer(texture, renderer, pixelBuffer);

	App app;
	app.SetPixelBuffer(pixelBuffer);
	app.Init();

	auto lastTime = Clock::now();

	bool running = true;
	while (running)
	{
		//get events
		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			switch (event.type)
			{
				//pressing the cross or pressing escape will quit the application
			case SDL_QUIT:
				running = false;
				break;

			case SDL_KEYDOWN:
				if (event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)
					running = false;
				break;

			default: //ignore other events for this demo
				break;
			}
		}

		//calculate delta time
		const auto now = Clock::now();
		const auto duration = duration_cast<microseconds>(now - lastTime);
		const float deltaTime = duration.count() / 1000000.0f;
		lastTime = now;

		//update the application
		app.Tick(deltaTime);

		//draw pixel buffer to the screen
		DrawPixelbuffer(texture, renderer, pixelBuffer);
	}


	//clean up
	app.Shutdown();
	SDL_DestroyTexture(texture);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	return 0;
}



