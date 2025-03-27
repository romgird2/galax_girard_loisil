#include <cmath>
#include <cstdio>
#ifdef GALAX_DISPLAY_SDL2

#include "Display_SDL2.hpp"

Display_SDL2
::Display_SDL2(Model_CPU_fast *model)
: Display(model)
{
	SDL_DisplayMode current;

	if (SDL_Init (SDL_INIT_EVERYTHING) < 0)
	{
		printf("error: unable to init sdl\n");
        // TODO throw exception
        exit(EXIT_FAILURE);
	}

	if (SDL_GetDesktopDisplayMode(0, &current))
	{
		printf("error: unable to get current display mode\n");
        // TODO throw exception
        exit(EXIT_FAILURE);
	}

	window = SDL_CreateWindow("SDL", 	SDL_WINDOWPOS_CENTERED,
										SDL_WINDOWPOS_CENTERED,
										width, height,
										SDL_WINDOW_OPENGL);

	glWindow = SDL_GL_CreateContext(window);

	GLenum status = glewInit();

	if (status != GLEW_OK)
	{
		printf("error: unable to init glew\n");
		// TODO throw exception
		exit(EXIT_FAILURE);
	}

	SDL_GL_SetSwapInterval(1);
}

Display_SDL2::~Display_SDL2()
{
	SDL_GL_DeleteContext(glWindow);
	SDL_DestroyWindow(window);
	SDL_Quit();
}
struct Color {
    float r, g, b;
};

Color generate_color(int id, int max_id) {
    float hue = (static_cast<float>(id) * 360.0f) / static_cast<float>(max_id);
    float saturation = 0.8f;
    float value = 0.95f;

    float c = value * saturation;
    float h_prime = hue / 60.0f;
    float x = c * (1.0f - std::abs(std::fmod(h_prime, 2.0f) - 1.0f));
    float m = value - c;

    int i = static_cast<int>(h_prime) % 6;
    float r, g, b;

    switch (i) {
    case 0: r = c; g = x; b = 0; break;
    case 1: r = x; g = c; b = 0; break;
    case 2: r = 0; g = c; b = x; break;
    case 3: r = 0; g = x; b = c; break;
    case 4: r = x; g = 0; b = c; break;
    default: r = c; g = 0; b = x; break;
    }

    r += m;
    g += m;
    b += m;

    return {r, g, b};
}

void display_particules(OctTree &tree,Color& color)
{
    if(tree.leaf)
    {
        for(int i = 0;i != tree.nb_particules;++i)
        {
            Particule &particule = OctTree::particules[tree.particules_index[i]];
            glBegin   (GL_POINTS);
            glColor3f (color.r,color.g,color.b);
            glVertex3f(particule.position.x, particule.position.y, particule.position.z);
            glEnd();

        }
    }
    else
    {
        for(int i = 0;i != 8;++i)
        {
            display_particules(tree.children[i],color);
        }

    }
}

int hash32shift(int key)
{
    key = ~key + (key << 15); // key = (key << 15) - key - 1;
    key = key ^ (key >> 12);
    key = key + (key << 2);
    key = key ^ (key >> 4);
    key = key * 2057; // key = (key + (key << 3)) + (key << 11);
    key = key ^ (key >> 16);
    return key;
}

void display_octree(OctTree &tree,int index,int profondeur,int delta,int max_index)
{
    if(profondeur == 0)
    {
        Color color = generate_color(hash32shift(index),max_index);
        display_particules(tree,color);
    }
    else
    {
        if(tree.leaf)
        {
            //Color color = generate_color(index,max_index);
            //display_particules(tree,color);
        }
        else
        {
            for(int i = 0;i != 8;++i)
            {
                display_octree(tree.children[i],index+i*delta,profondeur-1,delta>>3,max_index);
            }
        }
    }
}

int profondeur = 0;

void Display_SDL2
::update(bool& done)
{
	int i;

	while (SDL_PollEvent(&event))
	{

		unsigned int e = event.type;

		if (e == SDL_MOUSEMOTION)
		{
			mouseMoveX = event.motion.x;
			mouseMoveY = height - event.motion.y - 1;
		}
		else if (e == SDL_KEYDOWN)
		{
			if (event.key.keysym.sym == SDLK_F1)
				g_showGrid = !g_showGrid;
			else if (event.key.keysym.sym == SDLK_F2)
				g_showAxes = !g_showAxes;
			else if (event.key.keysym.sym == SDLK_ESCAPE)
				done = true;
            else if(event.key.keysym.sym == SDLK_UP)
                profondeur++;
            else if(event.key.keysym.sym == SDLK_DOWN)
                profondeur--;
        }

		if (e == SDL_QUIT)
		{
			printf("quit\n");
			done = true;
		}
	}

	mouseDeltaX = mouseMoveX - mouseOriginX;
	mouseDeltaY = mouseMoveY - mouseOriginY;

	if (SDL_GetMouseState(0, 0) & SDL_BUTTON_LMASK)
	{
		oldCamRot[ 0 ] += -mouseDeltaY / 5.0f;
		oldCamRot[ 1 ] += mouseDeltaX / 5.0f;
	}
	else if (SDL_GetMouseState(0, 0) & SDL_BUTTON_RMASK)
	{
		oldCamPos[ 2 ] += (mouseDeltaY / 100.0f) * 0.5 * fabs(oldCamPos[ 2 ]);
		oldCamPos[ 2 ]  = oldCamPos[ 2 ] > -5.0f ? -5.0f : oldCamPos[ 2 ];
	}

	mouseOriginX = mouseMoveX;
	mouseOriginY = mouseMoveY;

	glViewport     (0, 0, width, height);
	glClearColor   (0.2f, 0.2f, 0.2f, 1.0f);
	glClear        (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable       (GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc    (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable      (GL_TEXTURE_2D);
	glEnable       (GL_DEPTH_TEST);
	glMatrixMode   (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective (50.0f, (float)width / (float)height, 0.1f, 100000.0f);
	glMatrixMode   (GL_MODELVIEW);
	glLoadIdentity ();

	for (i = 0; i < 3; ++i)
	{
		newCamPos[i] += (oldCamPos[i] - newCamPos[i]) * g_inertia;
		newCamRot[i] += (oldCamRot[i] - newCamRot[i]) * g_inertia;
	}

	glTranslatef(newCamPos[0], newCamPos[1], newCamPos[2]);
	glRotatef   (newCamRot[0], 1.0f, 0.0f, 0.0f);
	glRotatef   (newCamRot[1], 0.0f, 1.0f, 0.0f);

	if (g_showGrid)
		DrawGridXZ(-100.0f, 0.0f, -100.0f, 20, 20, 10.0);

	if (g_showAxes)
		ShowAxes();
/*
    for(int i = 0;i != NB_TOTAL_CLUSTER;++i)
    {
        Cluster &cluster = clusters[i];
        float t = NB_TOTAL_CLUSTER > 1 ? (float) i / (NB_TOTAL_CLUSTER-1) : 0.5;
        Color color = generate_color(i,NB_TOTAL_CLUSTER);
        for(int j = 0;j != cluster.nb_particules;++j)
        {
            Particule &particule = particules[cluster.particules[j]];
            glBegin   (GL_POINTS);
            glColor3f (color.r, color.g, color.b);
            glVertex3f(particule.position.x, particule.position.y, particule.position.z);
            glEnd();

        }
    }*/




    /*for(int i = 0;i != NB_PARTICLES;++i)
    {
        Particule &particule = model->particules[i];
        glBegin   (GL_POINTS);
        glColor3f (1, 1, 1);
        glVertex3f(particule.position.x, particule.position.y, particule.position.z);
        glEnd();
    }*/


    display_octree(model->root,0,profondeur,(1<<(profondeur*3))>>3,1<<(profondeur*3));

	glMatrixMode  (GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D    (0, width, 0, height);
	glMatrixMode  (GL_MODELVIEW);
	glLoadIdentity();

	SDL_GL_SwapWindow(window);
	SDL_UpdateWindowSurface(window);
}




void Display_SDL2
::DrawPoint(float x, float y, float z) const
{
		glBegin(GL_POINTS);
		glColor3f(1.0f, 1.0f, 1.0f);
		glVertex3f(x, y, z);
		glEnd();
}

void Display_SDL2
::DrawGridXZ(float ox, float oy, float oz, int w, int h, float sz) const
{

	glLineWidth(1.0f);
	glBegin    (GL_LINES);
	glColor3f  (0.48f, 0.48f, 0.48f);

	for (auto i = 0; i <= h; ++i)
	{
		glVertex3f(ox, oy, oz + i * sz);
		glVertex3f(ox + w * sz, oy, oz + i * sz);
	}

	for (auto i = 0; i <= h; ++i)
	{
		glVertex3f(ox + i * sz, oy, oz);
		glVertex3f(ox + i * sz, oy, oz + h * sz);
	}

	glEnd();
}

void Display_SDL2
::ShowAxes() const
{
	glLineWidth(2.0f);
	glBegin(GL_LINES);

	glColor3f (1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(2.0f, 0.0f, 0.0f);

	glColor3f (0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 2.0f, 0.0f);

	glColor3f (0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 2.0f);

	glEnd();
}

#endif
