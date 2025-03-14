#include <string>
#include <iostream>
#include <vector>
#include <iomanip>

#include <CLI/CLI.hpp>

#include "Display/Display_NO/Display_NO.hpp"
#ifdef GALAX_DISPLAY_SDL2
#include "Display/Display_SDL2/Display_SDL2.hpp"
#endif
#include "Timing/Timing.hpp"
#include "Initstate.hpp"
#include "Model/Model_CPU/Model_CPU_naive/Model_CPU_naive.hpp"
#include "Model/Model_CPU/Model_CPU_fast/Model_CPU_fast.hpp"
#include "Model/Model_GPU/Model_GPU.hpp"

int main(int argc, char ** argv)
{
	// class for CLI (Command Line Instructions) management
	CLI::App app{"Galax"};

	// maximum number of particles to be simulated
	const int max_n_particles = 81920;

	// according to compile option (in cmake), use a graphical display or don't
#ifdef GALAX_DISPLAY_SDL2
	std::string  display_type = "SDL2";
#else
	std::string  display_type = "NO";
#endif

	// core version used by default : CPU
    std::string  core         = "CPU";

	// number of particles used by default : 2000
    unsigned int n_particles  = 10000;

        // decide wether to check particle position against the reference or not
        bool validatePositions = false;

	// define CLI arguments
	app.add_option("-c,--core"       , core       , "computing version")
	    ->check(CLI::IsMember({"CPU", "GPU", "CPU_FAST"}));
	app.add_option("-n,--n-particles", n_particles , "number of displayed particles")
	    ->check(CLI::Range(0,max_n_particles));
	app.add_option("--display"       , display_type, "disable graphical display")
	    ->check(CLI::IsMember({"SDL2", "NO"}));
	app.add_flag("--validate", validatePositions, "compute error in positions against the non-accelerated reference code");

	// parse arguments
	CLI11_PARSE(app, argc, argv);

        // No need to validate if we are using the ref code as main simulation
        //validatePositions = !(core == "CPU");

	// class used to measure timing and fps
	Timing timing;

	// load particles initial position into initstate
	Initstate initstate(n_particles);

	// particles positions
	Particles particles(n_particles);
	Particles particlesRef(n_particles);

	// init display


	// init models
    std::unique_ptr<Model_CPU_fast> model;
    std::unique_ptr<Model> referenceModel;

    if(validatePositions)
		referenceModel = std::make_unique<Model_CPU_naive>(initstate, particlesRef);

    model = std::make_unique<Model_CPU_fast>(initstate, particles);

    std::unique_ptr<Display> display;
    if (display_type == "NO")
        display = std::unique_ptr<Display>(new Display_NO(model->particules,model->clusters));
#ifdef GALAX_DISPLAY_SDL2
    else if (display_type == "SDL2")
        display = std::unique_ptr<Display>(new Display_SDL2(model->particules,model->clusters));
#endif
    else // TODO : add exception
        exit(EXIT_FAILURE);

	bool done = false;

	std::cout << std::setw(3);

    double sum_fps = 0;
    int nb_fps = 0;
	while (!done)
	{
		// display particles
		display->update(done);

		// We only want to time the computation of the model
		// not its display
		timing.sample_before();

		// update particles positions
        //model  ->step();

		timing.sample_after();
		float fps = timing.get_current_average_FPS();

		std::cout << "State updates per second: " << fps;
        sum_fps += fps;
        nb_fps++;
        std::cout << "  average:" << sum_fps/nb_fps;

		if(validatePositions)
		{
			referenceModel->step();
                        float average_error, error_min, error_max;
                        std::tie(error_min, error_max, average_error) = model->compareParticlesState(*referenceModel, /*returnRelativeDistances*/ true);
			std::cout << " ;               average distance vs reference: " << average_error
                                  << "; min error : " << error_min << "; max error : " << error_max;
		}
		std::cout << "\r" << std::flush;
	}

	return EXIT_SUCCESS;
}

