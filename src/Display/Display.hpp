#ifndef DISPLAY_HPP_
#define DISPLAY_HPP_

#include "../Particles.hpp"
#include "../Model/Model_CPU/Model_CPU_fast/Model_CPU_fast.hpp"

class Display
{
protected:
    Model_CPU_fast *model;

public:
    Display(Model_CPU_fast *model);
	~Display();

	virtual void update(bool& done) = 0;
};

#endif // DISPLAY_HPP_
