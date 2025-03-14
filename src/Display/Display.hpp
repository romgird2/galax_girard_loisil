#ifndef DISPLAY_HPP_
#define DISPLAY_HPP_

#include "../Particles.hpp"

class Display
{
protected:
    Particule *particules;
    Cluster *clusters;

public:
    Display(Particule *particules,Cluster *clusters);
	~Display();

	virtual void update(bool& done) = 0;
};

#endif // DISPLAY_HPP_
