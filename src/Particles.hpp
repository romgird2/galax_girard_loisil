#ifndef PARTICLES_HPP_
#define PARTICLES_HPP_

#include <vector>
#include <cmath>
#include <iostream>
#include <cstring>

struct Particles
{
	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> z;

	Particles(const int n_particles);
};

#define NB_PARTICLES 10000
#define NB_THREADS 8


#define PARTICULES_BEFORE_SEPARATION 2

class Vector3
{
public:
    float x,y,z;
    Vector3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}



    Vector3 operator+(const Vector3& other) const {
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    void operator+=(const Vector3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }

    void operator-=(const Vector3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }

    void operator*=(float scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
    }

    void operator/=(float scalar) {
        x /= scalar;
        y /= scalar;
        z /= scalar;
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
    }

    Vector3 operator/(float scalar) const {
        return Vector3(x / scalar, y / scalar, z / scalar);
    }

    void normalize() {
        float length = norm();
        if (length > 0) {
            x /= length;
            y /= length;
            z /= length;
        }
    }

    float normSquared() const {
        return x * x + y * y + z * z;
    }

    float norm() const {
        return std::sqrt(normSquared());
    }

    void set(float newX, float newY, float newZ) {
        x = newX;
        y = newY;
        z = newZ;
    }

    void set(const Vector3& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }

    void display() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")";
    }

};

class Particule {

public:
    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration;
    float mass;

    Particule()
        : position(), velocity(), acceleration() {}

    void resetAcceleration()
    {
        acceleration.set(0,0,0);
    }
};
#include <stdlib.h>
class OctTree
{
public:
    static Particule *particules;
    Vector3 center;
    float radius;
    Vector3 center_mass;
    float mass;
    OctTree *children;
    bool leaf;
    int particules_index[PARTICULES_BEFORE_SEPARATION-1];
    int nb_particules;


    OctTree() {

    }
    void init(Vector3 ncenter,float nradius)
    {
        center = ncenter;
        radius = nradius;
        leaf = true;
        nb_particules = 0;
    }

    void clear()
    {
        if(!leaf)
        {
            for(int i = 0;i != 8;++i)
                children[i].clear();
            delete children;
        }
    }

    void pre_compute()
    {
        mass = 0;
        center_mass.set(0,0,0);
        if(leaf)
        {
            for(int i = 0;i != nb_particules;++i)
            {
                float mass_particule = particules[particules_index[i]].mass;
                mass += mass_particule;
                center_mass += particules[particules_index[i]].position*mass_particule;
            }
            if(mass != 0)
                center_mass /= mass;
        }
        else
        {
            for(int i = 0;i != 8;++i)
            {
                children[i].pre_compute();
                float mass_children = children[i].mass;
                mass += mass_children;
                center_mass += children[i].center_mass*mass_children;
            }
            if(mass != 0)
                center_mass /= mass;
        }
    }

    #define THETA 0
    #define THETA_SQR (THETA*THETA)

    void compute_with(OctTree &target)
    {
        if((target.radius+radius) < (target.center-radius).norm()*THETA)
        {
            Vector3 diff = target.center_mass-center_mass;
            float dij = diff.normSquared();

            if (dij < 1.0)
            {
                dij = 10.0;
            }
            else
            {
                dij = std::sqrt(dij);
                dij = 10.0 / (dij * dij * dij);
            }

            Vector3 delta = diff * dij * target.mass;
            for(int i = 0;i != nb_particules;++i)
                particules[particules_index[i]].acceleration += delta;

        }
        else
        {
            if(target.leaf)
            {
                for(int i = 0;i != nb_particules;++i)
                {
                    for(int j = 0;j != target.nb_particules;++j)
                    {
                        int particule_index_i = particules_index[i];
                        int particule_index_j = target.particules_index[j];
                        if(i != j)
                        {
                            Particule &particulei = particules[particule_index_i];
                            Particule &particulej = particules[particule_index_j];

                            Vector3 diff = particulej.position-particulei.position;
                            float dij = diff.normSquared();

                            if (dij < 1.0)
                            {
                                dij = 10.0;
                            }
                            else
                            {
                                dij = std::sqrt(dij);
                                dij = 10.0 / (dij * dij * dij);
                            }

                            particulei.acceleration += diff * dij * particulej.mass;
                        }

                    }
                }
            }
            else
            {
                for(int i = 0;i != 8;++i)
                {
                    OctTree &new_target = target.children[i];
                    if(new_target.nb_particules != 0)
                        compute_with(new_target);
                }
            }
        }
    }

    void compute_acceleration(OctTree &root)
    {
        if(nb_particules == 0) return;
        if(leaf)
        {
            compute_with(root);
        }
        else
        {
            for(int i = 0;i != 8;++i)
            {
                children[i].compute_acceleration(root);
            }
        }
    }

    void print()
    {
        print(0);
    }

    void print(int profondeur)
    {
        if(leaf)
        {
            for(int i = 0;i != profondeur;++i)
                std::cout << "\t";
            std::cout << "[";
            for(int i = 0;i != nb_particules;++i)
            {
                std::cout << particules_index[i];
                if(i != nb_particules-1)
                    std::cout << ";";
            }
            std::cout << "]";
            std::cout << std::endl;
        }
        else
        {
            for(int i = 0;i != profondeur;++i)
                std::cout << "\t";
            std::cout << "[]";
            for(int i = 0;i != 8;++i)
            {
                children[i].print(profondeur+1);
            }
        }
    }

    void unleaf()
    {
        float half_radius = radius/2;
        float bottom_x = center.x-half_radius;
        float bottom_y = center.y-half_radius;
        float bottom_z = center.z-half_radius;
        float top_x = center.x+half_radius;
        float top_y = center.y+half_radius;
        float top_z = center.z+half_radius;


        children = new OctTree[8]; // make reinstantiate

        children[0].init(Vector3(bottom_x,bottom_y,bottom_z),half_radius);
        children[1].init(Vector3(top_x,bottom_y,bottom_z),half_radius);
        children[2].init(Vector3(bottom_x,top_y,bottom_z),half_radius);
        children[3].init(Vector3(top_x,top_y,bottom_z),half_radius);

        children[4].init(Vector3(bottom_x,bottom_y,top_z),half_radius);
        children[5].init(Vector3(top_x,bottom_y,top_z),half_radius);
        children[6].init(Vector3(bottom_x,top_y,top_z),half_radius);
        children[7].init(Vector3(top_x,top_y,top_z),half_radius);
        leaf = false;

    }


    void put_into_children(int &particule_index,Vector3 &position)
    {
        int children_index = (position.x > center.x) | ((position.y > center.y)<<1) | ((position.z > center.z)<<2);
        OctTree &children_target = children[children_index];
        children_target.particules_index[children_target.nb_particules] = particule_index;
        children_target.nb_particules++;
    }


    void insert(Vector3& position,int& particule)
    {
        if(leaf)
        {
            //std::cout << nb_particules << std::endl;
            if(nb_particules == PARTICULES_BEFORE_SEPARATION-1)
            {
                unleaf();
                //std::cout << "unleafing" << std::endl;
                for(int i = 0;i != nb_particules;++i)
                {
                    int particule_index = particules_index[i];
                    Vector3 &position_target = particules[particules_index[i]].position;
                    put_into_children(particule_index,position_target);
                }
                int children_index = (position.x > center.x) | ((position.y > center.y)<<1) | ((position.z > center.z)<<2);
                children[children_index].insert(position,particule);
            }
            else
            {
                //std::cout << "unleafing" << std::endl;
                particules_index[nb_particules] = particule;
                nb_particules++;
            }
        }
        else
        {
            int children_index = (position.x > center.x) | ((position.y > center.y)<<1) | ((position.z > center.z)<<2);
            children[children_index].insert(position,particule);
        }
    }
};

#endif // PARTICLES_HPP_
