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

#define NB_THREAD 8
#define NB_CLUSTER_PER_THREAD 64

#define NB_TOTAL_CLUSTER (NB_THREAD*NB_CLUSTER_PER_THREAD)
#define NB_PARTICLES 10000


#define MAX_PARTICULES_PER_CLUSTER (2*(NB_PARTICLES/NB_TOTAL_CLUSTER))
#define MIN_PARTICULES_PER_CLUSTER ((NB_PARTICLES/NB_TOTAL_CLUSTER)/2)

#define MIN_DISTANCE_FUSE 20

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

class Cluster
{
public:
    int complexity;
    int particules[MAX_PARTICULES_PER_CLUSTER];
    float distance_from_center[MAX_PARTICULES_PER_CLUSTER];
    int nb_particules;

    float total_mass;
    Vector3 center;
    float radius;

    void init_total_mass(Particule *particules_array)
    {
        total_mass = 0;
        for(int i = 0;i != nb_particules;++i)
            total_mass += particules_array[particules[i]].mass;
    }

    void calculate_properties(Particule *particules_array)
    {
        calculate_center(particules_array);
        radius = 0;
        for(int i = 0;i != nb_particules;++i)
        {
            Particule &current_particule = particules_array[particules[i]];
            float delta_x = std::abs(center.x-current_particule.position.x);
            float delta_y = std::abs(center.y-current_particule.position.y);
            float delta_z = std::abs(center.z-current_particule.position.z);
            distance_from_center[i] = std::max(delta_x,std::max(delta_y,delta_z));
            if(distance_from_center[i] > radius) radius = distance_from_center[i];
        }
    }


    void calculate_center(Particule *particules_array)
    {
        center.set(0,0,0);
        for(int i = 0;i != nb_particules;++i)
        {
            Particule &current_particule = particules_array[particules[i]];
            center += current_particule.position*current_particule.mass;
        }
        center /= total_mass;
    }

    void steal(Cluster &other_cluster,int index,Particule *particules_array)
    {
        int absolute_id = other_cluster.particules[index];
        Particule &stealed_particule = particules_array[absolute_id];
        other_cluster.total_mass -= stealed_particule.mass;
        total_mass += stealed_particule.mass;
        particules[nb_particules] = absolute_id;

        std::memmove(&other_cluster.particules[index],
                    &other_cluster.particules[index+1],
                    sizeof(int)*(other_cluster.nb_particules-index-1));


        nb_particules++;
        other_cluster.nb_particules--;
        calculate_center(particules_array);
        other_cluster.calculate_center(particules_array);

    }


    Cluster() {
        nb_particules = 0;
    }
};

#endif // PARTICLES_HPP_
