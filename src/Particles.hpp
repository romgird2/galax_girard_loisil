#ifndef PARTICLES_HPP_
#define PARTICLES_HPP_

#include <vector>
#include <cmath>
#include <iostream>

struct Particles
{
	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> z;

	Particles(const int n_particles);
};

#define NB_THREAD 8
#define NB_MAX_CLUSTER_PER_THREAD 16

#define NB_MAX_CLUSTER NB_THREAD*NB_MAX_CLUSTER_PER_THREAD
#define NB_PARTICLES 10000


#define MAX_PARTICULES_PER_CLUSTER 2*(NB_PARTICLES/NB_MAX_CLUSTER)
#define MIN_PARTICULES_PER_CLUSTER 3

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

    Vector3 operator-(const Vector3& other) const {
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    Vector3 operator*(float scalar) const {
        return Vector3(x * scalar, y * scalar, z * scalar);
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
    int nb_particules;
    float total_mass;

    Vector3 center;


    Cluster() {
        nb_particules = 0;
    }
};

#endif // PARTICLES_HPP_
