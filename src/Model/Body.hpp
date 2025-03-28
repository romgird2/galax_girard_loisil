#ifndef BODY_HPP
#define BODY_HPP

#include <immintrin.h>
constexpr double G = 6.67e-11;

/**
 * @brief The Vector3 class
 * Represents a vector in a 3D space
 */
class Vector3
{
public:
    union {
        struct { float x, y, z, w; }; // w for padding/alignment
        __m128 simd; // SIMD register type
    };

    // Calculate ||this-v||^2
    float dist_sq(Vector3 const& v) const;

    // Calculate (this+v)/2
    Vector3 med(Vector3 const& v);

    Vector3 operator-(Vector3 const& v) const;

    Vector3 operator+(Vector3 const& v) const;

    void operator+=(Vector3 const& v);

    Vector3 operator*(float f) const;

    void operator/=(float f);

    bool operator==(Vector3 const& v) const;

    // Calculate ||this||^2
    float norm_sqr() const;
};

/**
 * @brief The Body class
 * Represent a body with its mas
 */
struct Body
{
    Vector3 pos;
    Vector3 spd;
    Vector3 acceleration;

    double mass;

    double dist_sq(Body const& b) const;

    // Update the acceleration given the force between this and b
    void update_force(Body const& b);

    // Update the position by using the current acceleration
    void update_pos();

    // Transform this to the center of mass of b and this
    void operator+=(Body const& b);

    // Check if two bodies are at the same position
    bool operator==(Body const& b) const;
};

#endif // BODY_HPP
