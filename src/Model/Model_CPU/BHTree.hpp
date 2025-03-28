#ifndef BHTREE_H
#define BHTREE_H

#include <array>
#include <memory>
#include "../Body.hpp"

constexpr double theta = 0.5;

constexpr double theta_sqr = theta*theta;

// Direction. First bit x, second bit y, third bit z
enum Direction : int {
    UNW,
    UNE,
    USW,
    USE,
    LNW,
    LNE,
    LSW,
    LSE
};

/**
 * @brief The Region class
 * An octant of the tree
 */
struct Region
{
    Vector3 center;
    double width;
    double width_sqr;
    // Get the sub octant in the direction dir relative to the center of the region
    Region get_sub(Direction dir);
    void to_string();
};

/**
 * @brief The BHTree class
 * A octree for Barnes-Hut
 */
class BHTree
{
public:
    BHTree(Region region);
    // Build a tree from its children
    BHTree(Region region, std::array<std::unique_ptr<BHTree>,8>& children);
    void insert(Body b);
    bool is_leaf();
    void update_force(Body& b);
    int height();
    // Print the three in a file, for debug
    void print_tree();
    // Explore the tree with output in a file, for debug
    void explore(std::ofstream& fout, int depth);
    Body mass_center;
    void to_string();
private:
    Region region;
    std::array<std::unique_ptr<BHTree>,8> children;
    bool leaf;
    void add(Body b);
};

#endif // BHTREE_H
