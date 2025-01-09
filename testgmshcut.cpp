#include "/usr/local/include/gmsh.h"
#include <iostream>

int main() {
    gmsh::initialize();
    gmsh::model::add("cut_example");

    // Create a rectangle and a disk
    int rectTag = gmsh::model::occ::addRectangle(0, 0, 0, 10, 5);
    int diskTag = gmsh::model::occ::addDisk(5, 2.5, 0, 2);

    // Synchronize to update the model
    gmsh::model::occ::synchronize();

    // Perform the cut operation without removing objects or tools
    std::vector<std::pair<int, int>> outSurfaces;
    gmsh::model::occ::cut({{2, rectTag}}, {{2, diskTag}}, outSurfaces, {}, false, false);

    // Synchronize again to register the cut operation
    gmsh::model::occ::synchronize();

    // Print resulting surfaces
    std::cout << "Resulting surfaces from the cut:\n";
    for (const auto &surface : outSurfaces) {
        std::cout << "  Dimension: " << surface.first << ", Tag: " << surface.second << "\n";
    }

    // Get all entities in the model
    auto entities = gmsh::model::getEntities();
    std::cout << "\nAll entities in the model:\n";
    for (const auto &entity : entities) {
        std::cout << "  Dimension: " << entity.first << ", Tag: " << entity.second << "\n";
    }

    gmsh::finalize();
    return 0;
}