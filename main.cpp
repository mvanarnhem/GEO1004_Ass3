#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <cstdlib>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <limits>
#include <cmath>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef Kernel::Triangle_3                                  Triangle;
typedef CGAL::Bbox_3                                        Bbox_3;


struct Object {
    std::string name;
    std::vector<Triangle> shells;
    std::string material_type;
};

struct VoxelGrid {
    std::vector<unsigned int> voxels;
    unsigned int max_x, max_y, max_z;

    VoxelGrid(unsigned int x, unsigned int y, unsigned int z) {
        max_x = x;
        max_y = y;
        max_z = z;
        unsigned int total_voxels = x*y*z;
        voxels.reserve(total_voxels);
        for (unsigned int i = 0; i < total_voxels; ++i) voxels.push_back(0);
    }

    unsigned int &operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        assert(z >= 0 && z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }

    unsigned int operator()(const unsigned int &x, const unsigned int &y, const unsigned int &z) const {
        assert(x >= 0 && x < max_x);
        assert(y >= 0 && y < max_y);
        assert(z >= 0 && z < max_z);
        return voxels[x + y*max_x + z*max_x*max_y];
    }

};

void from_OBJ_to_Object(std::map<int, Object> &objects, std::vector<Point_3> &vertices, const std::string &input_file);
void get_extents(std::vector<Point_3> &vertices, std::map<std::string, int> &extent, double &voxel_size, Point_3 &origin);
void translate_RealWorldCoordinates_to_VoxelGridVoxel(Point_3 pt, double &voxel_size, Point_3 origin,
                                                    unsigned int &x, unsigned int &y, unsigned int &z);
Bbox_3 get_bbox_of_voxel(unsigned int &x, unsigned int &y, unsigned int &z,
                         double &voxel_size, Point_3 origin);
Point_3 translate_VoxelGridVoxel_to_RealWorldCoordinates(unsigned int &x, unsigned int &y,
                                                         unsigned int &z, double &voxel_size, Point_3 origin);
void generateOBJ(VoxelGrid voxelgrid, const std::string& objFilePath, Point_3 origin, double voxel_size, int value);
void intersection(std::map<int, Object> &objects, double &voxel_size, Point_3 &origin, VoxelGrid &my_building_grid);
void label_region(VoxelGrid &voxel_grid, unsigned int label, int start_voxel_index);

int main() {
    const std::string input_file = "../data/output_small_house.obj";
    std::map<int, Object> objects;
    std::vector<Point_3> vertices;
    from_OBJ_to_Object(objects, vertices, input_file);

    std::map<std::string, int> extent;
    double voxel_size = 0.5;
    Point_3 origin; // This is the point which is located left below in the voxelgrid. It contains the real coordinates
    get_extents(vertices, extent, voxel_size, origin);

    VoxelGrid my_building_grid(extent["x_range_VoxelGrid"], extent["y_range_VoxelGrid"], extent["z_range_VoxelGrid"]);

    intersection(objects, voxel_size, origin, my_building_grid);

    label_region(my_building_grid, 2, 0);

    generateOBJ(my_building_grid, "output_only_0.obj", origin, voxel_size, 0);
    generateOBJ(my_building_grid, "output_only_1.obj", origin, voxel_size, 1);
    generateOBJ(my_building_grid, "output_only_2.obj", origin, voxel_size, 2);
    return 0;
}

void label_region(VoxelGrid &voxel_grid, unsigned int label, int start_voxel_index) {
    // initialize the to_visit list of voxel indices we need to visit
    std::vector<int> to_visit;

    // define grid dimensions for further use
    int grid_size_x = voxel_grid.max_x;
    int grid_size_y = voxel_grid.max_y;
    int grid_size_z = voxel_grid.max_z;

    // extract coordinates of the starting index
    int start_x = start_voxel_index % grid_size_x;
    int start_y = (start_voxel_index / grid_size_x) % grid_size_y;
    int start_z = start_voxel_index / (grid_size_x * grid_size_y);

    std::cout << "\tGrowing region for label: " << label <<
              "\tstarting at: " << "(" << start_x << ", " << start_y << ", " << start_z << ")" << std::endl;

    to_visit.emplace_back(start_voxel_index); // place start voxel in to_visit
    voxel_grid(start_x, start_y, start_z) = label; // set it's value to the label

    // as long as there are voxels to visit, label their unlabelled neighbors and add them to to_visit
    while (!to_visit.empty()) {
        int current_voxel_id = to_visit.front();

        // find x, y, z of the current voxel
        int x = current_voxel_id % grid_size_x;
        int y = (current_voxel_id / grid_size_x) % grid_size_y;
        int z = current_voxel_id / (grid_size_x * grid_size_y);

        // loop over every neighboring voxel using 26-connectivity
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    // Skip the voxel itself
                    if (dx == 0 && dy == 0 && dz == 0)
                        continue;

                    // Calculate neighbor coordinates
                    int nx = x + dx;
                    int ny = y + dy;
                    int nz = z + dz;
                    int neighbor_id = nx + ny * grid_size_x + nz * grid_size_x * grid_size_y;

                    // if neighbor is witin grid bounds and is unlabelled
                    if (0 <= nx && nx < grid_size_x &&
                        0 <= ny && ny < grid_size_y &&
                        0 <= nz && nz < grid_size_z &&
                        voxel_grid(nx, ny, nz) == 0) {
                        // set the neighbor to the label and add it to the to_visit list
                        voxel_grid(nx, ny, nz) = label;
                        to_visit.push_back(neighbor_id);
                    }
                }
            }
        }
    std::cout << "to_visit length = " << to_visit.size() << std::endl;
    to_visit.erase(to_visit.begin());
    }
}




void from_OBJ_to_Object(std::map<int, Object> &objects, std::vector<Point_3> &vertices, const std::string &input_file){
    std::ifstream input_stream;
    input_stream.open(input_file);

    if (input_stream.is_open()) {
        std::string line;
        Object current_object;
        int index_object = -1;
        while (getline(input_stream, line)) {

            std::istringstream line_stream(line);
            std::string line_type;
            line_stream >> line_type;

            if (line_type == "g") {
                std::string name;
                if (index_object != -1){
                    objects[index_object] = current_object;
                    Object current_object;
                }
                line_stream >> current_object.name;
                index_object += 1;
            }
            else if (line_type == "v"){
                double x, y, z;
                line_stream >> x >> y >> z;
                Point_3 vert = Point_3(x, y, z);
                vertices.emplace_back(vert);
            }
            else if(line_type == "f"){
                std::string v1_str, v2_str, v3_str;
                int v1, v2, v3;
                line_stream >> v1_str >> v2_str >> v3_str;

                v1 = stoi(v1_str.substr(0, v1_str.find("//")));
                v2 = stoi(v2_str.substr(0, v2_str.find("//")));
                v3 = stoi(v3_str.substr(0, v3_str.find("//")));

                std::string remaining_part;
                line_stream >> remaining_part;
                if (!remaining_part.empty()){
                    std::cout << "The input face contains more than 3 vertices." << std::endl;
                }
                else{
                    current_object.shells.emplace_back(Triangle(vertices[v1-1], vertices[v2-1], vertices[v3-1]));
                }
            }
            else if(line_type == "usemtl"){
                line_stream >> current_object.material_type;
            }
        }
        objects[index_object] = current_object;
    }
}

void get_extents(std::vector<Point_3> &vertices, std::map<std::string, int> &extent, double &voxel_size, Point_3 &origin){
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::lowest();

    for(const auto& point : vertices) {
        if(point.x() < min_x) {min_x = point.x();}
        if(point.x() > max_x) {max_x = point.x();}
        if(point.y() < min_y) {min_y = point.y();}
        if(point.y() > max_y) {max_y = point.y();}
        if(point.z() < min_z) {min_z = point.z();}
        if(point.z() > max_z) {max_z = point.z();}
    }

    min_x = min_x - voxel_size;
    max_x = max_x + voxel_size;
    min_y = min_y - voxel_size;
    max_y = max_y + voxel_size;
    min_z = min_z - voxel_size;
    max_z = max_z + voxel_size;

    int x_range_VoxelGrid = static_cast<int>(std::ceil(std::abs( max_x - min_x) / voxel_size));
    int y_range_VoxelGrid = static_cast<int>(std::ceil(std::abs(max_y - min_y) / voxel_size));
    int z_range_VoxelGrid = static_cast<int>(std::ceil(std::abs(max_z - min_z) / voxel_size));

    extent["x_range_VoxelGrid"] = x_range_VoxelGrid;
    extent["y_range_VoxelGrid"] = y_range_VoxelGrid;
    extent["z_range_VoxelGrid"] = z_range_VoxelGrid;
    double x, y, z;
    x = min_x - 0.5 * (x_range_VoxelGrid * voxel_size - (max_x - min_x));
    y = min_y - 0.5 * (y_range_VoxelGrid * voxel_size - (max_y - min_y));
    z = min_z - 0.5 * (z_range_VoxelGrid * voxel_size - (max_z - min_z));
    origin = Point_3(x, y, z);
}

Point_3 translate_VoxelGridVoxel_to_RealWorldCoordinates(unsigned int &x, unsigned int &y, unsigned int &z,
                                                       double &voxel_size, Point_3 origin){
    double centroid_x_coordinate = origin.x() + (x + 0.5) * voxel_size;
    double centroid_y_coordinate = origin.y() + (y + 0.5) * voxel_size;
    double centroid_z_coordinate = origin.z() + (z + 0.5) * voxel_size;
    return Point_3(centroid_x_coordinate, centroid_y_coordinate, centroid_z_coordinate);
}

Bbox_3 get_bbox_of_voxel(unsigned int &x, unsigned int &y, unsigned int &z,
                                                       double &voxel_size, Point_3 origin){
    double x_min = origin.x() + x * voxel_size;
    double y_min = origin.y() + y * voxel_size;
    double z_min = origin.z() + z * voxel_size;
    double x_max = origin.x() + (x+1) * voxel_size;
    double y_max = origin.y() + (y+1) * voxel_size;
    double z_max = origin.z() + (z+1) * voxel_size;
    return Bbox_3(x_min, y_min, z_min, x_max, y_max, z_max);
}

void translate_RealWorldCoordinates_to_VoxelGridVoxel(Point_3 pt, double &voxel_size, Point_3 origin,
                                                    unsigned int &x, unsigned int &y, unsigned int &z){
    x = static_cast<unsigned int>(std::floor((pt.x() - origin.x()) / voxel_size));
    y = static_cast<unsigned int>(std::floor((pt.y() - origin.y()) / voxel_size));
    z = static_cast<unsigned int>(std::floor((pt.z() - origin.z()) / voxel_size));
}

void generateOBJ(VoxelGrid voxel_grid, const std::string& objFilePath, Point_3 origin, double voxel_size, int value) {
    std::ofstream objFile(objFilePath);

    if (!objFile.is_open()) {
        std::cerr << "Error: Unable to open file: " << objFilePath << std::endl;
        return;
    }

    // Define material library
    objFile << "mtllib colors.mtl\n\n";

    for (unsigned int x = 0; x < voxel_grid.max_x; ++x) {
        for (unsigned int y = 0; y < voxel_grid.max_y; ++y) {
            for (unsigned int z = 0; z < voxel_grid.max_z; ++z) {
                // Check if voxel value is 1
                if (voxel_grid(x, y, z) == value) {
                    Point_3 pt = translate_VoxelGridVoxel_to_RealWorldCoordinates(x, y, z, voxel_size, origin);
                    // Write vertex coordinates
                    objFile << "v " << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
                }
            }
        }
    }

    objFile.close();
}

void intersection(std::map<int, Object> &objects, double &voxel_size, Point_3 &origin, VoxelGrid &my_building_grid){
    for (const auto& entry : objects) {
        const Object &obj = entry.second;

        for (const auto &triangle: obj.shells) {
            Bbox_3 triangle_bbox = triangle.bbox();

            Point_3 min_corner = Point_3(triangle_bbox.xmin(), triangle_bbox.ymin(), triangle_bbox.zmin());
            Point_3 max_corner = Point_3(triangle_bbox.xmax(), triangle_bbox.ymax(), triangle_bbox.zmax());

            unsigned int x_min, y_min, z_min;
            unsigned int x_max, y_max, z_max;

            translate_RealWorldCoordinates_to_VoxelGridVoxel(min_corner, voxel_size, origin, x_min, y_min, z_min);
            translate_RealWorldCoordinates_to_VoxelGridVoxel(max_corner, voxel_size, origin, x_max, y_max, z_max);

            // loop over voxels within bbox of triangle
            for (unsigned int i = x_min; i <= x_max; i++) {
                for (unsigned int j = y_min; j <= y_max; j++) {
                    for (unsigned int k = z_min; k <= z_max; k++) {
                        //get bbox of the voxel
                        Bbox_3 voxel_bbox = get_bbox_of_voxel(i, j, k, voxel_size, origin);
                        // Check intersection
                        bool intersects = CGAL::do_intersect(voxel_bbox, triangle);
                        if (intersects) { // there is a triangle in this voxel, set voxel value to 1
                            my_building_grid(i,j,k) = 1;
                        }
                    }
                }
            }
        }
    }
}