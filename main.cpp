#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Plane_3.h>
#include <cstdlib>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <typeinfo>
#include <limits>
#include <cmath>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3                      Plane;
typedef Kernel::Point_3                      Point_3;
typedef Kernel::Triangle_3                   Triangle;

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
void get_extents(std::vector<Point_3> &vertices, std::map<std::string, int> &extent, double &voxel_size, Point_3 &model_coordinates);
Point_3 translate_from_XYZVoxelgrid_to_CentroidInModel(unsigned int &x, unsigned int &y,
                                                       unsigned int &z, double &voxel_size, Point_3 model_coordinates);
void translate_from_CentroidInModel_to_XYZVoxelgrid(Point_3 pt, double &voxel_size, Point_3 model_coordinates, unsigned int &x, unsigned int &y,
                                                                unsigned int &z);

int main() {
    const std::string input_file = "../data/output_small_house.obj";
    std::map<int, Object> objects;
    std::vector<Point_3> vertices;
    from_OBJ_to_Object(objects, vertices, input_file);

    std::map<std::string, int> extent;
    double voxel_size = 0.5;
    Point_3 model_coordinates; // This is the point which is located left below in the voxelgrid. It contains the real coordinates
    get_extents(vertices, extent, voxel_size, model_coordinates);
    VoxelGrid my_building_grid(extent["x_range_VoxelGrid"], extent["y_range_VoxelGrid"], extent["z_range_VoxelGrid"]);

//    BELOW IS JUST TO PRINT
    for (const auto& entry : objects) {
        if (entry.first < 10) {
            std::cout << "Object ID: " << entry.first << std::endl;
            const Object &obj = entry.second;
            std::cout << "Name: " << obj.name << std::endl;
            std::cout << "Material Type: " << obj.material_type << std::endl;
            std::cout << "Shells:" << std::endl;
            for (const auto &shell: obj.shells) {
                std::cout << "-------New Triangle------" << std::endl;
                std::cout << "\t Modelcoordinates " << "Vertex 1: " << shell.vertex(1) << ", Vertex 2: " << shell.vertex(2)
                          << ", Vertex 3: " << shell.vertex(3) << std::endl;
                std::cout << "\t Voxelgrid integers (xrows, yrows, zrows) ";
                unsigned int x1, y1, z1;
                translate_from_CentroidInModel_to_XYZVoxelgrid(shell.vertex(1), voxel_size, model_coordinates, x1, y1,z1);
                std::cout << "Vertex 1: (" << x1 << ", " << y1 << ", " << z1 << "), ";
                unsigned int x2, y2, z2;
                translate_from_CentroidInModel_to_XYZVoxelgrid(shell.vertex(2), voxel_size, model_coordinates, x2, y2, z2);
                std::cout << "Vertex 2: (" << x2 << ", " << y2 << ", " << z2 << "), ";
                unsigned int x3, y3, z3;
                translate_from_CentroidInModel_to_XYZVoxelgrid(shell.vertex(3), voxel_size, model_coordinates, x3, y3,z3);
                std::cout << "Vertex 3: (" << x3 << ", " << y3 << ", " << z3 << ")" << std::endl;
            }
        }
    }

    return 0;
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
                    current_object.shells.emplace_back(Triangle(vertices[v1], vertices[v2], vertices[v3]));
                }
            }
            else if(line_type == "usemtl"){
                line_stream >> current_object.material_type;
            }
        }
    }
}

void get_extents(std::vector<Point_3> &vertices, std::map<std::string, int> &extent, double &voxel_size, Point_3 &model_coordinates){
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
    model_coordinates = Point_3(x, y, z);
}

Point_3 translate_from_XYZVoxelgrid_to_CentroidInModel(unsigned int &x, unsigned int &y, unsigned int &z,
                                                       double &voxel_size, Point_3 model_coordinates){
    double centroid_x_coordinate = model_coordinates.x() + (x + 0.5) * voxel_size;
    double centroid_y_coordinate = model_coordinates.y() + (y + 0.5) * voxel_size;
    double centroid_z_coordinate = model_coordinates.z() + (z + 0.5) * voxel_size;
    return Point_3(centroid_x_coordinate, centroid_y_coordinate, centroid_z_coordinate);
}

void translate_from_CentroidInModel_to_XYZVoxelgrid(Point_3 pt, double &voxel_size, Point_3 model_coordinates, unsigned int &x, unsigned int &y,
                                                    unsigned int &z){
    x = static_cast<int>(std::floor((pt.x() - model_coordinates.x()) / voxel_size));
    y = static_cast<int>(std::floor((pt.y() - model_coordinates.y()) / voxel_size));
    z = static_cast<int>(std::floor((pt.z() - model_coordinates.z()) / voxel_size));
}