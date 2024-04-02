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
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3                                     Plane;
typedef Kernel::Point_3                                     Point_3;
typedef Kernel::Triangle_3                                  Triangle;
typedef CGAL::Bbox_3                                        Bbox_3;
typedef Kernel::Iso_cuboid_3                                Iso_cuboid_3;


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
void generateOBJ(VoxelGrid voxelgrid, const std::string& objFilePath, Point_3 origin, double voxel_size, bool visualize_walls);



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

    generateOBJ(my_building_grid, "output_all_points.obj", origin, voxel_size, false);

    std::cout << "Origin: " << origin.x() << ", " << origin.y() << ", " << origin.z() << ")" << std::endl;

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
//                            std::cout << "The triangle and the bounding box intersect." << std::endl;
                        }
                    }
                }
            }
        }
    }
    generateOBJ(my_building_grid, "output_only_0.obj", origin, voxel_size, false);
    generateOBJ(my_building_grid, "output_only_1.obj", origin, voxel_size, true);

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

void generateOBJ(VoxelGrid voxelgrid, const std::string& objFilePath, Point_3 origin, double voxel_size, bool visualize_walls) {
    std::ofstream objFile(objFilePath);
    int value;
    if (visualize_walls){
        value = 1;
    }
    else if (not visualize_walls){
        value = 0;
    }

    if (!objFile.is_open()) {
        std::cerr << "Error: Unable to open file: " << objFilePath << std::endl;
        return;
    }

    // Define material library
    objFile << "mtllib colors.mtl\n\n";

    for (unsigned int x = 0; x < voxelgrid.max_x; ++x) {
        for (unsigned int y = 0; y < voxelgrid.max_y; ++y) {
            for (unsigned int z = 0; z < voxelgrid.max_z; ++z) {
                // Check if voxel value is 1
                if (voxelgrid(x, y, z) == value) {
                    Point_3 pt = translate_VoxelGridVoxel_to_RealWorldCoordinates(x, y, z, voxel_size, origin);
                    // Write vertex coordinates
                    objFile << "v " << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
                }
            }
        }
    }

    objFile.close();
}