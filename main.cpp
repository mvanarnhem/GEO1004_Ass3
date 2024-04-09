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
#include <CGAL/Surface_mesh.h>
#include <set>
#include "../include/json.hpp"

using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Triangle_3                                              Triangle;
typedef CGAL::Bbox_3                                                    Bbox_3;
typedef Kernel::Vector_3                                                Vector_3;
typedef std::pair<Point_3, Vector_3>                                    Pwn;
typedef CGAL::Surface_mesh<Point_3>                                     Surface_mesh;

struct Point_3_Compare {
    bool operator()(const CGAL::Point_3<Kernel>& p1, const CGAL::Point_3<Kernel>& p2) const {
        if (p1.x() != p2.x())
            return p1.x() < p2.x();
        if (p1.y() != p2.y())
            return p1.y() < p2.y();
        return p1.z() < p2.z();
    }
};

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
void generateOBJ_from_VoxelGrid(VoxelGrid voxelgrid, const std::string& objFilePath, Point_3 origin, double voxel_size, int value);
void intersection(std::map<int, Object> &objects, double &voxel_size, Point_3 &origin, VoxelGrid &my_building_grid);
void label_region(VoxelGrid &voxel_grid, unsigned int label, int start_voxel_index);
void construct_square(const Pwn& point, double side_length, std::vector<Point_3>& vertices, std::vector<std::vector<int>>& faces);
void export_surfaces_as_OBJ(const std::vector<std::vector<Point_3>>& surface, const std::string& filename);
void extract_surfcaes(const VoxelGrid & building_grid, double &voxel_size, Point_3 &origin,
                      std::map<int, std::vector<std::vector<Point_3>>> &surfaces_assigned, unsigned int amount_of_rooms);
unsigned int label_all_regions(VoxelGrid &voxel_grid, unsigned int exterior_label);
void fill_holes_in_wall(VoxelGrid &voxel_grid, double &voxel_size, double threshold_volume);
void remove_parkingLots_from_Wellness(VoxelGrid &voxelgrid);
void exportToCitJSON(std::map<int, std::vector<std::vector<Point_3>>> &surfaces_assigned, std::string outFile);

int main() {
    const std::string input_file = "../data/output_small_house.obj";
    bool process_wellness = false;

    std::map<int, Object> objects;
    std::vector<Point_3> vertices;
    from_OBJ_to_Object(objects, vertices, input_file);

    std::map<std::string, int> extent;
    double voxel_size = 0.3;
    Point_3 origin; // This is the point which is located left below in the voxelgrid. It contains the real coordinates
    get_extents(vertices, extent, voxel_size, origin);

    VoxelGrid my_building_grid(extent["x_range_VoxelGrid"], extent["y_range_VoxelGrid"], extent["z_range_VoxelGrid"]);
    intersection(objects, voxel_size, origin, my_building_grid);
    unsigned int amount_of_rooms = label_all_regions(my_building_grid, 2);

    if (process_wellness){
        remove_parkingLots_from_Wellness(my_building_grid);
    }

    std::map<int, std::vector<std::vector<Point_3>>> surfaces_assigned;
    extract_surfcaes(my_building_grid, voxel_size, origin, surfaces_assigned, amount_of_rooms);

    exportToCitJSON(surfaces_assigned, "mybuilding.city.json");

    return 0;
}

Point_3 scale_point(Point_3 &vertex, json &j){
    double x = (vertex.x() - j["transform"]["translate"][0].get<double>())/
               (j["transform"]["scale"][0].get<double>());
    double y = (vertex.y() - j["transform"]["translate"][1].get<double>())/
               (j["transform"]["scale"][1].get<double>());
    double z = (vertex.z() - j["transform"]["translate"][2].get<double>())/
               (j["transform"]["scale"][2].get<double>());
    return Point_3(x, y, z);
}

void exportToCitJSON(std::map<int, std::vector<std::vector<Point_3>>> &surfaces_assigned, std::string outFile){
    nlohmann::json j;
    j["type"] = "CityJSON";
    j["version"] = "2.0";
    j["transform"] = nlohmann::json::object();
    j["transform"]["scale"] = nlohmann::json::array({1.0, 1.0, 1.0});
    j["transform"]["translate"] = nlohmann::json::array({0.0, 0.0, 0.0});
    j["CityObjects"] = nlohmann::json::object();
    j["vertices"] = json::array();

    for (const auto& entry : surfaces_assigned) {
        json room_json = json::array();
        int room_id = entry.first;
        json face_json = json::array();
        auto room_faces = entry.second;
        for (auto &face: room_faces) {
            face_json.clear();
            for (auto &pt: face) {
                Point_3 scaled_Point = scale_point(pt, j);
                int foundIndex = -1;
                for (size_t i = 0; i < j["vertices"].size(); ++i) {
                    auto &vertex = j["vertices"][i];
                    if (std::fabs(vertex[0].get<double>() - scaled_Point.x()) < 0.0001 &&
                        std::fabs(vertex[1].get<double>() - scaled_Point.y()) < 0.0001 &&
                        std::fabs(vertex[2].get<double>() - scaled_Point.z()) < 0.0001) {
                        foundIndex = i;
                        break;
                    }
                }
                if (foundIndex == -1) {
                    foundIndex = j["vertices"].size();
                    j["vertices"].push_back({scaled_Point.x(), scaled_Point.y(), scaled_Point.z()});
                }
                face_json.push_back(foundIndex);
            }
            room_json.push_back(json::array({face_json}));
        }

        if (room_id == 2) {
            j["CityObjects"]["ExteriorBuilding"]["type"]= "Building";
            j["CityObjects"]["ExteriorBuilding"]["children"] = json::array();
            j["CityObjects"]["ExteriorBuilding"]["geometry"] = json::array({
                                                                                   {
                                                                                           { "type", "MultiSurface"},
                                                                                           { "lod", "2.0"},
                                                                                           { "boundaries", room_json}
                                                                                   }});
        }

        else { // not 2 means its a room
            std::string room_name = "BuildingRoom";
            room_name += std::to_string(room_id - 2);
            j["CityObjects"]["ExteriorBuilding"]["children"].push_back(room_name);
            j["CityObjects"][room_name]["parents"]=json::array({"ExteriorBuilding"});
            j["CityObjects"][room_name]["type"]= "BuildingRoom";
            j["CityObjects"][room_name]["geometry"] = json::array({
                                                                          {
                                                                                  { "type", "MultiSurface"},
                                                                                  { "lod", "2.0"},
                                                                                  { "boundaries", room_json}
                                                                          }
                                                                  });
        }

    }

    std::string json_string = j.dump();
    std::ofstream out_stream(outFile);
    if (out_stream.is_open()) {
        out_stream << json_string;
        out_stream.close();
        std::cout << "CityJSON file generated: mybuilding.city.json" << std::endl;
    } else {
        std::cerr << "Error: Unable to open file for writing." << std::endl;
    }
}

void extract_surfcaes(const VoxelGrid & building_grid, double &voxel_size,  Point_3 &origin,
                      std::map<int, std::vector<std::vector<Point_3>>> &surfaces_assigned, unsigned int amount_of_rooms){
    std::map<int, std::vector<Pwn>> boundary_points_assigned;
    for (int i = 2; i < (2 + amount_of_rooms); ++i) {
        boundary_points_assigned[i] = std::vector<Pwn>();
    }
    unsigned int index_voxel = 0;
    for (auto const &voxel: building_grid.voxels) {
        unsigned int x = index_voxel % building_grid.max_x;
        unsigned int y = (index_voxel / building_grid.max_x) % building_grid.max_y;
        unsigned int z = index_voxel / (building_grid.max_x * building_grid.max_y);
        if (voxel == 1) {
            for (int dz = -1; dz <= 1; ++dz) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0)
                            continue;
                        if (abs(dx) + abs(dy) + abs(dz) != 1)
                            continue;
                        int nx = x + dx;
                        int ny = y + dy;
                        int nz = z + dz;
                        Point_3 centre_point = translate_VoxelGridVoxel_to_RealWorldCoordinates(x, y, z, voxel_size, origin);
                        if (building_grid(nx, ny, nz) != 1) {
                            double x_boundary = centre_point.x() + 0.5 * voxel_size * dx;
                            double y_boundary = centre_point.y() + 0.5 * voxel_size * dy;
                            double z_boundary = centre_point.z() + 0.5 * voxel_size * dz;
                            int neighbour_value = building_grid(nx, ny, nz);
                            boundary_points_assigned[neighbour_value].emplace_back(std::make_pair(
                                    Point_3(x_boundary, y_boundary, z_boundary),
                                    Vector_3(-dx, -dy, -dz)));
                        }
                    }
                }
            }
        }
        index_voxel += 1;
    }

    for (const auto& entry : boundary_points_assigned) {
        int key = entry.first;
        const auto& points = entry.second;
        std::vector<std::vector<Point_3>> current_surfaces;
        for (const auto& pwn : points) {
            std::vector<Point_3> vertices;
            std::vector<std::vector<int>> faces;
            construct_square(pwn, voxel_size, vertices, faces);
            current_surfaces.push_back(vertices);
        }
        surfaces_assigned[key] = current_surfaces;
    }
}


void export_surfaces_as_OBJ(const std::vector<std::vector<Point_3>>& surface, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
        return;
    }
    for (const auto& polygon : surface) {
        for (const auto& vertex : polygon) {
            outFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
        }
    }
    int offset = 1;
    for (size_t i = 0; i < surface.size(); ++i) {
        outFile << "f";
        for (size_t j = 0; j < surface[i].size(); ++j) {
            outFile << " " << j + offset;
        }
        offset += surface[i].size();
        outFile << std::endl;
    }
    outFile.close();
}


Vector_3 orthogonal_vector(const Vector_3& v) {
    if (abs(v.z()) == 1){
        return Vector_3(0, -v.z(), v.y());
    }
    else {
        return Vector_3(-v.y(), v.x(), 0);
    }
}

void construct_square(const Pwn& point, double side_length, std::vector<Point_3>& vertices, std::vector<std::vector<int>>& faces) {
    const Point_3& center = point.first;
    const Vector_3& normal = point.second;
    Vector_3 u = orthogonal_vector(normal);
    Vector_3 v = CGAL::cross_product(normal, u);
    double half_length = side_length / 2.0;
    std::vector<Point_3> square_vertices;
    square_vertices.push_back(center + half_length * u + half_length * v);
    square_vertices.push_back(center + half_length * u - half_length * v);
    square_vertices.push_back(center - half_length * u - half_length * v);
    square_vertices.push_back(center - half_length * u + half_length * v);
    for (const auto& vertex : square_vertices) {
        vertices.push_back(vertex);
    }
    int offset = vertices.size() - square_vertices.size();
    faces.push_back({offset, offset + 1, offset + 2});
    faces.push_back({offset + 2, offset + 3, offset});
}

void fill_holes_in_wall(VoxelGrid &voxel_grid, double &voxel_size, double threshold_volume) {
    std::map<unsigned int, int> region_size_map; // initiate a map to keep track of how many voxels per label
    for (unsigned int i = 0; i < voxel_grid.voxels.size(); ++i) {
        auto &voxel_value = voxel_grid.voxels[i];
        region_size_map[voxel_value]++;
    }
    std::set<int> holes_to_fill;
    for (const auto& pair : region_size_map){
        double region_volume = pair.second * std::pow(voxel_size, 3);
        if (region_volume < threshold_volume){
            holes_to_fill.insert(pair.first);
        }
    }

    for (unsigned int i = 0; i < voxel_grid.voxels.size(); ++i) {
        auto &voxel_value = voxel_grid.voxels[i];

        if (holes_to_fill.find(voxel_value) != holes_to_fill.end()){
            // std::cout<< "overwriting voxel id=" << i << "\t" << voxel_value << "--> 1" << std::endl;
            voxel_grid.voxels[i] = 1;
        }
    }
}

unsigned int label_all_regions(VoxelGrid &voxel_grid, unsigned int exterior_label) {
    std::cout << "Labelling spaces" << "-------------------" << std::endl;

    unsigned int label_number = exterior_label;

    for (unsigned int i = 0; i < voxel_grid.voxels.size(); ++i) {
        auto &voxel_value = voxel_grid.voxels[i];
        if (voxel_value == 0) {
            label_region(voxel_grid, label_number, i);
            ++label_number;
        }
    }
    std::cout << "-----------------------------------" << std::endl;
    unsigned int amount_of_rooms = label_number - exterior_label;
    return amount_of_rooms;
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

    to_visit.emplace_back(start_voxel_index); // place start voxel in to_visit
    voxel_grid(start_x, start_y, start_z) = label; // set it's value to this label
    int voxels_labelled = 0; // keep track of how many voxels are added to this label

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
        to_visit.erase(to_visit.begin());
        ++voxels_labelled;
    }
    std::cout << "\tlabel: " << label << "\tstarting at:\t" <<
              "(" << start_x << ", " << start_y << ", " << start_z << ")" << "\t room_size= " << voxels_labelled << " voxels" << std::endl;
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

void generateOBJ_from_VoxelGrid(VoxelGrid voxel_grid, const std::string& objFilePath, Point_3 origin, double voxel_size, int value) {
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

void remove_parkingLots_from_Wellness(VoxelGrid &voxelgrid){
    // First step: the red parts are removed, the parts that only contain 3 voxels assigned to 1 in every z-direction
    for (unsigned int x_fixed = 0; x_fixed < voxelgrid.max_x; ++x_fixed) {
        for (unsigned int y_fixed = 0; y_fixed < voxelgrid.max_y; ++y_fixed) {
            unsigned int onesCount = 0;
            bool includes_interior = false;
            std::vector<unsigned int> walls_z_index;
            for (unsigned int z = 0; z < voxelgrid.max_z; ++z) {
                if (voxelgrid(x_fixed, y_fixed, z) == 1) {
                    ++onesCount;
                    walls_z_index.emplace_back(z);
                }
                else if (voxelgrid(x_fixed, y_fixed, z) > 2){
                    includes_interior = true;
                }
            }
            if (not includes_interior && onesCount <= 3) {
                for (auto z: walls_z_index){
                    voxelgrid(x_fixed, y_fixed, z) = 2;
                }
            }
        }
    }
    // Second step: the yellow parts are removed, two voxels with value 1 in z-direction and with the lowest z-value
    // that have a lot of exterior voxels above them.
    for (unsigned int x_fixed = 0; x_fixed < voxelgrid.max_x; ++x_fixed) {
        for (unsigned int y_fixed = 0; y_fixed < voxelgrid.max_y; ++y_fixed) {
            std::vector<unsigned int> change_voxels_zIndex;
            std::vector<unsigned int> previous_ones;
            std::vector<unsigned int> keep_track_of_exterior;
            bool not_crossed_interior = true;
            bool removed_bottom = false;
            for (unsigned int z = 0; z < voxelgrid.max_z; ++z) {
                if (not removed_bottom) {
                    if (voxelgrid(x_fixed, y_fixed, z) == 1) {
                        if (keep_track_of_exterior.size() > 10 && not_crossed_interior) {
                            change_voxels_zIndex.insert(change_voxels_zIndex.end(), previous_ones.begin(),
                                                        previous_ones.end());
                            previous_ones.clear();
                            keep_track_of_exterior.clear();
                            removed_bottom = true;
                        }
                        previous_ones.emplace_back(z);
                    } else if (voxelgrid(x_fixed, y_fixed, z) == 2) {
                        if (previous_ones.size() == 2) {
                            keep_track_of_exterior.emplace_back(z);
                        } else {
                            previous_ones.clear();
                            keep_track_of_exterior.clear();
                        }
                    } else {
                        not_crossed_interior = false;
                        previous_ones.clear();
                        keep_track_of_exterior.clear();
                    }
                }
            }
            for (unsigned int z: change_voxels_zIndex) {
                voxelgrid(x_fixed, y_fixed, z) = 2;
            }
        }
    }
}