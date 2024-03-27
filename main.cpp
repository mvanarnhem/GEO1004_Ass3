#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Plane_3                      Plane;
typedef Kernel::Point_3                      Point_3;
typedef Kernel::Triangle_3                   Triangle;

struct Object {
    std::string name;
    std::vector<Triangle> shells;
    std::string material_type;
};

int main() {
    const std::string input_file = "../data/output_small_house.obj";
    std::map<int, Object> objects;
    std::vector<Point_3> vertices;

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

        for (const auto& entry : objects) {
            std::cout << "Object ID: " << entry.first << std::endl;
            const Object& obj = entry.second;
            std::cout << "Name: " << obj.name << std::endl;
            std::cout << "Material Type: " << obj.material_type << std::endl;
            std::cout << "Shells:" << std::endl;
            for (const auto& shell : obj.shells) {
                std::cout << "\t" <<"  Vertex 1: " << shell.vertex(1) << ", Vertex 2: " << shell.vertex(2) << ", Vertex 3: " << shell.vertex(3) << std::endl;
            }
            std::cout << std::endl;
        }
    }
        return 0;
}