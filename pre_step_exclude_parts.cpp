#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             Kernel;
typedef Kernel::Point_3                                                 Point_3;
typedef Kernel::Triangle_3                                              Triangle_3;

struct Object {
    std::string name;
    std::vector<Triangle_3> shells;
    std::string material_type;
};

void wrong_groups_from_DefaultMaterial(const std::string &input_file, std::vector<std::string> &remove_from_default);
void make_objects(const std::string &input_file, const std::vector<std::string> &groups_to_remove,
                   const std::vector<std::string> &usemtl_to_remove, std::map<int, Object> &objects,
                   std::vector<Point_3> &vertices);
void exportToObj(const std::map<int, Object>& objects, std::ofstream& outFile);

int main() {
    const std::string input_file = "../data/output_wellness.obj";

    std::vector<std::string> remove_from_default;
    // The following function will get the group names of the groups in Default which need to be removed
    // See the picture I have send you
    wrong_groups_from_DefaultMaterial(input_file, remove_from_default);

    std::vector<std::string> remove_from_floor = {"product-1889ce11-0dd2-4ce2-9da4-88b9677a679a-body",
                                                 "product-236da791-eb8e-4871-893b-88122d4465b8-body",
                                                 "product-1889ce11-0dd2-4ce2-9da4-88b9677a67f5-body"};
    std::vector<std::string> groups_to_remove;
    groups_to_remove.insert(groups_to_remove.end(), remove_from_floor.begin(), remove_from_floor.end());
    groups_to_remove.insert(groups_to_remove.end(), remove_from_default.begin(), remove_from_default.end());

    std::vector<std::string> usemtl_to_remove = {"surface-style-652-asphalt,-bitumen", "IfcSlab", "IfcRailing",
                                                 "IfcRamp", "IfcStairFlight", "IfcBuildingElementProxy",
                                                 "IfcFurnishingElement", "IfcFlowTerminal",
                                                 "surface-style-75901-metal---stainless-steel",
                                                 "IfcMember"};

    std::map<int, Object> objects;
    std::vector<Point_3> vertices;
    make_objects(input_file, groups_to_remove, usemtl_to_remove, objects, vertices);
    std::ofstream outFile("output_wellness_2.obj");
    if (outFile.is_open()) {
        exportToObj(objects, outFile);
        outFile.close();
        std::cout << "Objects saved to output.obj\n";
    } else {
        std::cerr << "Unable to create output.obj\n";
    }

}

void make_objects(const std::string &input_file, const std::vector<std::string> &groups_to_remove,
                   const std::vector<std::string> &usemtl_to_remove, std::map<int, Object> &objects,
                  std::vector<Point_3> &vertices){
        std::ifstream input_stream;
        input_stream.open(input_file);
        bool keep_faces = true;

        if (input_stream.is_open()) {
            std::string line;
            Object current_object;
            int index_object = -1;

            std::string current_group;
            std::string current_usemtl;

            while (getline(input_stream, line)) {
                std::istringstream line_stream(line);
                std::string line_type;
                line_stream >> line_type;
                if (line_type == "usemtl") {
                    line_stream >> current_usemtl;
                    if (std::count(usemtl_to_remove.begin(), usemtl_to_remove.end(), current_usemtl) != 0) {
                        keep_faces = false;
                    }
                    else{
                        current_object.material_type = current_usemtl;
                    }
                }
                else if (line_type == "v"){
                    double x, y, z;
                    line_stream >> x >> y >> z;
                    Point_3 vert = Point_3(x, y, z);
                    vertices.emplace_back(vert);
                }
                else if (line_type == "g") {
                    line_stream >> current_group;
                    if (index_object != -1) {
                        if (keep_faces) {
                            objects[index_object] = current_object;
                        }
                        Object current_object;
                    }
                    index_object += 1;
                    current_object.name = current_group;
                    keep_faces = true;
                    if (std::count(groups_to_remove.begin(), groups_to_remove.end(), current_group) != 0) {
                        keep_faces = false;
                    }
                }
                else if(line_type == "f"){
                    if (keep_faces) {
                        std::string v1_str, v2_str, v3_str;
                        int v1, v2, v3;
                        line_stream >> v1_str >> v2_str >> v3_str;

                        v1 = stoi(v1_str.substr(0, v1_str.find("//")));
                        v2 = stoi(v2_str.substr(0, v2_str.find("//")));
                        v3 = stoi(v3_str.substr(0, v3_str.find("//")));

                        std::string remaining_part;
                        line_stream >> remaining_part;
                        if (!remaining_part.empty()) {
                            std::cout << "The input face contains more than 3 vertices." << std::endl;
                        } else {
                            Triangle_3 triangle = Triangle_3(vertices[v1 - 1], vertices[v2 - 1], vertices[v3 - 1]);
                            current_object.shells.emplace_back(triangle);
                        }
                    }
                }
            }
            objects[index_object] = current_object;
        }
};

void exportToObj(const std::map<int, Object>& objects, std::ofstream& outFile) {
    std::map<Point_3, int> vertexIndices;
    for (const auto& entry : objects) {
        const Object& obj = entry.second;
        outFile << "g " << obj.name << "\n";
        for (const Triangle_3& triangle : obj.shells) {
            for (int i = 0; i < 3; ++i) {
                const Point_3& vertex = triangle.vertex(i);
                auto it = vertexIndices.find(vertex);
                if (it == vertexIndices.end()) {
                    int index = vertexIndices.size() + 1;
                    vertexIndices[vertex] = index;
                    outFile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
                }
            }
            outFile << "f ";
            for (int i = 0; i < 3; ++i) {
                const Point_3& vertex = triangle.vertex(i);
                outFile << vertexIndices[vertex] << " ";
            }
            outFile << "\n";
        }
        outFile << "usemtl " << obj.material_type << "\n";
    }
}

void wrong_groups_from_DefaultMaterial(const std::string &input_file, std::vector<std::string> &remove_from_default){
    // This function will save the names of the groups in Default that have more than 12 faces
    std::ifstream input_stream;
    input_stream.open(input_file);
    if (input_stream.is_open()) {
    std::string line;
    std::string current_group;
    std::string current_usemtl;
    int amount_of_faces = 0;
    while (getline(input_stream, line)) {
        std::istringstream line_stream(line);
        std::string line_type;
        line_stream >> line_type;
        if (line_type == "usemtl"){
            line_stream >> current_usemtl;
        }
        if (current_usemtl == "DefaultMaterial"){
            if (line_type == "f"){
                amount_of_faces += 1;
            }
            else if (line_type == "g"){
                if (amount_of_faces > 12){
                    remove_from_default.emplace_back(current_group);
                }
                amount_of_faces = 0;
            }
        }
        if (line_type == "g"){
            line_stream >> current_group;
            }
        }
    }
}