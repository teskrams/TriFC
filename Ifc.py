import ifcopenshell
import numpy as np
from itertools import combinations
import math
import csv
import os

class IfcFile():
    def __init__(self, file):
        self.settings = ifcopenshell.geom.settings()
        self.model = ifcopenshell.open(file)
        self.elements = None

    def material_kvalitet_is_betong(self, string, debug=False):
        """Return a boolean statement wheter or not the name idicates concrete

        debug option will anable the function to print
        """
        if (string.startswith("B1") 
            or string.startswith("B2") 
            or string.startswith("B3") 
            or string.startswith("B4") 
            or string.startswith("B5") 
            or string.startswith("B6") 
            or string.startswith("B7")) and not string.startswith("B500NC"):
            if debug:
                print("materialkvalitet string", string)
                print("Element is made of concrete")
            return True
        return False
    
    def set_elements(self, elements=[], printing=1, debug=False, academic=False):
        """Sets the elements in the ifcfile based on the file or a set of elements.

        debug option will anable the function to print
        printing option adjusts the amount of prining
        academic option for files with no extra information
        """
        finale_elements = []
        if len(elements) > 0:
            for i in range(len(elements)):
                if academic:
                    finale_elements.append(elements[i])
                    continue
                material_kvalitet = elements[i].info.get("01 Merknader").get("K01 Materialkvalitet")
                net_volume = float(elements[i].info.get("BaseQuantities").get("NetVolume"))
                if net_volume < 0.05:
                    continue
                if self.material_kvalitet_is_betong(material_kvalitet) or ("AVRETTING" in elements[i].name) or ("LAGER" in elements[i].name) or ("SPRENGNINGSNIVÅ" in elements[i].name) or ("FORSKALING" in elements[i].name):
                    finale_elements.append(elements[i])
            self.elements = finale_elements

        else:
            # Dette vil bare fungere på kvithammar-åsen
            # Hopper over alle som blir feil, kan være dårlig løsning, gjelder stort sett små elementer
            self.elements = []
            for storey in self.model.by_type("IfcBuildingStorey"):
                storey_elements_needed = []
                all_storey_elements = ifcopenshell.util.element.get_decomposition(storey)
                if debug or printing > 2:
                    print(f"There are {len(all_storey_elements)} elements located on storey {storey.Name}")
                for i in range(len(all_storey_elements)):
                    material_kvalitet = "None"
                    net_volume = 0
                    info = ifcopenshell.util.element.get_psets(all_storey_elements[i])
                    if academic:
                        storey_elements_needed.append(self.get_element(all_storey_elements[i].GlobalId, debug=debug, academic=academic))
                        if debug or printing > 3:
                            print("element ADDED")
                        continue
                    try:
                        #01 Merknader or E&KAA is the pset names used in the Kvithammar-Asen project
                        material_kvalitet = info.get("01 Merknader").get("K01 Materialkvalitet")
                    except:
                        if debug or printing > 3:
                            print("could not fetch 01 merknader")
                    try:
                        #01 Merknader or E&KAA is the pset names used in the Kvithammar-Asen project 
                        material_kvalitet = info.get("E6KAA").get("K01 Materialkvalitet")
                    except:
                        if debug or printing > 3:
                            print("could not fetch E6KAA")
                    try:
                        net_volume = float(info.get("BaseQuantities").get("NetVolume"))
                        if net_volume < 0.05:
                            if debug or printing > 3:
                                print("net_volume < 0.05")
                            continue
                    except:
                        if debug or printing > 3:
                            print("could not fetch basequantities")
                    try:
                        net_volume = float(info.get("Calculated Geometry Values").get("Volume"))
                        if net_volume < 0.05:
                            if debug or printing > 3:
                                print("net_volume < 0.05")
                            continue
                    except:
                        if debug or printing > 3:
                            print("could not fetch Calculated Geometry Values")

                    name = all_storey_elements[i].Name
                    if debug or printing > 3:
                        print("name:", name)
                    try:
                        if self.material_kvalitet_is_betong(material_kvalitet, debug=debug) or ("AVRETTING" in name) or ("LAGER" in name) or ("SPRENGNINGSNIVÅ" in name) or ("FORSKALING" in name):
                            if debug or printing > 3:
                                print("element NEEDED")
                                print("id:", all_storey_elements[i].GlobalId)
                            storey_elements_needed.append(self.get_element(all_storey_elements[i].GlobalId, debug=debug))
                            if debug or printing > 3:
                                print("element ADDED")
                    except:
                        continue
                self.elements = self.elements + storey_elements_needed
            if debug or printing > 1:
                print("Number of elements found in file:", len(self.elements))


    def get_element(self, guid, need_formwork=True, simplify=True, academic=False, printing=1, debug=False):
        """Fetches an element from a ifcfile based on the guid of the element

        debug option will anable the function to print
        printing option adjusts the amount of prining
        academic option for files with no extra information
        need_formwork option for manual tell the model to formwork or not
        simplify option for anabling simplifying of the element mesh
        """
        temp_element = self.model.by_guid(guid)
        shape = ifcopenshell.geom.create_shape(self.settings, temp_element)
        vertices = ifcopenshell.util.shape.get_vertices(shape.geometry)
        edges = ifcopenshell.util.shape.get_edges(shape.geometry)
        faces = ifcopenshell.util.shape.get_faces(shape.geometry)
        info = ifcopenshell.util.element.get_psets(temp_element)
        name = temp_element.Name
        edges_vertices = [[vertices[edge[0]], vertices[edge[1]]] for edge in edges]

        if academic:
            info = None
        element = IfcElement(faces, edges_vertices, vertices, guid=guid, name=name, info=info, need_formwork=need_formwork)

        if len(element.faces) > 200 and simplify:
            element.simplify_trimesh_element(target_face_count=200, debug=debug)

        return element

    def calculate_collisions(self, printing=1, debug=False, compute_new_triangles=False, compute_new_triangles_2d=False):
        """Calculaters all collisions between all elements in the elements list

        debug option will anable the function to print
        printing option adjusts the amount of prining
        compute_new_triangles option tells wheter or not new triangles in 3d should be computed
        compute_new_triangles_2d option tells wheter or not new triangles in 2d should be computed
        """
        for combo in combinations(self.elements, 2):
            test_needed = combo[0].need_collision_test(combo[1], debug=debug)
            if test_needed:
                if printing > 1:
                    print("combo:", combo[0].name, "and", combo[1].name) #Burde vel egentlig ikke være nødvendig med dobbel her? Blir litt mer med bare en
                combo[0].collisions_two_elements(combo[1], debug=debug, compute_new_triangles=compute_new_triangles, compute_new_triangles_2d=compute_new_triangles_2d)
                combo[1].collisions_two_elements(combo[0], debug=debug, compute_new_triangles=compute_new_triangles, compute_new_triangles_2d=compute_new_triangles_2d)
        if printing > 2:
            print("Updating verts and faces after new triangles in collisions")
        
        for element in self.elements:
            element.update_verts_and_faces(collisions=True)

    def calculate_formwork(self, printing=1, debug=False):
        """Classifies all triagles in each element for formwork on all elements in the elements list

        debug option will anable the function to print
        printing option adjusts the amount of prining
        """
        def update_triangle(element, triangle1_def):
            point11 = [element.x_coords[triangle1_def[0]], element.y_coords[triangle1_def[0]], element.z_coords[triangle1_def[0]]]
            point12 = [element.x_coords[triangle1_def[1]], element.y_coords[triangle1_def[1]], element.z_coords[triangle1_def[1]]]
            point13 = [element.x_coords[triangle1_def[2]], element.y_coords[triangle1_def[2]], element.z_coords[triangle1_def[2]]]
            triangle1 = IfcTriangle(point11, point12, point13)
            return triangle1

        for j in range(len(self.elements)):
            element = self.elements[j]
            element.formwork = [False] * len(element.faces)
            element.formwork_info = [None] * len(element.faces)

            if not element.need_formwork:
                continue

            for i in range(len(element.faces)):
                triangle1_def = element.faces[i]
                triangle1 = update_triangle(element, triangle1_def)
                # tried updating some here
                if triangle1.is_normal_vector_pointing_out(element):
                    normal = triangle1.normal_vector
                else:
                    normal = triangle1.normal_vector * -1

                degree = math.degrees(np.arcsin(normal[2]))
                collision = element.collisions[i]
                inside, inside_element = element.triangle_inside_element_in_list(triangle1, self.elements)

                if inside:
                    inside = triangle1.inside_element_v6(inside_element) # skal være v5

                decicion_variables = (degree, collision, inside, inside_element)
                element.formwork_info[i] = decicion_variables
                if debug or printing > 3:
                    print("(j, i)", (j, i))
                    print("(degree, collision, inside)", (degree, collision, inside))

                element.formwork[i] = triangle1.formwork_needed(degree, collision, inside, element, inside_element, debug=False)

    def update_formwork_area(self):
        """Updates the total formwork on each element in the elements list based on the formwork classification"""
        for i in range(len(self.elements)):
            self.elements[i].update_formwork_area()

    def make_result_csv(self, name=None, path=None, printing=1):
        """Returns a csv with the results of the formwork classification for all elements in the elements list

        printing option adjusts the amount of prining
        path option telles where the csv will be saved
        name option are the name of the result file
        """
        result = []
        for i in range(len(self.elements)):
            element = self.elements[i]
            element_result = [element.guid, element.name, element.need_formwork, element.get_surface_area(), element.formwork_area]
            result.append(element_result)

        file_path = "/Users/theodor/kode/master/master/results"

        if path is not None:
            file_path = path

        file_name = "results.csv"

        if name is not None:
            file_name = name

        result_filepath = os.path.join(file_path, file_name)

        with open(result_filepath, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["GUID", "Name", "Need_formwork", "Total_surface_area", "Formwork_area"])
            for line in result:
                writer.writerow(line) 

    def run_full_calculation(self, printing=1, name=None, compute_new_triangles_2d=False, compute_new_triangles=True, result_filepath=None, academic=False):
        """Runs the whole model on the file and makes the csv with results

        printing option adjusts the amount of prining
        result_filepath option telles where the csv will be saved
        name option are the name of the result file
        compute_new_triangles option tells wheter or not new triangles in 3d should be computed
        compute_new_triangles_2d option tells wheter or not new triangles in 2d should be computed
        academic option for running the model only based on the geometry
        """
        if printing > 0:
            print("Calculation started")

        self.set_elements(printing=printing, academic=academic)

        if printing > 0:
            print("All elements loaded")

        self.calculate_collisions(debug=False, printing=printing, compute_new_triangles=compute_new_triangles, compute_new_triangles_2d=compute_new_triangles_2d)

        if printing > 0:
            print("All collisions calculated")

        self.calculate_formwork(debug=False, printing=printing)

        if printing > 0:
            print("All triangles formwork classified")

        self.update_formwork_area()

        if printing > 0:
            print("Formwork area calculated for all elements")

        self.make_result_csv(name=name, path=result_filepath)

        if printing > 0:
            print("Results csv made")