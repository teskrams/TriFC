import ifcopenshell
import numpy as np
from itertools import combinations
import math
import csv
import os
import plotly.graph_objs as go
import matplotlib.pyplot as plt

class IfcPlotter():
    """Class for handling plotting and visualization of elements form IfcOpenShell, the IfcElements class and the IfcTriangle class. 
    Using Plotly and Matplotlib

    file is the path to the ifcfile that you want to load
    """

    def interactive_plot(self, verts_and_faces=[], elements=[], collisions=[], padding=1, check=None, random_color=False, origin_mode=False):
        """Initaites a plotly plot of the elements provided and adds collisions or other personal preferances 
        Using Plotly

        verts_and_faces if you want to use lists with vetices and faces instead of loading the geometry with the IfcElement class
        elements a list with all the elements you want to plot
        collisions boolean list corespoding with the index in elements list, for coloring the triangles
        padding the size of the area around the elements max and min
        check coloring triangles with specific index
        random_color for all triangles
        origin_mode if you use the verts_and_faces option 
        """
        objects = []
        settings = ifcopenshell.geom.settings()
        if origin_mode:
            if len(verts_and_faces) < 2:
                for i in range(len(elements)):
                    shape = ifcopenshell.geom.create_shape(settings, elements[i])
                    grouped_verts = ifcopenshell.util.shape.get_vertices(shape.geometry)
                    x_cords = np.transpose(grouped_verts)[0]
                    y_cords = np.transpose(grouped_verts)[1]
                    z_cords = np.transpose(grouped_verts)[2]
                    faces = ifcopenshell.util.shape.get_faces(shape.geometry)
                    objects.append([x_cords, y_cords, z_cords, faces])
            else:
                x_cords = np.transpose(verts_and_faces[0])[0]
                y_cords = np.transpose(verts_and_faces[0])[1]
                z_cords = np.transpose(verts_and_faces[0])[2]
                faces = verts_and_faces[1]
                objects.append([x_cords, y_cords, z_cords, faces])
        else:
            for i in range(len(elements)):
                objects.append([elements[i].x_coords, elements[i].y_coords, elements[i].z_coords, elements[i].faces])

        fig = go.Figure()

        fig.update_layout(title='Interactive 3D Geometries with Plotly',
                        scene=dict(xaxis_title='X Axis',
                                    yaxis_title='Y Axis',
                                    zaxis_title='Z Axis'))
        
        for i in range(len(objects)):
            x_cords = objects[i][0]
            y_cords = objects[i][1]
            z_cords = objects[i][2]
            faces = objects[i][3]

            if i == 0:
                x_min = min(x_cords) - padding
                x_max = max(x_cords) + padding
                y_min = min(y_cords) - padding
                y_max = max(y_cords) + padding
                z_min = min(z_cords) - padding
                z_max = max(z_cords) + padding

            else:
                x_min_temp = min(x_cords) - padding
                x_max_temp = max(x_cords) + padding
                y_min_temp = min(y_cords) - padding
                y_max_temp = max(y_cords) + padding
                z_min_temp = min(z_cords) - padding
                z_max_temp = max(z_cords) + padding

                if x_min > x_min_temp:
                    x_min = x_min_temp
                if x_max < x_max_temp:
                    x_max = x_max_temp
                if y_min > y_min_temp:
                    y_min = y_min_temp
                if y_max < y_max_temp:
                    y_max = y_max_temp
                if z_min > z_min_temp:
                    z_min = z_min_temp
                if z_max < z_max_temp:
                    z_max = z_max_temp


            for j in range(len(faces)):
                triangle = faces[j]
                x_coords = [x_cords[triangle[0]], x_cords[triangle[1]], x_cords[triangle[2]]]
                y_coords = [y_cords[triangle[0]], y_cords[triangle[1]], y_cords[triangle[2]]]
                z_coords = [z_cords[triangle[0]], z_cords[triangle[1]], z_cords[triangle[2]]]
                
                if random_color:
                    face_color = f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})'
                    opacity = 0.9
                else:
                    face_color = "blue"
                    opacity = 0.9

                if len(collisions) > 0:
                    if collisions[i][j]:
                        face_color = 'red'
                        opacity = 1

                if j == check:
                    face_color = 'black'
                    opacity = 1

                fig.add_trace(go.Mesh3d(
                    x=x_coords,
                    y=y_coords,
                    z=z_coords,
                    i=[0], 
                    j=[1], 
                    k=[2],
                    opacity=opacity,
                    color=face_color,
                    name=j               
                ))

                for k in range(-1, 2):
                    fig.add_trace(go.Scatter3d(
                    x=[x_coords[k], x_coords[k+1]],
                    y=[y_coords[k], y_coords[k+1]],
                    z=[z_coords[k], z_coords[k+1]],
                    mode='lines',
                    line=dict(color="rgb(0,0,0)", width=5),
                    showlegend=False
            ))
                    
        x_range = [x_min, x_max]
        y_range = [y_min, y_max]
        z_range = [z_min, z_max]

        fig.update_layout(scene=dict(
            xaxis=dict(range=x_range),
            yaxis=dict(range=y_range),
            zaxis=dict(range=z_range)
        ))

        fig.update_layout(scene=dict(aspectmode='data'))
        fig.show()

    def plot_triangles(self, triangles, random_color=False):
        """Initaites a plotly plot of all IfcTriangles in the triangles list

        triangles a lisrt of IfcTriangle objects
        random_color for all triangles
        """
        fig = go.Figure()
        for i in range(len(triangles)):
            if random_color:
                color = f'rgb({np.random.randint(0,256)}, {np.random.randint(0,256)}, {np.random.randint(0,256)})'
            else:
                color= "blue"
            fig.add_trace(go.Mesh3d(
                            x=triangles[i].x_coords,
                            y=triangles[i].y_coords,
                            z=triangles[i].z_coords,
                            i=[0, 1, 2],
                            j=[1, 2, 0], 
                            k=[2, 0, 1],  
                            opacity=0.9,
                            color=color,               
                        ))
            for j in range(3): 
                fig.add_trace(go.Scatter3d(
                    x=[triangles[i].x_coords[j], triangles[i].x_coords[(j + 1) % 3]],
                    y=[triangles[i].y_coords[j], triangles[i].y_coords[(j + 1) % 3]], 
                    z=[triangles[i].z_coords[j], triangles[i].z_coords[(j + 1) % 3]],
                    mode='lines',
                    line=dict(color="rgb(0,0,0)", width=5),
                    showlegend=False
                ))

        fig.update_layout(title='Interactive 3D Geometries with Plotly',
                        scene=dict(xaxis_title='X Axis',
                                    yaxis_title='Y Axis',
                                    zaxis_title='Z Axis'))
        
        fig.update_layout(scene=dict(aspectmode='data'))

        fig.show() 

    def plot_triangles_2d(self, triangles): 
        """Initaites a Matplotlib plot of all triangles in the triangles list

        triangles: [[[x, y], [x1, y1], [....]], ..]
        """   
        plt.ion()
        plt.figure()
        for i in range(len(triangles)):
            plt.plot(*zip(*triangles[i]), marker='o', color='r', label='Triangle 1')
            plt.plot(*zip(*triangles[i], triangles[i][0]), color='r')

        plt.xlabel('Axis 1')
        plt.ylabel('Axis 2')
        plt.title('Projected Triangles')
        plt.legend()
        plt.show()

    def plot_intersection_test(self, triangles, points):
        """Initaites a Matplotlib plot of all triangles in the triangles list and points in the points list

        triangles a list of IfcTriangle objects
        points a list of 3d_point
        """
        fig = go.Figure()
        for i in range(len(triangles)):
            color = "blue"
            if i == 0:
                color = "green"
            fig.add_trace(go.Mesh3d(
                            x=triangles[i].x_coords,
                            y=triangles[i].y_coords,
                            z=triangles[i].z_coords,
                            i=[0, 1, 2],
                            j=[1, 2, 0], 
                            k=[2, 0, 1],  
                            opacity=0.9,
                            color=color,               
                        ))
            for j in range(3): 
                fig.add_trace(go.Scatter3d(
                    x=[triangles[i].x_coords[j], triangles[i].x_coords[(j + 1) % 3]],
                    y=[triangles[i].y_coords[j], triangles[i].y_coords[(j + 1) % 3]], 
                    z=[triangles[i].z_coords[j], triangles[i].z_coords[(j + 1) % 3]],
                    mode='lines',
                    line=dict(color="rgb(0,0,0)", width=5),
                    showlegend=False
                ))

        x_list = []
        y_list = []
        z_list = []

        for i in range(len(points)):
            if points[i] is None:
                a = 0
            if np.any(points[i] is None):
                a = 1
            else:
                x_list.append(points[i][0])
                y_list.append(points[i][1])
                z_list.append(points[i][2])

        fig.add_trace(go.Scatter3d(
        x=x_list,
        y=y_list,
        z=z_list,
            ))

        fig.update_layout(title='Interactive 3D Geometries with Plotly',
                        scene=dict(xaxis_title='X Axis',
                                    yaxis_title='Y Axis',
                                    zaxis_title='Z Axis'))

        fig.update_layout(scene=dict(aspectmode='data'))
        fig.show()

class IfcFile():
    def __init__(self, file):
        """Class for handling ifcfiles. Calculates collisions between elements and classifies formwork. 
        Using IfcOpenShell

        file is the path to the ifcfile that you want to load
        """
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