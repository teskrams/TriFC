import ifcopenshell
import numpy as np
from itertools import combinations
import math
import csv
import os
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import trimesh


#########################################################################################################################
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



#########################################################################################################################
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



#########################################################################################################################
class IfcElement():
    def __init__(self, faces, edges, verticies, name=None, guid=None, info=None, need_formwork=True):
        """The IfcElement class saves elements from the IFC-file strucutue. The class saves information on formwork and collisions.

        faces from IfcOpenShells loading of the geometry
        edges from IfcOpenShells loading of the geometry
        vertices from IfcOpenShells loading of the geometry
        name from IfcOpenShell or of own preferance
        guid from IfcOpenShells loading of the geometry
        need_formwork if you want to manual tell if the element should have formwork
        info from the pset, containg all element information in the IFC-file
        """
        self.name = name
        self.guid = guid
        self.info = info
        if info is not None:
            self.need_formwork = self.check_need_formwork(info)
        else:
            self.need_formwork = need_formwork
        self.faces = faces
        self.faces_copy = faces
        self.faces_org = faces
        self.edges = edges
        self.edges_copy = edges
        self.verticies = verticies
        self.verticies_copy = verticies
        self.faces_changed = [False] * len(self.faces)

        self.x_coords = np.transpose(verticies)[0]
        self.y_coords = np.transpose(verticies)[1]
        self.z_coords = np.transpose(verticies)[2]
        self.x_min = min(self.x_coords)
        self.x_max = max(self.x_coords)
        self.y_min = min(self.y_coords)
        self.y_max = max(self.y_coords)
        self.z_min = min(self.z_coords)
        self.z_max = max(self.z_coords)

        self.x_coords_copy = np.transpose(verticies)[0]
        self.y_coords_copy = np.transpose(verticies)[1]
        self.z_coords_copy = np.transpose(verticies)[2]
        self.x_min_copy = min(self.x_coords_copy)
        self.x_max_copy = max(self.x_coords_copy)
        self.y_min_copy = min(self.y_coords_copy)
        self.y_max_copy = max(self.y_coords_copy)
        self.z_min_copy = min(self.z_coords_copy)
        self.z_max_copy = max(self.z_coords_copy)

        self.x_coords_org = np.transpose(verticies)[0]
        self.y_coords_org = np.transpose(verticies)[1]
        self.z_coords_org = np.transpose(verticies)[2]
        self.x_min_org = min(self.x_coords_org)
        self.x_max_org = max(self.x_coords_org)
        self.y_min_org = min(self.y_coords_org)
        self.y_max_org = max(self.y_coords_org)
        self.z_min_org = min(self.z_coords_org)
        self.z_max_org = max(self.z_coords_org)

        self.collisions = [False] * len(faces)
        self.collisions_copy = [False] * len(faces)
        self.formwork = [False] * len(faces)
        self.formwork_info = None
        self.formwork_area = 0
        self.triangle_areas = []

    def get_surface_area(self):
        """Returns the size of the outer surface of the element"""
        surface_area = 0
        for i in range(len(self.faces)):
            surface_area = surface_area + self.get_triangle(i).area()
        return surface_area
    
    def get_triangle(self, id):
        """Returns the IfcTriangle defined by the vertices and face on the given index"""
        triangle_def = self.faces[id]
        point1 = [self.x_coords[triangle_def[0]], self.y_coords[triangle_def[0]], self.z_coords[triangle_def[0]]]
        point2 = [self.x_coords[triangle_def[1]], self.y_coords[triangle_def[1]], self.z_coords[triangle_def[1]]]
        point3 = [self.x_coords[triangle_def[2]], self.y_coords[triangle_def[2]], self.z_coords[triangle_def[2]]]
        triangle = IfcTriangle(point1, point2, point3)

        return triangle 

    def get_org_triangle(self, id):
        """Returns the orgininal IfcTriangle defined by the vertices and face on the given index needs index from the orginal lists"""
        triangle_def = self.faces[id]
        point1 = [self.x_coords[triangle_def[0]], self.y_coords[triangle_def[0]], self.z_coords[triangle_def[0]]]
        point2 = [self.x_coords[triangle_def[1]], self.y_coords[triangle_def[1]], self.z_coords[triangle_def[1]]]
        point3 = [self.x_coords[triangle_def[2]], self.y_coords[triangle_def[2]], self.z_coords[triangle_def[2]]]
        triangle = IfcTriangle(point1, point2, point3)

        return triangle 
    
    def get_temp_triangle(self, id):
        """Returns the copy IfcTriangle defined by the vertices and face on the given index needs index from the copy lists"""
        triangle_def = self.faces[id]
        point1 = [self.x_coords[triangle_def[0]], self.y_coords[triangle_def[0]], self.z_coords[triangle_def[0]]]
        point2 = [self.x_coords[triangle_def[1]], self.y_coords[triangle_def[1]], self.z_coords[triangle_def[1]]]
        point3 = [self.x_coords[triangle_def[2]], self.y_coords[triangle_def[2]], self.z_coords[triangle_def[2]]]
        triangle = IfcTriangle(point1, point2, point3)

        return triangle
    
    def set_edges(self, edges):
        """Updates the edges attribute in the element to the new edges"""
        vertices = self.verticies
        edges_vertices = [[vertices[edge[0]], vertices[edge[1]]] for edge in edges]
        self.edges = edges_vertices
    
    def set_verticies(self, verticies):
        """Updates all the data linked with the vertices when you have a new vertices list"""
        self.verticies = verticies
        self.x_coords = np.transpose(verticies)[0]
        self.y_coords = np.transpose(verticies)[1]
        self.z_coords = np.transpose(verticies)[2]
        self.x_min = min(self.x_coords)
        self.x_max = max(self.x_coords)
        self.y_min = min(self.y_coords)
        self.y_max = max(self.y_coords)
        self.z_min = min(self.z_coords)
        self.z_max = max(self.z_coords)

    def set_copy_verticies(self, verticies):
        """Updates all the data linked with the copy vertices when you have a new copy vertices list"""
        self.verticies_copy = verticies
        self.x_coords_copy = np.transpose(verticies)[0]
        self.y_coords_copy = np.transpose(verticies)[1]
        self.z_coords_copy = np.transpose(verticies)[2]
        self.x_min_copy = min(self.x_coords_copy)
        self.x_max_copy = max(self.x_coords_copy)
        self.y_min_copy = min(self.y_coords_copy)
        self.y_max_copy = max(self.y_coords_copy)
        self.z_min_copy = min(self.z_coords_copy)
        self.z_max_copy = max(self.z_coords_copy)

    def update_verts_and_faces(self, collisions=False):
        """Updates the orginal vertices and face list to be similar to the vertices_copy and faces_copy"""
        self.set_verticies(self.verticies_copy)
        self.faces = self.faces_copy
        if collisions:
            self.collisions = self.collisions_copy

    def new_verts_and_faces(self, triangles_obj, debug=False):
        """Adds the new vertices and faces to the end of the existing lists based on the triangles in the list 
        
        triangles_obj a list of IfcTriangles
        """
        input_is_list = isinstance(self.verticies_copy, list)

        if input_is_list:
            temp_verts = np.array(self.verticies_copy)
            temp_faces = np.array(self.faces_copy)
        else:
            temp_verts = self.verticies_copy
            temp_faces = self.faces_copy

        triangles = []
        for i in range(len(triangles_obj)):
            triangle = [triangles_obj[i].point1, triangles_obj[i].point2, triangles_obj[i].point3]
            triangles.append(triangle)

        triangles = np.array(triangles)
        new_points = np.unique(triangles.reshape(-1, 3), axis=0)
        
        temp_verts = np.vstack((temp_verts, new_points))
        
        new_faces = []
        for tri in triangles:
            face = [np.where(np.all(temp_verts == point, axis=1))[0][0] for point in tri]
            new_faces.append(face)

        if debug:
            print("len old faces:", len(temp_faces))
            print("len new faces:", len(temp_faces) + len(new_faces))

        temp_faces = np.vstack((temp_faces, new_faces))

        temp_collisions = [False] * len(new_faces)
        temp_collisions = self.collisions_copy + temp_collisions

        temp_verts = temp_verts.tolist()
        temp_faces = temp_faces.tolist()

        self.set_copy_verticies(temp_verts)
        self.faces_copy = temp_faces
        self.collisions_copy = temp_collisions

    def update_formwork_area(self):
        """Updates the formwork amount for the element based on the formwork list"""
        self.triangle_areas = [0] * len(self.faces)
        for i in range(len(self.faces)):
            triangle_area = self.get_triangle(i).area()
            self.triangle_areas[i] = triangle_area
            if self.formwork[i]:
                self.formwork_area = self.formwork_area + triangle_area

    def check_need_formwork(self, info):
        """Returns True if a element need formwork based on the production method and else False

        info the pset with the production method attribute
        """
        try: 
            prod_metode = info.get("01 Merknader").get("K11 Produksjonsmetode") # Set to work for the kvithammar aasen project
        except:
            a = "error"
        try: 
            prod_metode = info.get("E6KAA").get("K11 Produksjonsmetode") # Set to work for the kvithammar aasen project
        except:
            a = "error"
        if prod_metode == "Plasstøpt":
            return True
        return False
    
    def triangle_inside_element_in_list(self, triangle, elements, debug=False, margin=0.0001):
        """Returns True and the element if the triangle is inside any of the elements in the elements list else returns False and None

        triangle: IfcTriangle
        elements: list with IfcElements
        debug ooption enables printing
        margin adjusts the margin of a triangle to be classified as inside
        """
        for element in elements:
            jump = False
            if element is self:
                jump = True
            if triangle.inside_element(element, debug=debug, margin=margin) and not jump:
                return True, element
        return False, None
    
    def simplify_trimesh_element(self, target_face_count=None, debug=False):
        """Updates the vertices, faces and edges based the simplification by the trimesh library

        target_face_count: number of triangles in the finished element
        debug ooption enables printing
        """
        mesh = trimesh.Trimesh(vertices=self.verticies, faces=self.faces)

        if target_face_count is None:
            target_face_count = len(mesh.faces) // 2
            
        simplified_mesh = mesh.simplify_quadric_decimation(target_face_count)

        simplified_vertices = simplified_mesh.vertices.tolist()
        simplified_faces = simplified_mesh.faces.tolist()
        simplified_edges = simplified_mesh.edges_unique.tolist()

        if debug:
            print("len(faces)", len(self.faces))
            print("len(simplified_faces)", len(simplified_faces))
            print("len(edges)", len(self.edges))
            print("len(simplified_edges)", len(simplified_edges))

        self.set_verticies(simplified_vertices)
        self.set_copy_verticies(self.verticies)
        self.faces = simplified_faces
        self.faces_copy = self.faces
        self.set_edges(simplified_edges)
    
    def need_collision_test(self, element, debug=False):
        """Checks if the element needs collision_test with the element

        element: IfcElement
        debug option enables printing
        """
        if (not self.need_formwork and not element.need_formwork):
            return False

        if (self.x_min > element.x_max or self.x_max < element.x_min) or (self.y_min > element.y_max or self.y_max < element.y_min) or (self.z_min > element.z_max or self.z_max < element.z_min):
            if debug:
                print("collision test not needed")
            return False
        
        threshold = 0.01
        safe = 0
        if (self.x_min >= element.x_max - threshold or self.x_max <= element.x_min + threshold):
            safe += 1
        if (self.y_min >= element.y_max - threshold or self.y_max <= element.y_min + threshold):
            safe += 1
        if (self.z_min >= element.z_max - threshold or self.z_max <= element.z_min + threshold):
            safe += 1

        if safe > 1:
            if debug:
                print("collision test not needed")
            return False

        return True    

    def collisions_two_elements(self, element2, compute_new_triangles=False, compute_new_triangles_2d=False, debug=False):
        """calculates the collisions between two elements and updates the triangles in each element based on the collisions

        element2: IfcElement
        compute_new_triangles: boolean (are we computing new triangles based on 3d collisons?)
        compute_new_triangles_2d: boolean (are we computing new triangles based on 2d collisons?)
        debug option enables printing
        """
        i = 0
        while i < len(self.faces_copy):
            collision = False
            triangle1 = self.get_temp_triangle(i)
            j = 0
            while j < len(element2.faces):
                temp_collision = False
                triangle2 = element2.get_triangle(j)
                if debug:
                    print("(self, element)", self.guid, element2.guid)
                    print("(i, j)", (i, j))
                
                ### Not the best solution, but for now
                compute_new_triangles_2d_1 = compute_new_triangles_2d
                if compute_new_triangles_2d and (i >= len(self.faces_changed) or j >= len(element2.faces_changed)):
                    compute_new_triangles_2d_1 = False
                elif (i < len(self.faces_changed) and j < len(element2.faces_changed)):
                    if compute_new_triangles_2d and (self.faces_changed[i] or element2.faces_changed[j]):
                        compute_new_triangles_2d_1 = False
                ### needs development

                temp_collision, new_triangles_1, new_triangles_2, changed_2d, collision_2d = triangle1.collision_test(triangle2, compute_new_triangles=compute_new_triangles, compute_new_triangles_2d=compute_new_triangles_2d, debug=debug)

                ### needs to make new function for computing new triangles in 2d!!
                if collision_2d and compute_new_triangles_2d_1:
                    if debug:
                        print("start compute new_triangles in 2d")
                    new_triangles_1 = triangle1.compute_new_triangles_2d_v4(element2, debug=debug)
                    if new_triangles_1 == None:
                        new_triangles_2 = triangle2.compute_new_triangles_2d_v4(self, debug=debug)
                ### This does not work

                if temp_collision:
                    collision = True
                if new_triangles_1 != None:
                    self.new_verts_and_faces(new_triangles_1, debug=debug)
                    if debug:
                        print("New_triangles != None")
                        print("triangle1_def", self.faces_copy[i])
                    self.faces_copy.remove(self.faces_copy[i])
                    self.collisions_copy.pop(i)
                    i = i - 1
                    if debug:
                        print("3dcollision, triangles updated")
                    if changed_2d:
                        self.faces_changed[i] = True
                    break
                j = j + 1
            if new_triangles_1 is None and new_triangles_2 is None and collision:
                self.collisions_copy[i] = True
            i = i + 1



#########################################################################################################################
class IfcTriangle():
    def __init__(self, point1, point2, point3):
        """The IfcTriangle class represents a triangle and takes care of all calculations and classification on triangles

        point1: [x, y, z]
        point2: [x, y, z]
        point3: [x, y, z]
        """
        self.point1 = self.make_list_float(point1)
        self.point2 = self.make_list_float(point2)
        self.point3 = self.make_list_float(point3)
        self.x_coords = [self.point1[0], self.point2[0], self.point3[0]]
        self.y_coords = [self.point1[1], self.point2[1], self.point3[1]]
        self.z_coords = [self.point1[2], self.point2[2], self.point3[2]]
        self.x_min = min(self.x_coords)
        self.x_max = max(self.x_coords)
        self.y_min = min(self.y_coords)
        self.y_max = max(self.y_coords)
        self.z_min = min(self.z_coords)
        self.z_max = max(self.z_coords)
        self.normal_vector = self.get_normal_vector() # Should not compute this every time, not always needed

    def get_point_near_middle_of_triangle(self):
        """Returns point in the middle of the triangle"""
        middle_point = (np.array(self.point1) + np.array(self.point2) + np.array(self.point3)) / 3
        return middle_point

    def get_point_inside_triangle(self):
        """Returns random point inside the triangle"""
        r, s = np.random.rand(2)
        if r + s >= 1:
            r = 1 - r
            s = 1 - s
        new_point = np.array(self.point1) + r * (np.array(self.point2) - np.array(self.point1)) + s * (np.array(self.point3) - np.array(self.point1))
        return new_point.tolist()

    def get_min_max_2d_triangle(self, triangle):
        """Returns the x_min, x_max, y_min, y_max from triangle defined in 2d
        
        triangle: [[x, y], [x1, y1], [x2, y2]]
        """
        x_coords = [triangle[0][0], triangle[1][0], triangle[2][0]]
        y_coords = [triangle[0][1], triangle[1][1], triangle[2][1]]
        x_min = min(x_coords)
        x_max = max(x_coords)
        y_min = min(y_coords)
        y_max = max(y_coords)
        return x_min, x_max, y_min, y_max
    
    def get_normal_vector(self):
        """Returns the normal vector to the triangle"""
        v1 = np.array(self.point2) - np.array(self.point1)
        v2 = np.array(self.point3) - np.array(self.point1)
        normal_vector = np.cross(v1, v2)
        return normal_vector / np.linalg.norm(normal_vector)

    def area(self):
        """Returns the area of the triangle"""
        return 0.5 * np.linalg.norm(np.cross(np.array(self.point2) - np.array(self.point1), np.array(self.point3) - np.array(self.point1)))

    def print(self):
        """Prints the three points of the triangle"""
        print("Triangle points:", self.point1, self.point2, self.point3)

    def is_close(self, value1, value2, threshold):
        """Returns True if the values are closer than the threshold else False"""
        if abs(value1-value2) < threshold:
            return True
        return False
    
    def point_is_close(self, point1, point2, tolerance):
        """Returns True if the 2d points are closer than the threshold else False"""
        distance = np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
        return distance < tolerance
    
    def are_edges_collinear_2d(self, triangle1, triangle2, threshold=0.001):
        """Returns the two edges that are coliniar, else None, None"""
        def vector(p1, p2):
            return np.array(p2) - np.array(p1)
        
        def is_collinear(v1, v2, threshold=0.001):
            cross_product = np.cross(v1, v2)
            norm = np.linalg.norm(cross_product)
            return abs(norm) < threshold

        # Create edges from triangles
        edges1 = [(triangle1[i], triangle1[(i+1) % 3]) for i in range(3)]
        edges2 = [(triangle2[i], triangle2[(i+1) % 3]) for i in range(3)]
        
        # Check each pair of edges
        for edge1 in edges1:
            for edge2 in edges2:
                v1 = vector(*edge1)
                v2 = vector(*edge2)
                if is_collinear(v1, v2, threshold=threshold):
                    return edge1, edge2  # Found collinear edges
        return None, None
    
    def orientation(self, a, b, c):
        """Returns 1 if a are on one side of b-c, 0 if they are colinear, else 2"""
        val = (b[1] - a[1]) * (c[0] - b[0]) - (b[0] - a[0]) * (c[1] - b[1])
        if abs(val) < 0.001: return 0 
        return 1 if val > 0 else 2 
    
    def make_list_float(self, numbers):
        """Taks in numbers and returns the list as float"""
        numbers = list(numbers)
        for i in range(len(numbers)):
            numbers[i] = float(numbers[i])
        return numbers
    
    def rotate_vector(self, vector):
        """Rotates the vector"""
        return [vector[2], vector[0], vector[1]]
    
    def change_order(self, vector):
        """Changes order on the vector"""
        return [vector[0], vector[2], vector[1]]
    
    def rotate(self):
        """Rotates the points that defines the triangle"""
        temp = self.point1
        self.point1 = self.point3
        self.point3 = self.point2
        self.point2 = temp

    def ensure_counterclockwise_2d(self, triangle):
        """Returns the 2d triangle defined counterclockwise"""
        area = 0.5 * (triangle[0][0]*(triangle[1][1]-triangle[2][1]) + triangle[1][0]*(triangle[2][1]-triangle[0][1]) + triangle[2][0]*(triangle[0][1]-triangle[1][1]))
        if area < 0:
            triangle[1], triangle[2] = triangle[2], triangle[1]

        return triangle

    def fix_permutation(self, det_values_1, det_values_2, triangle1, triangle2, depth=0, debug=False):
        """Returns the signs for the triangles and the triangles with fixed permutation. Mening the vertices with oposite sign first
        
        det_values_1: The determinant sign for triangle 1 compared with triangle2
        det_values_2: The determinant sign for triangle 2 compared with triangle1
        triangle1: IfcTriangle
        triangle2: IfcTriangle
        """
        if depth == 6:
            if debug:
                print("depth_reached")
            return det_values_1, det_values_2, triangle1, triangle2

        if(det_values_1[0] >= 0 and det_values_1[1] <= 0 and det_values_1[2] <= 0) or (det_values_1[0] <= 0 and det_values_1[1] >= 0 and det_values_1[2] >= 0):
            if(det_values_2[0] >= 0 and det_values_2[1] <= 0 and det_values_2[2] <= 0) or (det_values_2[0] <= 0 and det_values_2[1] >= 0 and det_values_2[2] >= 0):
                return det_values_1, det_values_2, triangle1, triangle2
            
            det_values_2 = self.rotate_vector(det_values_2)
            if debug:
                print("triangle2 before rotate")
                triangle2.print()
            triangle2.rotate()
            if debug:
                print("triangle2 after rotate")
                triangle2.print()
            return self.fix_permutation(det_values_1, det_values_2, triangle1, triangle2, depth + 1)

        det_values_1 = self.rotate_vector(det_values_1)
        if debug:
            print("triangle1 before rotate")
            triangle1.print()
        triangle1.rotate()
        if debug:
            print("triangle1 after rotate")
            triangle1.print()
        return self.fix_permutation(det_values_1, det_values_2, triangle1, triangle2, depth + 1)
    

    def det_sign_2d_v2(self, triangle1, triangle2):
        """Returns the sign of the determinant of the vertices in triangle 1 compared with triangle 2 (2D)
        
        triangle1: [[x1, y1], ...]
        triangle2: [[x1, y1], ...]
        """
        det_values_finale = []
        for i in range(3):
            det_values = []
            det = self.point_orientation_threshold(triangle2[0], triangle2[1], triangle1[i])
            det_values.append(det)
            det = self.point_orientation_threshold(triangle2[1], triangle2[2], triangle1[i])
            det_values.append(det)
            det = self.point_orientation_threshold(triangle2[2], triangle2[0], triangle1[i])
            det_values.append(det)
            det_values_finale.append(det_values)
        return det_values_finale
    
    def need_collision_test(self, triangle, debug=False):
        """Returns True if the triangles are within each others bounding box, else False"""
        if (self.x_min > triangle.x_max or self.x_max < triangle.x_min) or (self.y_min > triangle.y_max or self.y_max < triangle.y_min) or (self.z_min > triangle.z_max or self.z_max < triangle.z_min):
                if debug:
                    print("collisions test not needed")
                return False   
        
        threshold = 0.01
        safe = 0
        if (self.x_min >= triangle.x_max - threshold or self.x_max <= triangle.x_min + threshold):
            safe += 1
        if (self.y_min >= triangle.y_max - threshold or self.y_max <= triangle.y_min + threshold):
            safe += 1
        if (self.z_min >= triangle.z_max - threshold or self.z_max <= triangle.z_min + threshold):
            safe += 1

        if safe > 1:
            if debug:
                print("collision test not needed")
            return False

        return True
    
    def is_normal_vector_pointing_out(self, element, debug=False):
        """Returns True if the triangles normal_vector is pointing out of the element it belongs to, else False

        element: IfcElement
        debug option enables printing
        """       
        def ray_intersects_triangle(ray_origin, ray_vector, triangle):
            
            ray_origin = np.array(ray_origin)
            ray_vector = np.array(ray_vector)
            vertex0, vertex1, vertex2 = np.array(triangle.point1), np.array(triangle.point2), np.array(triangle.point3)

            edge1 = vertex1 - vertex0
            edge2 = vertex2 - vertex0

            pvec = np.cross(ray_vector, edge2)

            det = edge1.dot(pvec)
            if abs(det) < 0.1:
                return False
            
            inv_det = 1.0 / det
            tvec = ray_origin - vertex0
            u = tvec.dot(pvec) * inv_det
            if u < 0 or u > 1:
                return False

            qvec = np.cross(tvec, edge1)
            v = ray_vector.dot(qvec) * inv_det
            if v < 0 or u + v > 1:
                return False

            t = edge2.dot(qvec) * inv_det

            return t >= 0.01
        
        ray_origin = self.get_point_near_middle_of_triangle()
        ray_vector = self.normal_vector

        if debug:
            print("ray_origin:", ray_origin)
            print("ray_vector:", ray_vector)
        
        intersections = 0
        for i in range(len(element.faces)):
            triangle = element.get_triangle(i)
            intersects = ray_intersects_triangle(ray_origin, ray_vector, triangle)
            if intersects:
                if debug:
                    print("Normal vector intersects triangle")
                    print("Triangle:", triangle)
                intersections += 1
        
        if debug:
            print("Number of intersections", intersections)

        if intersections % 2 == 1:
            return False

        return True 
    
    def inside_element(self, element, margin=0.0001, debug=False):
        """Returns True if the triangle is inside the bounding box of the element, else False

        element: IfcElement
        margin: the margin for if the point is inside
        debug option enables printing
        """    
        if (self.x_max <= element.x_max + margin and self.x_min >= element.x_min - margin) and (self.y_max <= element.y_max + margin and self.y_min >= element.y_min - margin) and (self.z_max <= element.z_max + margin and self.z_min >= element.z_min - margin):
            if debug:
                print("triangle is inside element")
            return True
        return False
    
    def inside_element_v5(self, element, debug=False):
        """Returns True if the triangle is inside the element, else False
        using ray-casting

        element: IfcElement
        debug option enables printing
        """   
        if self.is_point_inside_element(self.point1, element, debug=debug) or self.is_point_inside_element(self.point2, element, debug=debug) or self.is_point_inside_element(self.point2, element, debug=debug) or self.is_point_inside_element(self.get_point_near_middle_of_triangle(), element, debug=debug):
            return True
        return False

    ##########  Should change it use the old version and using the original triangles from the element instead
    ### NB ###
    ##########
    def is_point_inside_element(self, point, element, debug=False, new_version=False):
        """Returns True if the point is inside the element, else False

        element: IfcElement
        point: [x, y, z]
        debug option enables printing
        """    
        def ray_intersects_triangle(ray_origin, direction, triangle, y_direction=0, z_direction=0, new_version=False):
            if (not new_version and ((ray_origin[1] > triangle.y_max or ray_origin[1] < triangle.y_min) or (ray_origin[2] > triangle.z_max or ray_origin[2] < triangle.z_min))):
                return False
            
            ray_origin = np.array(ray_origin)
            ray_vector = np.array([direction, y_direction, z_direction])
            vertex0, vertex1, vertex2 = np.array(triangle.point1), np.array(triangle.point2), np.array(triangle.point3)

            edge1 = vertex1 - vertex0
            edge2 = vertex2 - vertex0

            pvec = np.cross(ray_vector, edge2)

            det = edge1.dot(pvec)
            if abs(det) < 0.001:
                return False
            
            inv_det = 1.0 / det
            tvec = ray_origin - vertex0
            u = tvec.dot(pvec) * inv_det
            if u < 0 or u > 1:
                return False

            qvec = np.cross(tvec, edge1)
            v = ray_vector.dot(qvec) * inv_det
            if v < 0 or u + v > 1:
                return False

            t = edge2.dot(qvec) * inv_det

            return t >= 0.01
        
        if new_version:
            intersection_list = []
            for i in range(2):
                intersections = 0
                direction = 1 - 2*i
                y_direction = 0 - 1*i
                z_direction = 0
                if debug:
                    print("x, y, z", direction, y_direction, z_direction)
                for i in range(len(element.faces_org)):
                    triangle = element.get_org_triangle(i)
                    intersects = ray_intersects_triangle(point, direction, triangle, y_direction=y_direction, z_direction=z_direction, new_version=True)
                    if intersects:
                        if debug:
                            print("Point intersects triangle")
                            print("triangle:", triangle.print())
                            print("point", point)
                        intersections += 1
                if debug:
                    print("Number of intersections:", intersections)
                intersection_list.append(intersections % 2 == 1)

            if debug:
                print("intersection list:", intersection_list)
            return sum(1 for ct in intersection_list if ct) > 1
        
        intersection_list = []
        for i in range(2):
            intersections = 0
            direction = 1 - 2*i
            for i in range(len(element.faces)):
                triangle = element.get_triangle(i)
                intersects = ray_intersects_triangle(point, direction, triangle)
                if intersects:
                    if debug:
                        print("Point intersects triangle")
                    intersections += 1
            if debug:
                print("Number of intersections:", intersections)
            intersection_list.append(intersections % 2 == 1)

        return intersection_list[0] and intersection_list[1]
    
    def bigger(self, triangle1_2d, triangle2_2d):
        """Returns True if the triangle1_2d is bigger than triangle2_2d

        triangle1_2d: [[x, y], [....], [....]]
        triangle2_2d: [[x, y], [....], [....]]
        """    
        x_min, x_max, y_min, y_max = self.get_min_max_2d_triangle(triangle1_2d)
        x_min_2, x_max_2, y_min_2, y_max_2 = self.get_min_max_2d_triangle(triangle2_2d)

        if (x_max >= x_max_2 and x_min <= x_min_2) and (y_max >= y_max_2 and y_min <= y_min_2):
                return True
        if (x_max_2 >= x_max and x_min_2 <= x_min) and (y_max_2 >= y_max and y_min_2 <= y_min):
                return True   
        
        return False
    
    def triangle_similarity_v2(self, triangle1, triangle2, threshold=0.001, debug=False):
        """Returns True if the triangle1 is similar to triangle2, mening they have colinear and overlapping edges and that the last vertices ly on the same side.

        triangle1_2d: [[x, y], [....], [....]]
        triangle2_2d: [[x, y], [....], [....]]
        threshold: float
        debug enables printing
        """  
        edge1, edge2 = self.are_edges_collinear_2d(triangle1, triangle2, threshold=threshold)

        if edge1 is None or edge2 is None:
            if debug:
                print("No edges are coliniar")
            return False
        if debug:
            print("edge_1:", edge1)
            print("edge_2:", edge2)

        edge_1_bounds = [[min(edge1[0][i], edge1[1][i]), max(edge1[0][i], edge1[1][i])] for i in range(2)]
        edge_2_bounds = [[min(edge2[0][i], edge2[1][i]), max(edge2[0][i], edge2[1][i])] for i in range(2)]

        if debug:
            print("edge_1_bounds:", edge_1_bounds)
            print("edge_2_bounds:", edge_2_bounds)
        
        overlap = all(edge_1_bounds[dim][0] < edge_2_bounds[dim][1] + threshold and edge_1_bounds[dim][1] > edge_2_bounds[dim][0] - threshold for dim in range(2))

        if not overlap:
            if debug:
                print("The colinear edges are not overlapping")
            return False

        if debug:
            print("The colinear edges are overlapping")
        
        o1 = self.point_orientation_threshold(edge1[0], edge1[1], triangle1[0])
        o2 = self.point_orientation_threshold(edge1[0], edge1[1], triangle1[1])
        o3 = self.point_orientation_threshold(edge1[0], edge1[1], triangle1[2])

        tri_1_or = o1 + o2 + o3

        o4 = self.point_orientation_threshold(edge1[0], edge1[1], triangle2[0])
        o5 = self.point_orientation_threshold(edge1[0], edge1[1], triangle2[1])
        o6 = self.point_orientation_threshold(edge1[0], edge1[1], triangle2[2])

        tri_2_or = o4 + o5 + o6

        if debug:
            print("o1, o2, o3:", o1, o2, o3)
            print("o4, o5, o6:", o4, o5, o6)

        if (tri_1_or > 0 and tri_2_or > 0) or (tri_1_or < 0 and tri_2_or < 0):
            if debug:
                print("The thrid point lies on the same side")
            return True
        
        return False
    
    def point_orientation_threshold(self, p1, p2, p, threshold=0.001):
        """
        Determines if point p is to the left, right, or on the line defined by segment p1p2, with a threshold for points
        close to the line.
        :param p1: First point of the line segment (x1, y1).
        :param p2: Second point of the line segment (x2, y2).
        :param p: The point to check (x, y).
        :param threshold: Threshold for considering points close to the line.
        :return: 1 if point p is to the left, -1 if it's to the right, 0 if it's within the threshold of the line.
        """
        vector_p1p2 = (p2[0] - p1[0], p2[1] - p1[1])
        vector_p1p = (p[0] - p1[0], p[1] - p1[1])

        determinant = vector_p1p2[0] * vector_p1p[1] - vector_p1p2[1] * vector_p1p[0]

        if determinant > threshold:
            return 1  
        elif determinant < -threshold:
            return -1  
        else:
            return 0 
        
    def check_line_intersection(self, p1, p2, q1, q2, debug=False):
        """Returns True if the two lines intersects, else False"""
        def on_segment(a, b, c):
            return min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])
        
        def distance_point_to_line(px, py, x1, y1, x2, y2):
            norm = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
            distance = abs((x2 - x1)*(y1 - py) - (x1 - px)*(y2 - y1)) / norm
            return distance

        def is_close(p, q1, q2, tolerance):
            d = distance_point_to_line(p[0], p[1], q1[0], q1[1], q2[0], q2[1])
            return d < tolerance

        tolerance = 0.00000001

        o1 = self.orientation(p1, p2, q1)
        o2 = self.orientation(p1, p2, q2)
        o3 = self.orientation(q1, q2, p1)
        o4 = self.orientation(q1, q2, p2)

        if o1 == 0 or o2 == 0 or o3 == 0 or o4 == 0:
            return False

        if o1 != o2 and o3 != o4:
            if (((is_close(q1, p1, p2, tolerance)) and on_segment(p1, p2, q1)) or
                ((is_close(q2, p1, p2, tolerance)) and on_segment(p1, p2, q2)) or
                ((is_close(p1, q1, q2, tolerance)) and on_segment(q1, q2, p1)) or
                ((is_close(p2, q1, q2, tolerance)) and on_segment(q1, q2, p2))):
                return False
            if debug:
                print("Genral case collision")
            return True

        return False
    
    def check_triangle_intersection(self, triangle1, triangle2):
        """Returns True if any of the lines in the triangles intersects, else False
        
        triangle1: [[x, y], ...]
        triangle2: [[x, y], ...]
        """
        for i in range(3):
            for j in range(3):
                if self.check_line_intersection(triangle1[i], triangle1[(i+1)%3], triangle2[j], triangle2[(j+1)%3]):
                    return True
        return False