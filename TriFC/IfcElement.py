import numpy as np
import trimesh
from TriFC.IfcTriangle import IfcTriangle

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
        if info is not None and need_formwork is None:
            self.need_formwork = self.check_need_formwork(info)
        elif need_formwork is None:
            self.need_formwork = True
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
        self.half_formwork = [False] * len(faces)
        self.formwork_info = None
        self.formwork_area = 0
        self.triangle_areas = []
        self.normal_vectors = [None] * len(faces)

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
        triangle_def = self.faces_org[id]
        point1 = [self.x_coords_org[triangle_def[0]], self.y_coords_org[triangle_def[0]], self.z_coords_org[triangle_def[0]]]
        point2 = [self.x_coords_org[triangle_def[1]], self.y_coords_org[triangle_def[1]], self.z_coords_org[triangle_def[1]]]
        point3 = [self.x_coords_org[triangle_def[2]], self.y_coords_org[triangle_def[2]], self.z_coords_org[triangle_def[2]]]
        triangle = IfcTriangle(point1, point2, point3)

        return triangle 
    
    def get_temp_triangle(self, id):
        """Returns the copy IfcTriangle defined by the vertices and face on the given index needs index from the copy lists"""
        triangle_def = self.faces_copy[id]
        point1 = [self.x_coords_copy[triangle_def[0]], self.y_coords_copy[triangle_def[0]], self.z_coords_copy[triangle_def[0]]]
        point2 = [self.x_coords_copy[triangle_def[1]], self.y_coords_copy[triangle_def[1]], self.z_coords_copy[triangle_def[1]]]
        point3 = [self.x_coords_copy[triangle_def[2]], self.y_coords_copy[triangle_def[2]], self.z_coords_copy[triangle_def[2]]]
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
        self.half_formwork = [False]*len(self.faces_copy)
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

    def update_formwork_area(self, wall_formwork=False):
        """Updates the formwork amount for the element based on the formwork list
        
        wall_formwork: boolean 
        """
        self.triangle_areas = [0] * len(self.faces)
        for i in range(len(self.faces)):
            triangle_area = self.get_triangle(i).area()
            self.triangle_areas[i] = triangle_area
            if self.formwork[i]:
                if wall_formwork and self.half_formwork[i]:
                    self.formwork_area = self.formwork_area + (triangle_area/2)
                else:
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
                
                # ### Not the best solution, but for now
                # compute_new_triangles_2d_1 = compute_new_triangles_2d
                # if compute_new_triangles_2d and (i >= len(self.faces_changed) or j >= len(element2.faces_changed)):
                #     compute_new_triangles_2d_1 = False
                # elif (i < len(self.faces_changed) and j < len(element2.faces_changed)):
                #     if compute_new_triangles_2d and (self.faces_changed[i] or element2.faces_changed[j]):
                #         compute_new_triangles_2d_1 = False
                # ### needs development

                temp_collision, new_triangles_1, new_triangles_2, changed_2d, collision_2d = triangle1.collision_test(triangle2, compute_new_triangles=compute_new_triangles, compute_new_triangles_2d=compute_new_triangles_2d, debug=debug)

                # ### needs to make new function for computing new triangles in 2d!!
                # if collision_2d and compute_new_triangles_2d_1:
                #     if debug:
                #         print("start compute new_triangles in 2d")
                #     new_triangles_1 = triangle1.compute_new_triangles_2d_v4(element2, debug=debug)
                #     if new_triangles_1 == None:
                #         new_triangles_2 = triangle2.compute_new_triangles_2d_v4(self, debug=debug)
                # ### This does not work

                if temp_collision:
                    collision = True
                if new_triangles_1 != None:
                    self.new_verts_and_faces(new_triangles_1, debug=debug)
                    if debug:
                        print("New_triangles != None")
                        print("triangle1_def", self.faces_copy[i])
                        print("new_triangles", new_triangles_1[0].point1, new_triangles_1[1].point1)
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