import numpy as np
import math
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
    
    def ensure_counterclockwise(self):
        """Makes sure that the triangle is defined counterclockwise"""
        edge1 = np.array(self.point2) - np.array(self.point1)
        edge2 = np.array(self.point3) - np.array(self.point1)
        normal = np.cross(edge1, edge2)
        dot_product = np.dot(normal, np.array(self.point1))
        if dot_product < 0:
            temp = self.point3
            self.point3 = self.point2
            self.point3 = temp

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
    
    def project_to_plane(self, points, normal, origin=None, debug=False):
        """Projects the points in the list to a 2d-plane defines by origin and a normal vector, 
        returns origin, axis1, axis2, projected"""
        if origin is None:
            origin = np.array(points[0])

        axis1 = np.cross(normal, [1, 0, 0]) if normal[0] == 0 else np.cross(normal, [0, 1, 0])
        axis1 = axis1 / np.linalg.norm(axis1)
        axis2 = np.cross(normal, axis1)

        projected = [[np.dot(p - origin, axis1), np.dot(p - origin, axis2)] for p in points]

        return origin, axis1, axis2, projected

    def project_from_plan(self, triangles, origin, axis1, axis2, debug=False):
        """Projects triangles back to 3d given their projection origin and axis,
        returns the triangles as IfcTriangles"""
        new_triangles = []
        for triangle in triangles:
            points_3d = []
            for point in triangle:
                point_3d = origin + point[0] * axis1 + point[1] * axis2
                points_3d.append(point_3d)
            new_triangles.append(IfcTriangle(points_3d[0], points_3d[1], points_3d[2]))
        
        return new_triangles
    
    def det_sign_check(self, triangle1, triangle2):
        """Returns the sign of the determinant of the vertices in triangle 1 compared with triangle 2 (3D)
        
        triangle1: IfcTriangle
        triangle2: IfcTriangle
        """
        det_values = []
        triangle_points = [triangle1.point1, triangle1.point2, triangle1.point3]
        for i in range(len(triangle_points)):
            matrix = np.array([triangle2.point1 + [1], triangle2.point2 + [1], triangle2.point3 + [1], triangle_points[i] + [1]])
            det = np.linalg.det(matrix)
            det_values.append(det)
        
        return det_values

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
    
    def is_normal_vector_pointing_out(self, element, normal_vector, debug=False):
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
        ray_vector = normal_vector

        if debug:
            print("ray_origin:", ray_origin)
            print("ray_vector:", ray_vector)
        
        intersections = 0
        for i in range(len(element.faces)):
            triangle = element.get_triangle(i)
            intersects = ray_intersects_triangle(ray_origin, ray_vector, triangle)
            if intersects:
                if debug:
                    print("Normal vector intersects triangle in element")
                intersections += 1
        
        if debug:
            print("Number of intersections", intersections)

        if intersections % 2 == 1:
            return False

        return True   
    
    def is_point_coplanar(self, point, debug=False): 
        """Returns True if the point is coplanar with the triangle"""
        a1 = self.point2[0] - self.point1[0]
        b1 = self.point2[1] - self.point1[1]
        c1 = self.point2[2] - self.point1[2]
        a2 = self.point3[0] - self.point1[0]
        b2 = self.point3[1] - self.point1[1]
        c2 = self.point3[2] - self.point1[2]
        a = b1 * c2 - b2 * c1
        b = a2 * c1 - a1 * c2
        c = a1 * b2 - b1 * a2
        d = (- a * self.point1[0] - b * self.point1[1] - c * self.point1[2])
        value = abs(a * point[0] + b * point[1] + c * point[2] + d)
        if(value < 0.001):
            if debug:
                print("point is coplanar")
            return True
        else:
            return False

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
    
    def inside_element_v6(self, element, debug=False):
        """Returns True if the triangle is inside the element, else False
        using ray-casting

        element: IfcElement
        debug option enables printing
        """   
        if self.is_point_inside_element(self.get_point_near_middle_of_triangle(), element, debug=debug, new_version=False):
            return True
        return False

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
            for i in range(len(element.faces_org)):
                triangle = element.get_org_triangle(i)
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

    def add_v3v3(self, v0, v1):
        return (
            v0[0] + v1[0],
            v0[1] + v1[1],
            v0[2] + v1[2],
        )

    def sub_v3v3(self, v0, v1):
        return (
            v0[0] - v1[0],
            v0[1] - v1[1],
            v0[2] - v1[2],
        )

    def dot_v3v3(self, v0, v1):
        return (
            (v0[0] * v1[0]) +
            (v0[1] * v1[1]) +
            (v0[2] * v1[2])
        )

    def len_squared_v3(self, v0):
        return self.dot_v3v3(v0, v0)

    def mul_v3_fl(self, v0, f):
        return (
            v0[0] * f,
            v0[1] * f,
            v0[2] * f,
        )
    
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
    
    def isect_line_plane_v3(self, p0, p1, p_co, p_no, epsilon=1e-6, debug=False):
        """Returns a point if the intersection between line and plan exists and else None
        
        p0: point on line
        p1: point on line
        p_c0: point on plane
        p_no: normal vector plane
        epsilon: threshold for parallel decision 
        debug: enables printing
        """
        u = self.sub_v3v3(p1, p0)
        dot = self.dot_v3v3(p_no, u)
        if debug:
            print("dot:", dot)
            print("epsilon: ", epsilon)

        if abs(dot) > epsilon:
            w = self.sub_v3v3(p0, p_co)
            fac = -self.dot_v3v3(p_no, w) / dot
            u = self.mul_v3_fl(u, fac)
            return self.add_v3v3(p0, u)

        return None
    
    def isect_triangle_edges_plane_v3(self, triangle, p_co, p_no, epsilon=1e-6, debug=False):
        """Returns a list with all the intersection points between the plane and the triangle edges, none for no intersection
        
        triangle: IfcTriangle
        p_c0: point on plane
        p_no: normal vector plane
        epsilon: threshold for parallel decision 
        debug: enables printing
        """
        # use permutation!! 
        point_pairs = [[triangle.point2, triangle.point1], [triangle.point3, triangle.point1], [triangle.point3, triangle.point2]]
        u1, u2 = self.sub_v3v3(triangle.point2, triangle.point1), self.sub_v3v3(triangle.point3, triangle.point1)
        edges = [u1, u2]
        intersections = []
        for i in range(len(edges)):
            u = edges[i]
            dot = self.dot_v3v3(p_no, u)
            parallel_check = self.dot_v3v3(p_no / np.linalg.norm(p_no), u / np.linalg.norm(u)) # test
            if debug:
                print("dot:", dot)
                print("epsilon: ", epsilon)
            if abs(parallel_check) > epsilon:
                w = self.sub_v3v3(point_pairs[i][0], p_co)
                fac = -self.dot_v3v3(p_no, w) / dot
                u = self.mul_v3_fl(u, fac)
                intersections.append(self.add_v3v3(point_pairs[i][0], u))
            else:
                intersections.append(None)
        return intersections
    
    def is_point_in_triangle(self, pt, v1, v2, v3, debug=False):
        """Returns True if the point is inside the triangle, else False
        
        pt: point
        v1: point
        v2: point
        v3: point
        debug: enables printing
        """
        if pt == None:
            return False
        
        pt, v1, v2, v3 = np.array(pt), np.array(v1), np.array(v2), np.array(v3)
        vec_v1pt = pt - v1
        vec_v2pt = pt - v2
        vec_v3pt = pt - v3

        tri_normal = np.cross(v2 - v1, v3 - v1)

        normal1 = np.cross(vec_v1pt, vec_v2pt)
        normal2 = np.cross(vec_v2pt, vec_v3pt)
        normal3 = np.cross(vec_v3pt, vec_v1pt)

        if debug:
            print(np.dot(normal1, tri_normal))
            print(np.dot(normal2, tri_normal))
            print(np.dot(normal3, tri_normal))

        inside = (np.dot(normal1, tri_normal) >= -0.00001 and
                np.dot(normal2, tri_normal) >= -0.00001 and
                np.dot(normal3, tri_normal) >= -0.00001)

        return inside

    def is_point_in_triangle_v2(self, pt, debug=False):
        """Returns True if the point is inside the triangle, else False
        
        pt: point
        debug: enables printing
        """
        v1 = self.point1
        v2 = self.point2
        v3 = self.point3

        if pt is None:
            if debug:
                print("The point is None")
            return False
        
        pt, v1, v2, v3 = np.array(pt), np.array(v1), np.array(v2), np.array(v3)
        vec_v1pt = pt - v1
        vec_v2pt = pt - v2
        vec_v3pt = pt - v3

        tri_normal = np.cross(v2 - v1, v3 - v1)

        normal1 = np.cross(vec_v1pt, vec_v2pt)
        normal2 = np.cross(vec_v2pt, vec_v3pt)
        normal3 = np.cross(vec_v3pt, vec_v1pt)

        if debug:
            print(np.dot(normal1, tri_normal))
            print(np.dot(normal2, tri_normal))
            print(np.dot(normal3, tri_normal))

        inside = (np.dot(normal1, tri_normal) >= -0.000001 and
                np.dot(normal2, tri_normal) >= -0.000001 and
                np.dot(normal3, tri_normal) >= -0.000001)

        return inside

    def order_points_on_line(self, points, debug=False):
        '''Returns a orderd list based on the vector p0-p1 from the list'''
        if debug:
            print("line_points in order points", points)
        points = np.array(points)
        direction_vector = points[1] - points[0]
        scalar_projections = np.dot(points - points[0], direction_vector) / np.dot(direction_vector, direction_vector)
        ordered_indices = np.argsort(scalar_projections)
        ordered_points = points[ordered_indices]

        return ordered_points.tolist()
    
    def order_points_along_vector(self, points, vector_start, vector_end):
        '''Returns a orderd list based on the vector created from the start, end from the points list'''
        direction_vector = np.array(vector_end) - np.array(vector_start)
        direction_vector /= np.linalg.norm(direction_vector) 
        
        distances = [np.dot(np.array(point) - np.array(vector_start), direction_vector) for point in points]
        
        sorted_indices = sorted(range(len(distances)), key=lambda k: distances[k])
        sorted_points = [points[i] for i in sorted_indices]
        
        return sorted_points
    
    def get_unique_points_from_list(self, list, list2=None):
        '''Returns a list with no duplicated 3d-points'''
        unique_list = []
        if list2 != None:
            list = list+list2

        for i in range(len(list)):
            unique = True
            for j in range(len(list)):
                if j == i:
                    unique = True
                else:
                    x = abs(list[i][0] - list[j][0])
                    y = abs(list[i][1] - list[j][1])
                    z = abs(list[i][2] - list[j][2])
                    if x < 0.001 and y < 0.001 and z < 0.001:
                        unique = False
            if unique:
                unique_list.append(list[i])

        return unique_list
    
    def get_unique_points_from_list_2d(self, list, list2=None):
        '''Returns a list with no duplicated 2d-points'''
        unique_list = []
        if list2 != None:
            list = list+list2

        for i in range(len(list)):
            unique = True
            for j in range(len(list)):
                if j == i:
                    unique = True
                else:
                    x = abs(list[i][0] - list[j][0])
                    y = abs(list[i][1] - list[j][1])
                    if x < 0.001 and y < 0.001:
                        unique = False
            if unique:
                unique_list.append(list[i])

        return unique_list
    
    def list_contains_point(self, point, list):
        '''Returns True if the 3d-point is in the list, else False'''
        unique = False
        for i in range(len(list)):
            x = abs(list[i][0] - point[0])
            y = abs(list[i][1] - point[1])
            z = abs(list[i][2] - point[2])
            if x < 0.001 and y < 0.001 and z < 0.001:
                unique = True
                return unique
        return unique
    
    def list_contains_point_2d(self, point, list):
        '''Returns True if the 2d-point is in the list, else False'''
        unique = False
        for i in range(len(list)):
            x = abs(list[i][0] - point[0])
            y = abs(list[i][1] - point[1])
            if x < 0.001 and y < 0.001:
                unique = True
                return unique
        return unique

    def get_edge_intersection_points(self, triangle2, debug=False):
        """Returns the intersection points between the edges of the triangles"""
        det_values_1 = self.det_sign_check(self, triangle2)
        det_values_2 = triangle2.det_sign_check(triangle2, self)
        det_values_1, det_values_2, self, triangle2 = self.fix_permutation(det_values_1, det_values_2, self, triangle2, debug=debug)
        if debug:
            print("Triangle 1 after perm")
            self.print()
            print("det_values_1:", det_values_1)
            print("Triangle 2 after perm")
            triangle2.print()
            print("det_values_2:", det_values_2)

        i, j = self.isect_triangle_edges_plane_v3(self, triangle2.point2, triangle2.get_normal_vector(), debug=debug)
        if debug:
            print("intersection points [i, j]", [i,j])
        return [i, j]

    def create_non_intersecting_triangles(self, line_points, debug=False):
        '''Returns new triangles defineing the triangle based on a intersection line between to triangles
        
        line_points: colinear points defing the intersection line and conating the intersections with the triangle edges
        '''
        if debug:
            print("line_points", line_points)
        
        if (line_points[0] == None or line_points[1] == None):
            if debug:
                print("line points are None")
            return None
        
        self.ensure_counterclockwise()

        unique_line_points = self.get_unique_points_from_list(list(line_points))

        if debug:
            print("unique line points", unique_line_points)

        inside_points = [pt for pt in unique_line_points if self.is_point_in_triangle(pt, self.point1, self.point2, self.point3)]

        if debug:
            print("unique and inside points", inside_points)

        ordered_line_points = self.order_points_along_vector(inside_points, self.point3, self.point2)

        new_triangles = []
        if debug:
            print("len(ordered_line_points)", len(ordered_line_points))
            print("ordered_line_points", ordered_line_points)

        if len(ordered_line_points) < 2:
            if debug:
                print("len(ordered_line_points) less than 2")
            return None

        if len(ordered_line_points) >= 2:
            normal = True
            i_value, j_value = ordered_line_points[0], ordered_line_points[1]
            if self.list_contains_point(self.point2, ordered_line_points):
                normal = False
                if debug:
                    print("[i, j, k, l] contains point2")
                new_triangles.append(IfcTriangle(i_value, j_value, self.point3))
            if self.list_contains_point(self.point3, ordered_line_points):
                normal = False
                if debug:
                    print("[i, j, k, l] contains point3")
                new_triangles.append(IfcTriangle(i_value, j_value, self.point2))
            if normal:
                new_triangles.append(IfcTriangle(i_value, self.point2, self.point3))
                new_triangles.append(IfcTriangle(i_value, j_value, self.point2))

        last_value = ordered_line_points[-1]
        first_value = ordered_line_points[0]
        new_triangles.append(IfcTriangle(first_value, last_value, self.point1))

        return new_triangles
    
    def new_2d_triangles_needed(self, degree):
        """Returns False if the triangle slope is more than the degree"""
        z_component_normal = self.get_normal_vector()[2]
        degree = math.degrees(np.arcsin(z_component_normal))
        if abs(degree) > degree:
            return False
        return True

    def create_non_intersecting_triangles_2d(self, line_points, triangle, debug=False):
        """Returns new 2d-triangles based on the line splitting the original 2d-triangle"""
        def fix_permutation_2d_line(line, triangle):
            sign_1 = self.point_orientation_threshold(line[0], line[1], triangle[0])
            sign_2 = self.point_orientation_threshold(line[0], line[1], triangle[1])
            sign_3 = self.point_orientation_threshold(line[0], line[1], triangle[2])
            det_values_1 = [sign_1, sign_2, sign_3]

            if(det_values_1[0] >= 0 and det_values_1[1] <= 0 and det_values_1[2] <= 0) or (det_values_1[0] <= 0 and det_values_1[1] >= 0 and det_values_1[2] >= 0):
                return triangle
            
            point1 = triangle[0]
            triangle[0] = triangle[2]
            triangle[2] = triangle[1]
            triangle[1] = point1

            return fix_permutation_2d_line(line, triangle)

        if (line_points[0] == None or line_points[1] == None):
            if debug:
                print("line points are None")
            return None

        triangle = fix_permutation_2d_line(line_points, triangle)

        ordered_line_points = self.order_points_along_vector(line_points, triangle[1], triangle[2]) # funker dette uten permutation?

        new_triangles = []
        if debug:
            print("len(ordered_line_points)", len(ordered_line_points))
            print("ordered_line_points", ordered_line_points)

        if len(ordered_line_points) < 2:
            if debug:
                print("len(ordered_line_points) less than 2")
            return None

        if len(ordered_line_points) >= 2:
            normal = True
            i_value, j_value = ordered_line_points[0], ordered_line_points[1]
            if normal:
                new_triangles.append([i_value, triangle[1], triangle[2]])
                new_triangles.append([i_value, j_value, triangle[2]])

        last_value = ordered_line_points[-1]
        first_value = ordered_line_points[0]
        new_triangles.append([first_value, last_value, triangle[0]])

        return new_triangles

    def line_intersection(self, line1, line2, debug=False):
        """Returns true if the 2d-lines are intersecting, else False"""
        if debug:
            print("2d_line1:", line1)
            print("2d_line2:", line2)
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]

        div = det(xdiff, ydiff)
        if div == 0:
            return None

        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return [x, y]
    
    def compute_new_triangles_2d_v4(self, collison_element, debug=False):
        """Returns new triangles if the triangle is intersecting with the boundries of the element"""
        def point_inside_line(point, line):
            line_x_min = min(line[0][0], line[1][0])
            line_y_min = min(line[0][1], line[1][1])
            line_x_max = max(line[0][0], line[1][0])
            line_y_max = max(line[0][1], line[1][1])
            if point[0] > line_x_min and point[0] < line_x_max and point[1] > line_y_min and point[1] < line_y_max:
                return True
            return False
        
        origin, axis1, axis2, point1_2d = self.project_to_plane([self.point1], self.get_normal_vector(), debug=debug)
        origin, axis1, axis2, point2_2d = self.project_to_plane([self.point2], self.get_normal_vector(), origin=origin, debug=debug)
        origin, axis1, axis2, point3_2d = self.project_to_plane([self.point3], self.get_normal_vector(), origin=origin, debug=debug)
        triangle_2d = [point1_2d[0], point2_2d[0], point3_2d[0]]
        triangle_edges_2d = [[point1_2d[0], point2_2d[0]], [point1_2d[0], point3_2d[0]], [point2_2d[0], point3_2d[0]]]

        if debug:
            print("len(edges)", len(collison_element.edges))

        for edge in collison_element.edges:
            if not self.overlapping_with_edge(edge):
                continue
            edge_point1 = edge[0]
            edge_point2 = edge[1]
            if self.is_point_coplanar(edge_point1, debug=debug) and self.is_point_coplanar(edge_point2, debug=debug):
                if debug:
                    print("Both points are coplanar")
                intersection_points = []
                origin, axis1, axis2, edge_points_2d = self.project_to_plane([edge_point1, edge_point2], self.get_normal_vector(), origin=origin, debug=debug)
                for t_edge_2d in triangle_edges_2d:
                    intersection_point = self.line_intersection(t_edge_2d, edge_points_2d, debug=debug)
                    if intersection_point is not None and point_inside_line(intersection_point, t_edge_2d):
                        if debug:
                            print("Intersection point inside line found")
                            print("Intersection point: ", intersection_point)
                        intersection_points.append(intersection_point)
                if debug:
                    print("intersection_points: ", intersection_points)
                if len(intersection_points) == 2:
                    if debug:
                        print("NEW triangles are being computed")
                    new_triangles = self.create_non_intersecting_triangles_2d(intersection_points, triangle_2d, debug=debug)
                    if new_triangles is None:
                        if debug:
                            print("new_triangles are None for edge:", edge)
                        continue
                    return self.project_from_plan(new_triangles, origin, axis1, axis2, debug=debug)
        return None
    
    def collision_test(self, triangle2, compute_new_triangles=False, compute_new_triangles_2d=False, debug=False):
        '''Returns True if the two triangles are intersecting either in 2d or 3d, else False
        returns: boolean_collision, new_triangles1, new_triangles2, changed_2d, collision_2d
        
        triangle2: IfcTriangle
        compute_new_triangles: If we want new non intersecting triangles for 3d-collisions
        compute_new_triangles_2d: If we want new non intersecting triangles for 2d-collisions
        debug: enabling printing
        '''
        result = False
        collision_2d = False

        if not self.need_collision_test(triangle2, debug=debug):
            return result, None, None, False, collision_2d
        
        if debug:
            print("Triangle points inside collision_test")
            print("Triangle 1")
            self.print()
            print("Triangle 2")
            triangle2.print()

        det_values_1 = self.det_sign_check(self, triangle2)
        
        if debug:
            print("det_values_1: ", det_values_1)

        if det_values_1[0] > -0.0001 and det_values_1[0] < 0.0001 and det_values_1[1] > -0.0001 and det_values_1[1] < 0.0001 and det_values_1[2] > -0.0001 and det_values_1[2] < 0.0001:
            if debug:
                print("coplanar")
                print("2d-check needed")

            normal = self.get_normal_vector()

            def project_to_plane(points, normal, debug=False):
                axis1 = np.cross(normal, [1, 0, 0]) if normal[0] == 0 else np.cross(normal, [0, 1, 0])

                axis1 = axis1 / np.linalg.norm(axis1)
                axis2 = np.cross(normal, axis1)

                projected = [[np.dot(p, axis1), np.dot(p, axis2)] for p in points]
                return projected

            projected_triangle1 = project_to_plane([self.point1, self.point2, self.point3], normal)
            projected_triangle2 = project_to_plane([triangle2.point1, triangle2.point2, triangle2.point3], normal)

            projected_triangle1 = self.ensure_counterclockwise_2d(projected_triangle1)
            projected_triangle2 = self.ensure_counterclockwise_2d(projected_triangle2)

            if debug:
                print("2d determinant test")

            det_values = self.det_sign_2d_v2(projected_triangle1, projected_triangle2)
            det_values_2 = self.det_sign_2d_v2(projected_triangle2, projected_triangle1)
            
            if debug:
                print("projected_triangle1", projected_triangle1)
                print("projected_triangle2", projected_triangle2)
                print("2d_determinant: ", det_values)

            collision = False
            for i in range(3):
                if (det_values[i][0] > 0 and det_values[i][1] > 0 and det_values[i][2] > 0):
                    collision = True

            if debug:
                print("2d_determinant_2: ", det_values_2)
            
            if not collision:
                for i in range(3):
                    if (det_values_2[i][0] > 0.001 and det_values_2[i][1] > 0.001 and det_values_2[i][2] > 0.001):
                        collision = True

            if not collision:
                if debug:
                    print("edge intersection test")
                collision = self.check_triangle_intersection(projected_triangle1, projected_triangle2)
            if not collision:
                collision = self.bigger(projected_triangle1, projected_triangle2)
                if debug:
                    print("triangle inside check")
            if not collision:
                if debug:
                    print("similarity_test")
                collision = self.triangle_similarity_v2(projected_triangle1, projected_triangle2, threshold=0.00001, debug=debug)
            if collision:
                if debug:
                    print("collision")
                result = True
                collision_2d = True
                return result, None, None, False, collision_2d
            else:
                if debug:
                    print("no collision")
                result = False
                return result, None, None, False, collision_2d
            
        if (det_values_1[0] > -0.001 and det_values_1[1] > -0.001 and det_values_1[2] > -0.001) or (det_values_1[0] < 0.001 and det_values_1[1] < 0.001 and det_values_1[2] < 0.001):
            if debug:
                print("no_collision")
            result = False
            return result, None, None, False, collision_2d
            
        else:
            if debug:
                print("3d-check needed")

            det_values_2 = self.det_sign_check(triangle2, self)

            if (det_values_2[0] > -0.001 and det_values_2[1] > -0.001 and det_values_2[2] > -0.001) or (det_values_2[0] < 0.001 and det_values_2[1] < 0.001 and det_values_2[2] < 0.001):
                if debug:    
                    print("no_collision")
                result = False
                return result, None, None, False, collision_2d  

            det_values_1, det_values_2, self, triangle2 = self.fix_permutation(det_values_1, det_values_2, self, triangle2, debug=False)
            if debug:
                print("det_values_1", det_values_1)
                print("det_values_2", det_values_2)

            p1 = np.array(self.point1)
            q1 = np.array(self.point2)
            r1 = np.array(self.point3)
            p2 = np.array(triangle2.point1)
            q2 = np.array(triangle2.point2)
            r2 = np.array(triangle2.point3)

            one = np.array([1])

            matrix_p1q1p2q2 = np.array([np.append(p1, one), np.append(q1, one), np.append(p2, one), np.append(q2, one)])
            if debug:
                print("matrix_p1q1p2q2", matrix_p1q1p2q2)
            det_p1q1p2q2 = np.linalg.det(matrix_p1q1p2q2)

            matrix_p1r1r2p2 = np.array([np.append(p1, one), np.append(r1, one), np.append(r2, one), np.append(p2, one)])
            if debug:
                print("matrix_p1r1r2p2", matrix_p1r1r2p2)
            det_p1r1r2p2 = np.linalg.det(matrix_p1r1r2p2)

            if debug:
                print("det_p1q1p2q2: ", det_p1q1p2q2)
                print("det_p1r1r2p2: ", det_p1r1r2p2)

            if (det_p1q1p2q2 <= 0.0001 and det_p1r1r2p2 <= 0.0001) or (det_p1q1p2q2 > -0.0001 and det_p1r1r2p2 > -0.0001):
                if debug:
                    print("collision")
                result = True
                if compute_new_triangles:
                    points_on_line = self.get_edge_intersection_points(triangle2, debug=debug)
                    new_triangles = self.create_non_intersecting_triangles(points_on_line, debug=debug)
                    #new_triangles_2 = triangle2.create_non_intersecting_triangles(points_on_line)
                    return result, new_triangles, None, False, collision_2d
                return result, None, None, False, collision_2d 

            if debug:
                print("no_collision")
            result = False
            return result, None, None, False, collision_2d 
        
    def collisions_triangle_element(self, element2, debug=False):
        '''Returns a list with True and False based on the collisions with the elements in the list
        
        element: IfcElement
        debug: enabling printing
        '''
        def update_triangle(element, triangle1_def):
            point11 = [element.x_coords[triangle1_def[0]], element.y_coords[triangle1_def[0]], element.z_coords[triangle1_def[0]]]
            point12 = [element.x_coords[triangle1_def[1]], element.y_coords[triangle1_def[1]], element.z_coords[triangle1_def[1]]]
            point13 = [element.x_coords[triangle1_def[2]], element.y_coords[triangle1_def[2]], element.z_coords[triangle1_def[2]]]
            triangle1 = IfcTriangle(point11, point12, point13)
            return triangle1
        
        collisions = []

        j = 0
        while j < len(element2.faces):
            triangle2_def = element2.faces[j]
            triangle2 = update_triangle(element2, triangle2_def)
            if debug:
                print("(i)", (j))
            collision, new_triangles1, new_triangles2, changed_2d, collision_2d = self.collision_test(triangle2, debug=debug)
            j = j + 1
            collisions.append(collision)

        return collisions
    
    def formwork_needed(self, degree, collision, inside, element, inside_element, debug=False, slope_degree=60):
        """Returns true if the triangle needs formwork, else False
        
        degree: float
        collision: boolean
        inside: boolean
        element: IfcElement
        inside_element: IfcElement
        debug: enabling printing 
        """
        formwork = False
        half_formwork = False

        if not element.need_formwork:
            return False

        if not collision:
            if degree < slope_degree:
                if inside:
                    if element.z_min < inside_element.z_min:
                        if debug:
                            print("[1], formwork")
                        formwork = True
                        return formwork, half_formwork
                    else:
                        if debug:
                            print("[2], no formwork")
                        formwork = False
                        return formwork, half_formwork
                else:
                    if abs(degree) > 80 and abs(degree) < 95:
                            if debug:
                                print("[3], formwork")
                            formwork = True
                            return formwork, half_formwork
                    else:
                        if debug:
                            print("[3], formwork")
                        formwork = True
                        return formwork, half_formwork
            else:
                if debug:
                    print("[4], no formwork")
                formwork = False
                return formwork, half_formwork                

        if collision:
            if abs(degree) > 80 and abs(degree) < 95:
                if debug:
                    print("[5]")
                formwork = False
                return formwork, half_formwork
            if degree < slope_degree:
                if degree > -80:
                    half_formwork = True

                if debug:
                    print("[6]")
                formwork = True
                return formwork, half_formwork
            else:
                formwork = False
                if debug:
                    print("[7]")
                return formwork, half_formwork
        
        if debug:
            print("[11]")
        return formwork, half_formwork