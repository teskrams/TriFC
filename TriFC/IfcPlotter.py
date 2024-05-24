import ifcopenshell.geom
import ifcopenshell.util.shape
import numpy as np
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