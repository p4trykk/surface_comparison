import pathlib
import numpy as np
import pandas as pd
import plotly.express as px
import skspatial.objects
from plotly import graph_objs as go
from sympy import Plane, Point3D
import pymol
from sympy.abc import u,t,v,r, x,y,z
from  skspatial.objects import Plane
import matplotlib.pyplot as plt
import matplotlib
import plotly



def parseShellFile(pathToFolder):

    Filepath = next((pathlib.Path(pathToFolder).joinpath('shells')).glob('*.csv'))

    df1 = pd.read_csv(Filepath, header=None)
    df1.columns = ['x', 'y', 'z']

    return df1



def plane_points(a, b, c, d, size=60, res=20, center=None):
    """
    Return a grid of points on the plane  a*x + b*y + c*z + d = 0.

    center : optional (3,) point that the patch should be centered on.
             If None, uses the point on the plane closest to the origin.
    """
    n = np.array([a, b, c], dtype=float)
    n /= np.linalg.norm(n)

    # A point on the plane
    if center is None:
        P0 = -(d / np.dot(n, n)) * n          # closest point to origin
    else:
        P0 = np.asarray(center, dtype=float)

    # Two orthonormal in-plane basis vectors
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, helper); u /= np.linalg.norm(u)
    v = np.cross(n, u)

    # Grid of parameters
    s = np.linspace(-size, size, res)
    S, T = np.meshgrid(s, s)

    X = P0[0] + S*u[0] + T*v[0]
    Y = P0[1] + S*u[1] + T*v[1]
    Z = P0[2] + S*u[2] + T*v[2]
    return X, Y, Z

#
# def projectPoint(plane, point):





if __name__ == "__main__":
    #loading shell file
    pathToFolder = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/nylon/output_NYL_rad_1'
    df1 = parseShellFile(pathToFolder)

    #making plane and moving it
    plane = Plane.best_fit([np.array([28.841,65.113,58.105]), np.array([15.231,30.442,49.824]), np.array([3.007,4.288,41.563]), np.array([-2.145,65.211,39.752]), np.array([31.409,1.629,58.441])])
    a,b,c,d = plane.cartesian()
    normal = plane.normal

    normal /= np.linalg.norm(normal)
    normal = normal*50

    plane.point = np.array([a + normal[0],b + normal[1],c + normal[2]])
    a,b,c,d = plane.cartesian()



####    #extraacting points certain distance from plane
    points = []
    for row in df1.iterrows():
        points.append(np.array([row[1][0], row[1][1], row[1][2]]))

    points = np.array(points)


    distanceToPlane = np.array([plane.distance_point(point) for point in points])
    print(f'mean distance to plane: {np.mean(distanceToPlane)}    mindistance: {np.min(distanceToPlane)}, maxdistance: {np.max(distanceToPlane)}' )


    df1['distanceToPlane'] = distanceToPlane

    df2 = df1[df1['distanceToPlane'] < 20]




#       projecting points onto a plane
    points = []
    for row in df2.iterrows():
        points.append(np.array([row[1][0], row[1][1], row[1][2]]))

    pointProjection = [plane.project_point(point) for point in points]



    projx = []
    projy = []
    projz = []
    print(f'df2 len {len(df2["distanceToPlane"])}')
    print(f'pointProjection: {len(pointProjection)}')
    for point in pointProjection:

        projx.append(point[0])
        projy.append(point[1])
        projz.append(point[2])

    # calculating local coordinates

    normalized = plane.normal
    normalized /= np.linalg.norm(normalized)
    xvec = np.array([1,0,0])
    referenceZero = pointProjection[0]

    localX = (np.cross(xvec, normalized)) / np.linalg.norm(np.cross(xvec, normalized))
    localY = (np.cross(normalized, localX))

    xcoords = [np.dot((projectedPoint - referenceZero) , localX) for projectedPoint in pointProjection]
    ycoords = [np.dot((projectedPoint - referenceZero) , localY) for projectedPoint in pointProjection]


    print(f'normalized: {normalized}, localx: {localX} localy: {localY}')
    print(f'len(xcoords): {len(xcoords)}')



    #3d plotting for visuals

    palneX,planeY,planeZ = plane_points(a, b, c, d, size=500)



    fig = px.scatter_3d(df2, x='x', y='y', z='z', color='distanceToPlane')
    fig.add_trace(go.Scatter3d( x=df1['x'], y=df1['y'], z=df1['z']))
    fig.update_traces(marker=dict(size=5))
    fig.add_trace(go.Surface(x=palneX, y=planeY, z=planeZ, opacity=0.6, showscale=False))
    fig.add_trace(go.Scatter3d(x=projx, y=projy, z=projz ,opacity=0.6))
    fig.update_layout(scene_aspectmode='data')
    fig.show()

# 2d plotting of a plane
#
    twodfig = (px.scatter(pd.DataFrame({'x':xcoords, 'y':ycoords, 'distanceToPlane':df2["distanceToPlane"].to_list()}) ,x='x', y='y', color='distanceToPlane'))
    twodfig.show()
