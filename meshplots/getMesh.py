import pathlib
from typing import Union

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



def parseShellFile(pathToFolder: str) -> pd.DataFrame:

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

def makeFirstPlane(pointsArr:np.ndarray, vector=None) -> Plane:

    if len(pointsArr) >= 3:
        return Plane.best_fit(pointsArr)

    elif len(pointsArr) == 1 and vector != None:
        return Plane(pointsArr[0], vector)


    else:
        raise Exception(f'First plane couldnt be made with parameters: points:{pointsArr}  vector:{vector}')

def movePlaneToFurthestPoint(plane:Plane, pointsDf:pd.DataFrame) -> None:

    normal = plane.normal
    normalized = normal / np.linalg.norm(normal)
    planePoint = plane.point
    maxDistance = max( plane.distance_point_signed(np.array(vert)) for vert in pointsDf.values.tolist())

    newVec = normalized * maxDistance
    newPoint = planePoint + newVec

    plane = Plane(newPoint, normalized)
    newDist = [plane.distance_point_signed(np.array(vert)) for vert in pointsDf.values.tolist()]


    if min( newDist)  > 1:

        raise Exception(f'moving Plane has failed, closest distance larger than 1 ({min(newDist)})')

    return Plane(newPoint, newVec)


def getClosestPointsToPlane(plane:Plane, pointsDf:pd.DataFrame, cutoff:int=8) -> pd.DataFrame:

    distanceToPlane = [plane.distance_point(np.array(vert)) for vert in pointsDf.values.tolist()]
    print(min(distanceToPlane))


    return pointsDf[(distanceToPlane < np.linspace(cutoff,cutoff,len(pointsDf)))].copy()

def getPointPlaneCoordinates(plane:Plane, pointsDf:pd.DataFrame) -> pd.DataFrame:

    pointProjection = [plane.project_point(vert) for vert in pointsDf.values.tolist()]
    normalized = plane.normal / np.linalg.norm(plane.normal)
    xvec = np.array([1,0,0])
    referenceZero = pointProjection[0]


    localX = (np.cross(xvec, normalized)) / np.linalg.norm(np.cross(xvec, normalized))
    localY = (np.cross(normalized, localX))

    xcoords = [np.dot((projectedPoint - referenceZero) , localX) for projectedPoint in pointProjection]
    ycoords = [np.dot((projectedPoint - referenceZero) , localY) for projectedPoint in pointProjection]

    return pd.DataFrame({'x': xcoords, 'y': ycoords})

def computeDihedral(plane1:Plane, plane2:Plane) -> np.float64:
    return np.arccos(np.dot(plane1.normal, plane2.normal) / (np.linalg.norm(plane1.normal) * np.linalg.norm(plane2.normal)))


def exampleNylon():
    # loading shell file
    pathToFolder = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/nylon/output_NYL_rad_1'
    df1 = parseShellFile(pathToFolder)

    # making plane

    plane = makeFirstPlane(
        [np.array([28.841, 65.113, 58.105]), np.array([15.231, 30.442, 49.824]), np.array([3.007, 4.288, 41.563]),
         np.array([-2.145, 65.211, 39.752]), np.array([31.409, 1.629, 58.441])])

    # moving plane to furthest point
    #
    movePlaneToFurthestPoint(plane, df1)
    points = []
    for row in df1.iterrows():
        points.append(np.array([row[1][0], row[1][1], row[1][2]]))

    ####    #extraacting points certain distance from plane

    distanceToPlane = np.array([plane.distance_point(point) for point in points])
    print(
        f'mean distance to plane: {np.mean(distanceToPlane)}    mindistance: {np.min(distanceToPlane)}, maxdistance: {np.max(distanceToPlane)}')

    df2 = getClosestPointsToPlane(plane, df1)

    localCoords = getPointPlaneCoordinates(plane, df2)

    # 3d plotting for visuals
    df2['distanceToPlane'] = [plane.distance_point(vert) for vert in df2.values.tolist()]

    a, b, c, d = plane.cartesian()
    palneX, planeY, planeZ = plane_points(a, b, c, d, size=500)

    fig = px.scatter_3d(df2, x='x', y='y', z='z', color='distanceToPlane')
    fig.add_trace(go.Scatter3d(x=df1['x'], y=df1['y'], z=df1['z']))
    fig.update_traces(marker=dict(size=5))
    fig.add_trace(go.Surface(x=palneX, y=planeY, z=planeZ, opacity=0.6, showscale=False))
    # fig.add_trace(go.Scatter3d(x=projx, y=projy, z=projz ,opacity=0.6))
    fig.update_layout(scene_aspectmode='data')
    fig.show()

    # 2d plotting of a plane
    #
    twodfig = (px.scatter(pd.DataFrame(
        {'x': localCoords['x'], 'y': localCoords['y'], 'distanceToPlane': df2["distanceToPlane"].to_list()}), x='x',
                          y='y', color='distanceToPlane'))
    twodfig.show()






def exampleProtein():
    pathToFolder = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET/output_SURFMAP_1CEX_cutfv_4rev_aligned_bfactor'


    alignPoint = [-25.77, 1.56, -0.93]
    df1 = parseShellFile(pathToFolder)
    plane = makeFirstPlane([[df1['x'].mean(),df1['y'].mean(),df1['z'].mean()]], alignPoint)

    plane = movePlaneToFurthestPoint(plane, df1)

    df2 = getClosestPointsToPlane(plane, df1)
    newPlane = makeFirstPlane(df2.values.tolist())
    newPlane = movePlaneToFurthestPoint(newPlane, df2)

    print('dihedral:')
    print((computeDihedral(plane, newPlane)*180) / np.pi)

    a,b,c,d = plane.cartesian()
    x,y,z = newPlane.to_mesh([-5,5], [-5,5])


    fig = px.scatter_3d(df1, x='x', y='y', z='z')
    fig.add_trace(go.Scatter3d(x=[alignPoint[0]], y=[alignPoint[1]], z=[alignPoint[2]]))
    fig.add_trace(go.Surface(x=x, y=y, z=z, opacity=0.6, showscale=False))
    fig.add_trace(go.Scatter3d(x=df2['x'], y=df2['y'], z=df2['z']))

    fig.show()




if __name__ == "__main__":
    exampleProtein()