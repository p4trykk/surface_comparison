import pathlib
from enum import unique
from typing import Union

import numpy as np
import pandas as pd
import plotly.express as px
import skspatial.objects
from plotly import graph_objs as go
from sympy.abc import u,t,v,r, x,y,z
from  skspatial.objects import Plane
import matplotlib.pyplot as plt
import matplotlib
import plotly
import trimesh
from trimesh.ray.ray_pyembree import RayMeshIntersector



def parseShellFile(pathToFolder: str) -> pd.DataFrame:

    Filepath = next((pathlib.Path(pathToFolder).joinpath('shells')).glob('*.csv'))

    df1 = pd.read_csv(Filepath, header=None)
    df1.columns = ['x', 'y', 'z']

    return df1



import numpy as np

def plane_mesh_from_points(x, y, z, margin=0.1, n=30):
    """
    Fit a plane to scattered 3D points and return a Plotly-ready mesh grid.

    Parameters
    ----------
    x, y, z : array-like
        Coordinates of points lying on (or near) the same plane.
    margin : float
        Fractional padding added to the grid extent (0.1 = 10%).
    n : int
        Grid resolution per axis.

    Returns
    -------
    X, Y, Z : np.ndarray, each (n, n)
        Grids for go.Surface(x=X, y=Y, z=Z).
    """
    pts = np.column_stack([x, y, z]).astype(float)

    # 1. Fit the plane via SVD (total least squares).
    #    Centroid is a point on the plane; smallest singular vector is the normal.
    centroid = pts.mean(axis=0)
    _, _, Vt = np.linalg.svd(pts - centroid)
    normal = Vt[-1]                      # unit normal
    normal /= np.linalg.norm(normal)

    # 2. Build an in-plane orthonormal basis (e1, e2, normal).
    #    Pick any vector not parallel to the normal.
    a = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(a, normal)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(normal, a); e1 /= np.linalg.norm(e1)
    e2 = np.cross(normal, e1)            # already unit

    # 3. Project points onto the local (e1, e2) frame to find bounds.
    rel = pts - centroid
    u = rel @ e1
    v = rel @ e2
    du = (u.max() - u.min()) * margin
    dv = (v.max() - v.min()) * margin
    u_lim = (u.min() - du, u.max() + du)
    v_lim = (v.min() - dv, v.max() + dv)

    # 4. Sample a grid in local coords and map back to 3D.
    U, V = np.meshgrid(np.linspace(*u_lim, n), np.linspace(*v_lim, n))
    X = centroid[0] + U * e1[0] + V * e2[0]
    Y = centroid[1] + U * e1[1] + V * e2[1]
    Z = centroid[2] + U * e1[2] + V * e2[2]

    return X, Y, Z

def plane_points(a, b, c, d, size=60, res=20, center=None):
    """
    Return 2D grids of x, y, z coordinates for points on the plane
        a*x + b*y + c*z + d = 0

    Suitable for plotly's go.Surface(x=..., y=..., z=...).

    center : optional (3,) point that the patch should be centered on.
             If None, uses the point on the plane closest to the origin.

    Returns
    -------
    X, Y, Z : 2D numpy arrays of shape (res, res)
    """
    n_raw = np.array([a, b, c], dtype=float)
    norm = np.linalg.norm(n_raw)

    if norm == 0:
        raise ValueError("Plane normal (a, b, c) must be nonzero.")

    n = n_raw / norm
    c0 = d / norm  # normalized plane: n·x + c0 = 0

    # A point on the plane
    if center is None:
        P0 = -c0 * n  # closest point to origin
    else:
        center = np.asarray(center, dtype=float)
        # Project the given center onto the plane
        P0 = center - (np.dot(n, center) + c0) * n

    # Two orthonormal in-plane basis vectors
    helper = (
        np.array([1.0, 0.0, 0.0])
        if abs(n[0]) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )

    u = np.cross(n, helper)
    u /= np.linalg.norm(u)

    v = np.cross(n, u)
    v /= np.linalg.norm(v)

    # Grid of parameters
    s = np.linspace(-size / 2.0, size / 2.0, res)
    t = np.linspace(-size / 2.0, size / 2.0, res)

    S, T = np.meshgrid(s, t, indexing="xy")

    points = P0 + S[..., None] * u + T[..., None] * v

    X = points[..., 0].ravel().tolist()
    Y = points[..., 1].ravel().tolist()
    Z = points[..., 2].ravel().tolist()

    return X, Y, Z



def makeFirstPlane(pointsArr:np.ndarray, vector=None) -> Plane:

    if len(pointsArr) >= 3:
        return Plane.best_fit(pointsArr)

    elif len(pointsArr) == 1 and vector != None:
        return Plane(pointsArr[0], vector)


    else:
        raise Exception(f'First plane couldnt be made with parameters: points:{pointsArr}  vector:{vector}')

def movePlaneToFurthestPoint(plane:Plane, pointsDf:pd.DataFrame) -> Plane:

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


def getClosestPointsToPlane(plane:Plane, pointsDf:pd.DataFrame, cutoff:int=12) -> pd.DataFrame:

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

def computeAllDihedrals():

    folderPath = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET'

    path = pathlib.Path(folderPath).glob('output*')
    alignPoint = [-25.77, 1.56, -0.93]
    computedDihedrals = []
    protnames = []
    plane1list = []
    plane2list = []
    verts = []
    chosenPoints = []

    for folder in path:
        df1 = parseShellFile(str(folder))

        plane1 = makeFirstPlane([df1.mean().values.tolist()], alignPoint)
        plane1 = movePlaneToFurthestPoint(plane1, df1)

        df2 = getClosestPointsToPlane(plane1, df1)
        plane2 = makeFirstPlane(df2.values.tolist())
        plane2 = movePlaneToFurthestPoint(plane2, df2)

        verts.append(df1)
        chosenPoints.append(df2)
        plane1list.append(plane1)
        plane2list.append(plane2)
        protnames.append(folder.name[15:])
        computedDihedrals.append((computeDihedral(plane1, plane2)*180)/ np.pi)

    dihedraldf = pd.DataFrame({'protein': protnames, 'dihedral': computedDihedrals, 'plane1': plane1list, 'plane2': plane2list, 'verticies': verts, 'closestPoints':chosenPoints})

    print(f'max dihedral: {max(computedDihedrals)}, min dihedral: {min(computedDihedrals)}, avg: {sum(computedDihedrals)/len(computedDihedrals)}')
    print(f"planes dihedral over 10: {dihedraldf[dihedraldf['dihedral'] >= 10]['protein']}")

    weirddf = dihedraldf[dihedraldf['dihedral'] <= 10]
    iterator = 0
    for row in weirddf.iterrows():
        fig = px.scatter_3d(row[1]['verticies'], x='x', y='y', z='z', title=f"{row[1]['protein']}  theta: {row[1]['dihedral']}")
        a,b,c,d = row[1]['plane1'].cartesian()
        a2,b2,c2,d2 = row[1]['plane2'].cartesian()

        x1,y1,z1 = plane_points(a,b,c,d, size=20)
        x2,y2,z2 = plane_points(a2,b2,c2,d2, size=20)
        print(x1,y1,z1)

        fig.add_trace(go.Scatter3d(x=row[1]['closestPoints']['x'], y= row[1]['closestPoints']['y'], z= row[1]['closestPoints']['z'], mode='markers'))
        fig.add_trace(go.Surface(x=x1, y=y1, z=z1, name='first Plane'))
        fig.add_trace(go.Surface(x=x2, y=y2, z=z2, opacity=0.6, showscale=False, name=f'second Plane, theta= {row[1]["dihedral"]}'))
        fig.show()

        if iterator == 10:
            break
        else :
            iterator += 1


def raytracing(offFilePath=None):
    file = '/home/bob/PyCharmMiscProject/surface_comparison/meshplots/xyzFiles/1CEX_cutfv_4rev_aligned.off'
    mesh = trimesh.load_mesh(file)
    # trimesh.Scene(mesh).show()

    pathToFolder = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET/output_SURFMAP_1CEX_cutfv_4rev_aligned_bfactor'
    alignPoint = [-25.77, 1.56, -0.93]
    df1 = parseShellFile(str(pathToFolder))

    plane = makeFirstPlane([[df1['x'].mean(),df1['y'].mean(),df1['z'].mean()]], alignPoint)

    plane = movePlaneToFurthestPoint(plane, df1)

    a,b,c,d = plane.cartesian()
    rayx,rayy,rayz = plane_points(a,b,c,d, size=100, res=110)
    rayDirection = plane.normal * -1


    intersector = RayMeshIntersector(mesh)
    rayCoords = [[x,y,z] for x,y,z in zip(rayx,rayy,rayz)]
    dirList = [rayDirection for i in range(len(rayCoords))]

    faces = intersector.intersects_first(ray_origins=rayCoords, ray_directions=dirList)
    faces = pd.Series(faces)[faces >= 1]
    meshdict = mesh.to_dict()

    verticeID = []

    for hit in faces.values.tolist():
        verticeID += meshdict['faces'][hit]


    verticeID = pd.Series(verticeID).unique().tolist()

    df2 = df1.iloc[verticeID]


    x,y,z = plane_mesh_from_points(rayx, rayy, rayz)
    fig = px.scatter_3d(df1, x='x',y='y',z='z', title='Ray Direction')
    fig.add_trace(go.Scatter3d(x=df2['x'],y=df2['y'],z=df2['z'], mode='markers'))
    fig.add_trace(go.Surface(x=x,y=y,z=z, showscale=False, name='Ray Direction'))
    fig.add_trace(go.Scatter3d(x=rayx,y=rayy,z=rayz, mode='markers'))
    fig.show()



def makeOFFfile(xyzPath:str, facePath:str, outputPath:str) -> None:

    name = pathlib.Path(xyzPath).name

    with open(xyzPath, 'r') as f:
        df = pd.read_csv(f, sep=',', header=None)


    newlines = df.values.tolist()
    newlines = [f'{line[0]} {line[1]} {line[2]}' + '\n' for line in newlines]

    verticeNum = len(df)

    with open(facePath, 'r') as f:

        lines = f.readlines()

        newline = []
        for line in lines[3:]:
            nowline = line.split(' ')

            vertid = [3]
            i=1
            for item in nowline:

                if item != '' and i <=3:
                    vertid.append(int(item) -1 )
                    i +=1

            newline.append(vertid)

        faceNum = len(newline)
        OFFfaceLines = [rf'{item[0]} {item[1]} {item[2]} {item[3]}' + '\n' for item in newline]

        newfile =  ['OFF\n', f'{verticeNum} {faceNum} 0' + '\n'] + newlines + OFFfaceLines

    with open(outputPath + '/' + name[:-4] + '.off', 'w') as f:
        f.writelines(newfile)



if __name__ == "__main__":
    xyzfile = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET/output_SURFMAP_1CEX_cutfv_4rev_aligned_bfactor/shells/1CEX_cutfv_4rev_aligned.csv'
    facefile = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET/output_SURFMAP_1CEX_cutfv_4rev_aligned_bfactor/shells/1CEX_cutfv_4rev_aligned.face'
    makeOFFfile(xyzfile, facefile, '/home/bob/PyCharmMiscProject/surface_comparison/meshplots/xyzFiles')
    offFile = '/home/bob/PyCharmMiscProject/surface_comparison/meshplots/xyzFiles/1CEX_cutfv_4rev_aligned.off'

    raytracing(offFile)




