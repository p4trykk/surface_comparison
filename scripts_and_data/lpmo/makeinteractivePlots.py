import plotly
import pathlib
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.express as px

def parseMatrix(path):

    with open(path) as f:
        parsedFile = [line.rstrip('\n').split('\t') for line in f.readlines()]

    df = pd.DataFrame(parsedFile[1::], columns=['x', 'y', 'value', 'residue'])
    df['x'] = df['x'].astype(int) -180
    df['y'] = df['y'].astype(int) -90
    df['value'] = [float(x) if x not in ['NA', 'Inf'] else np.nan for x in df['value'] ]

    return df

base_dir = '/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/lpmo/out'

for dir in pathlib.Path(base_dir).iterdir():

    dirname = dir.name
    print(dirname)
    if dirname.startswith('output'):
        mapMatrix = base_dir + f'/{dirname}/matrices/{dirname[15::]}_matrix.txt'
        print(str(mapMatrix) + '\n\n\n')

        #CHANGE HERE to matrix path in order to show 1 plot!!!!
        matrixDf = parseMatrix('/home/bob/PyCharmMiscProject/surface_comparison/scripts_and_data/output_fixedNEW_PROTEIN_SET/output_SURFMAP_1CEX_cutfv_4rev_aligned_bfactor/matrices/1CEX_cutfv_4rev_aligned_bfactor_matrix.txt')
        # matrixDf.dropna(subset=['value'], inplace=True)

        plot = px.scatter(matrixDf, x='x', y='y', color='value', color_continuous_scale=[
                      (0, 'rgb(255, 11, 11)'),      # Lowest values = Red (position 0)
                      (0.5, 'rgb(254, 254, 255)'),  # Middle = White (position 0.5)
                      (1, 'rgb(8, 8, 255)')         # Highest values = Blue (position 1)
                  ], hover_data=['residue'], title=dirname)
        plot.update_traces(marker=dict(size=4, symbol="square"))
        plot.update_layout(plot_bgcolor='gray', xaxis=dict(showgrid=False),
    yaxis=dict(showgrid=False))
        plot.show(config={'showAxisDragHandles': False})

    break

