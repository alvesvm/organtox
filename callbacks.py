import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from app import app
from dash.dependencies import Input, Output, State

from results import layout_results
from organtox_backend import *

###################################################### Callbacks ######################################################

#### OrganTox Profile ####
@app.callback(
    Output('graph-profile', 'figure'),
    [Input('data-source-profile', 'value'), Input('descriptor-type-profile', 'value'),
     Input('slider-profile', 'value'), Input('similarity_df', 'data')])
def update_figure(source_profile, descriptor_profile, slider_value_profile, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.source.isin(source_profile)]
        df_1 = df_0[df_0.descriptor == descriptor_profile]
        df_2 = df_1[df_1.time == slider_value_profile]

        df_3 = df_2.groupby(['compound_name', 'parameter'], as_index=False).agg(
            {'outcome': 'max', 'similarity': 'median'})
        df_4 = df_3.sort_values(by='similarity', ascending=False).head(1000)

        df_pivot = df_4.groupby(['compound_name', 'parameter'],
                                sort=False,
                                as_index=True)['outcome'].median().unstack()
        sort_list = [k for k in param_list if k in df_pivot.columns.values]
        df_pivot = df_pivot[sort_list]

        df_to_plotly = {'x': df_pivot.columns.tolist(), 'y': df_pivot.index.tolist(), 'z': df_pivot.values}

        hovertext = list()
        z_text = np.where(df_to_plotly['z'] >= 0.5, 'Observed', df_to_plotly['z'])
        z_text = np.where(z_text == '0.0', 'Not observed', z_text)
        z_text = np.where(z_text == 'nan', 'No data', z_text)
        for yi, yy in enumerate(df_to_plotly['y']):
            hovertext.append(list())
            for xi, xx in enumerate(df_to_plotly['x']):
                hovertext[-1].append('Parameter: {}<br />Compound: {}<br />Effect: {}'.format(xx, yy, z_text[yi][xi]))

        fig = go.Figure(data=go.Heatmap(df_to_plotly,
                                        colorscale=[[0, '#f0e442'],
                                                    [0.25, '#f0e442'],
                                                    [0.5, '#f0e442'],
                                                    [0.5, '#093a75'],
                                                    [0.75, '#093a75'],
                                                    [1, '#093a75']],
                                        hoverinfo='text',
                                        text=hovertext))

        fig.update_layout(height=700,
                          showlegend=False,
                          transition_duration=500,
                          xaxis_type='category',
                          xaxis=dict(
                              tickangle=30,
                              tickmode='array',
                              tickvals=[x for x in sorted(df_pivot.columns.tolist())]),
                          margin=dict(t=30, b=30),
                          yaxis=dict(
                              title='Nearest Neighbors',
                              tickmode='array',
                              tickvals=[y for y in sorted(df_pivot.index.tolist())]),
                          yaxis_nticks=len(df_pivot.index.tolist())
                          )

        fig.update_traces(colorbar_title='Effect',
                          colorbar_tickmode='array',
                          colorbar_tickvals=[1, 0.75, 0.5, 0.25, 0],
                          colorbar_ticktext=['', 'Observed', '', 'Not observed', ''],
                          colorbar_len=0.15,
                          selector=dict(type='heatmap'))

        return fig

    except:
        fig = go.Figure().add_annotation(text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(height=700,
                         xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                         yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                         margin=dict(t=30,b=30),
                         plot_bgcolor='white')
        return fig


    #### Tissue Profile ####
@app.callback(
    Output('graph-weight', 'figure'),
    [Input('data-source-tissue', 'value'), Input('descriptor-type-tissue', 'value'),
     Input('tissue-type', 'value'), Input('slider-similarity-tissue', 'value'),
     Input('similarity_df', 'data')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Weight']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]
        df_5 = df_4.round(2)

        fig = px.scatter(df_5,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='RdYlBu',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'],
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': f'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )

        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )

        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0',
                                                '0.1',
                                                '0.2',
                                                '0.3',
                                                '0.4',
                                                '0.5',
                                                '0.6',
                                                '0.7',
                                                '0.8',
                                                '0.9',
                                                '1.0']
                            )

        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig

    #### Tissue MuDRA plot ####
@app.callback(
    Output('graph-mudra-tissue', 'figure'),
    [Input('data-source-tissue', 'value'), Input('tissue-type', 'value'),
     Input('slider-mudra-tissue', 'value'), Input('similarity_df', 'data')])
def update_mudra(source, tissue, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Weight']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == tissue]

        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                          ascending=False).head(slider_value)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                                  ascending=False).head(slider_value)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)

        fig = go.Figure()
        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )

        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )

        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )
                        ), margin=dict(t=30,b=30),

        )


        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig


    #### Clinical Hemistry Profile ####
@app.callback(
    Output('graph-clin-chem', 'figure'),
    [Input('data-source-clin-chem', 'value'), Input('descriptor-type-clin-chem', 'value'),
     Input('clin-chem-type', 'value'), Input('slider-similarity-clin-chem', 'value'),
     Input('similarity_df', 'data')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Clinical chemistry']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]

        fig = px.scatter(df_4,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='RdYlBu',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'],
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': 'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )

        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )

        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0',
                                                '0.1',
                                                '0.2',
                                                '0.3',
                                                '0.4',
                                                '0.5',
                                                '0.6',
                                                '0.7',
                                                '0.8',
                                                '0.9',
                                                '1.0']
                            )

        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig


    #### Clinical Chemistry MuDRA plot ####
@app.callback(
    Output('graph-mudra-clin-chem', 'figure'),
    [Input('data-source-clin-chem', 'value'), Input('clin-chem-type', 'value'),
     Input('slider-mudra-clin-chem', 'value'), Input('similarity_df', 'data')])
def update_mudra(source, clinchem, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Clinical chemistry']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == clinchem]

        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                          ascending=False).head(slider_value)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                                  ascending=False).head(slider_value)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)

        fig = go.Figure()

        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )

        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )

        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )
                        ), margin=dict(t=30,b=30),

        )


        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig


        #### Hematology Profile ####
@app.callback(
    Output('graph-hemato', 'figure'),
    [Input('data-source-hemato', 'value'), Input('descriptor-type-hemato', 'value'), Input('hemato-type', 'value'),
     Input('slider-similarity-hemato', 'value'), Input('similarity_df', 'data')])
def update_figure(source_tissue, descriptor_tissue, tissue, slider_value_tissue, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df_0 = df[df.parameter_type == 'Hematology']
        df_1 = df_0[df_0.source.isin(source_tissue)]
        df_2 = df_1[df_1.descriptor == descriptor_tissue]
        df_3 = df_2[df_2.parameter == tissue]
        df_4 = df_3[df_3.similarity.between(slider_value_tissue[0], slider_value_tissue[1])]

        fig = px.scatter(df_4,
                         x='time',
                         y='value',
                         color='similarity',
                         color_continuous_scale='RdYlBu',
                         symbol=f'compound_name',
                         range_color=[slider_value_tissue[0], slider_value_tissue[1]],
                         hover_data=['value', 'value_unit', 'dose'],
                         labels={'compound_name': 'Compound Name',
                                 'time': 'Days',
                                 'dose': 'Dose (mg/kg/day)',
                                 'value': 'Experiment Value',
                                 'value_unit': 'Value Unit',
                                 'similarity': 'Tanimoto'}
                        )

        fig.update_layout(showlegend=False,
                          transition_duration=500,
                          xaxis = dict(
                              tickmode = 'array',
                              tickvals = [x for x in sorted(df_4['time'].unique())]),
                          margin=dict(t=30,b=30)
                         )

        fig.update_coloraxes(colorbar_tickmode='array',
                             colorbar_tickvals=np.arange(0,1.1,0.1),
                             colorbar_ticktext=['0.0',
                                                '0.1',
                                                '0.2',
                                                '0.3',
                                                '0.4',
                                                '0.5',
                                                '0.6',
                                                '0.7',
                                                '0.8',
                                                '0.9',
                                                '1.0']
                            )

        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig



        #### Hematology MuDRA plot ####
@app.callback(
    Output('graph-mudra-hemato', 'figure'),
    [Input('data-source-hemato', 'value'), Input('hemato-type', 'value'),
     Input('slider-mudra-hemato', 'value'), Input('similarity_df', 'data')])
def update_mudra(source, hemato, slider_value, json_df):
    try:
        df = pd.read_json(json_df, orient='split')
        df0 = df[df.parameter_type == 'Hematology']
        df1 = df0[df0.source.isin(source)]
        df2 = df1[df1.parameter == hemato]

        df_morgan = df2[df2.descriptor == 'morgan'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)
        df_maccs = df2[df2.descriptor == 'maccs'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                          ascending=False).head(slider_value)
        df_atom_pair = df2[df2.descriptor == 'atom_pair'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                                  ascending=False).head(slider_value)
        df_avalon = df2[df2.descriptor == 'avalon'].groupby('InChIKey').first().sort_values(by='similarity',
                                                                                            ascending=False).head(slider_value)

        fig = go.Figure()
        fig.add_trace(go.Scatterpolargl(r = df_morgan.similarity,
                                        theta = np.arange(92, 179, 87/df_morgan.shape[0]),
                                        name = 'Morgan',
                                        marker=dict(size=10, color='mediumseagreen'),
                                        customdata=df_morgan['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r= df_maccs.similarity,
                                        theta = np.arange(2, 89, 87/df_maccs.shape[0]),
                                        name = 'MACCS',
                                        marker=dict(size=10, color='darkorange'),
                                        customdata=df_maccs['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_atom_pair.similarity,
                                        theta = np.arange(182, 269, 87/df_avalon.shape[0]),
                                        name = 'AtomPair',
                                        marker=dict(size=10, color='mediumpurple'),
                                        customdata=df_avalon['compound_name']
                                       )
                     )

        fig.add_trace(go.Scatterpolargl(r = df_avalon.similarity,
                                        theta = np.arange(272, 359, 87/df_atom_pair.shape[0]),
                                        name = 'Avalon',
                                        marker=dict(size=10, color = 'magenta'),
                                        customdata=df_atom_pair['compound_name']
                                       )
                     )

        fig.update_traces(mode='markers',
                          marker=dict(line_color='white', opacity=0.7),
                          hovertemplate='<br>'.join(['Compound Name: %{customdata}',
                                                     'Similarity: %{r}',
                                                    ])
                         )

        fig.update_layout(
            transition_duration=500,
            font_size = 16,
            showlegend = True,
            polar = dict(angularaxis = dict(showticklabels=True,
                                            showgrid=True,
                                            ticks='',
                                            direction = 'counterclockwise',
                                            tickmode='array',
                                            tickvals=[0, 90, 180, 270],
                                            ticktext=['','', '', ''],
                                            linewidth=1,
                                            showline=True,
                                            linecolor='black'
                                           ),
                         radialaxis = dict(tickfont_size=12,
                                           tickangle=30,
                                           ticks='outside',
                                           showline=True,
                                           range=[0, 1]
                                          )
                        ), margin=dict(t=30,b=30),

        )

        return fig
    except:
        fig = go.Figure().add_annotation(x=2.5, y=1.5, text='No data to display', font=dict(color='#EAEAEA', size=25),
                                         showarrow=False)

        fig.update_layout(xaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          yaxis={'showgrid': False, 'zeroline': False, 'showticklabels':False},
                          margin=dict(t=30,b=30),
                          plot_bgcolor='white')
        return fig