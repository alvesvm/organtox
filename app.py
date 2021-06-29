import dash

import dash_bootstrap_components as dbc

external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__,
                suppress_callback_exceptions = True,
                prevent_initial_callbacks=True,
                external_stylesheets = external_stylesheets)
server = app.server
