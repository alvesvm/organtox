import dash
import dash_bootstrap_components as dbc
#from jupyter_dash import JupyterDash

#JupyterDash.infer_jupyter_proxy_config()

external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, suppress_callback_exceptions=True, external_stylesheets = external_stylesheets)
server = app.server