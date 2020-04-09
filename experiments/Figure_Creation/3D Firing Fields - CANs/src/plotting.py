import plotly.graph_objects as go

def scatter(X, w=1., mode="markers", ranges=None):
    fig = go.Figure(data=[go.Scatter3d(
        x=X[:,0], 
        y=X[:,1], 
        z=X[:,2],
        mode=mode,
        marker=dict(
            size=5,
            color= w,
            colorscale='Jet',   
            opacity=1.,
            line=dict(
                    color='black',
                    width=1
                )     
        ),
        line=dict(
            color='black',
            width=2
        ),
        opacity=1.)])

    fig.update_layout(title=go.layout.Title(
            text="",
            xref="paper",
            x=0
        ),
        scene_aspectmode='cube')
    
    if ranges is not None:
        fig.update_layout(scene = dict(
                     xaxis = dict(nticks=10, range=ranges[0],),
                     yaxis = dict(nticks=10, range=ranges[1],),
                     zaxis = dict(nticks=10, range=ranges[2],),))

    
    fig.show()



def plot_3d_iso(X, Y, Z, values, c, eps):

    fig = go.Figure(data=go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        lighting=dict(roughness=0.3, specular=.9, diffuse=1, ambient=0.1),
        lightposition=dict(x=45,y=15,z=25),
        value=values.flatten(),
        isomin=c,
        isomax=c+eps,
        colorscale='Jet',
        opacity=1.,
        caps=dict(x_show=True, y_show=True, z_show=True)
        ))


    fig.update_layout(scene = dict(
                        xaxis = dict(
                             backgroundcolor="rgba(0, 0,0,0.1)",
                             gridcolor="white",
                             showbackground=True,
                             zerolinecolor="white",),
                        yaxis = dict(
                            backgroundcolor="rgba(0, 0,0,0.1)",
                             gridcolor="white",
                            showbackground=True,
                            zerolinecolor="white"),
                        zaxis = dict(
                            backgroundcolor="rgba(0, 0,0,0.2)",
                             gridcolor="white",
                            showbackground=True,
                            zerolinecolor="white",),),
                        width=700,
                        margin=dict(
                        r=10, l=10,
                        b=10, t=10),
                        
                      )

    fig['layout'].update(
        scene=dict(camera= dict(
                            up=dict(x=0, y=0, z=1),
                            center=dict(x=0, y=0, z=0),
                            eye=dict(x=1.5, y=2, z=0.6)))
    )

    fig.show()