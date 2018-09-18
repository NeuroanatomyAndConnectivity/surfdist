def viz(coords, faces, stat_map=None,
        elev=0, azim=0, cmap='coolwarm',
        threshold=None, alpha='auto',
        bg_map=None, bg_on_stat=False,
        figsize=None,
        **kwargs):

    ''' Visualize results on cortical surface using matplotlib.

    Inputs
    -------
    coords : numpy array of shape (n_nodes,3), each row specifying the x,y,z
            coordinates of one node of surface mesh
    faces : numpy array of shape (n_faces, 3), each row specifying the indices
            of the three nodes building one node of the surface mesh
    stat_map : numpy array of shape (n_nodes,) containing the values to be
               visualized for each node.
    elev, azim : integers, elevation and azimuth parameters specifying the view
                 on the 3D plot. For Freesurfer surfaces elev=0, azim=0 will
                 give a lateral view for the right and a medial view for the
                 left hemisphere, elev=0, azim=180 will give a medial view for
                 the right and lateral view for the left hemisphere.
    cmap : Matplotlib colormap, the color range will me forced to be symmetric.
           Colormaps can be specified as string or colormap object.
    threshold : float, threshold to be applied to the map, will be applied in
                positive and negative direction, i.e. values < -abs(threshold)
                and > abs(threshold) will be shown.
    alpha : float, determines the opacity of the background mesh, in'auto' mode
            alpha defaults to .5 when no background map is given, to 1 otherwise.
    bg_map : numpy array of shape (n_nodes,) to be plotted underneath the
             statistical map. Specifying a sulcal depth map as bg_map results
             in realistic shadowing of the surface.
    bg_on_stat : boolean, specifies whether the statistical map should be
                 multiplied with the background map for shadowing. Otherwise,
                 only areas that are not covered by the statsitical map after
                 thresholding will show shadows.
    figsize : tuple of intergers, dimensions of the figure that is produced.


    Output
    ------
    Matplotlib figure object
    '''

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    from mpl_toolkits.mplot3d import Axes3D

    # load mesh and derive axes limits
    faces = np.array(faces, dtype=int)
    limits = [coords.min(), coords.max()]

    # set alpha if in auto mode
    if alpha == 'auto':
        if bg_map is None:
            alpha = .5
        else:
            alpha = 1

    # if cmap is given as string, translate to matplotlib cmap
    if type(cmap) == str:
        cmap = plt.cm.get_cmap(cmap)

    # initiate figure and 3d axes
    if figsize is not None:
        fig = plt.figure(figsize=figsize)
    else:
        fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', xlim=limits, ylim=limits)
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()

    # plot mesh without data
    p3dcollec = ax.plot_trisurf(coords[:, 0], coords[:, 1], coords[:, 2],
                                triangles=faces, linewidth=0.,
                                antialiased=False,
                                color='white')

    # If depth_map and/or stat_map are provided, map these onto the surface
    # set_facecolors function of Poly3DCollection is used as passing the
    # facecolors argument to plot_trisurf does not seem to work
    if bg_map is not None or stat_map is not None:

        face_colors = np.ones((faces.shape[0], 4))
        face_colors[:, :3] = .5*face_colors[:, :3]

        if bg_map is not None:
            bg_data = bg_map
            if bg_data.shape[0] != coords.shape[0]:
                raise ValueError('The bg_map does not have the same number '
                                 'of vertices as the mesh.')
            bg_faces = np.mean(bg_data[faces], axis=1)
            bg_faces = bg_faces - bg_faces.min()
            bg_faces = bg_faces / bg_faces.max()
            face_colors = plt.cm.gray_r(bg_faces)

        # modify alpha values of background
        face_colors[:, 3] = alpha*face_colors[:, 3]

        if stat_map is not None:
            stat_map_data = stat_map
            stat_map_faces = np.mean(stat_map_data[faces], axis=1)

            # Ensure symmetric colour range, based on Nilearn helper function:
            # https://github.com/nilearn/nilearn/blob/master/nilearn/plotting/img_plotting.py#L52
            vmax = max(-np.nanmin(stat_map_faces), np.nanmax(stat_map_faces))
            vmin = -vmax

            if threshold is not None:
                kept_indices = np.where(abs(stat_map_faces) >= threshold)[0]
                stat_map_faces = stat_map_faces - vmin
                stat_map_faces = stat_map_faces / (vmax-vmin)
                if bg_on_stat:
                    face_colors[kept_indices] = cmap(stat_map_faces[kept_indices]) * face_colors[kept_indices]
                else:
                    face_colors[kept_indices] = cmap(stat_map_faces[kept_indices])
            else:
                stat_map_faces = stat_map_faces - vmin
                stat_map_faces = stat_map_faces / (vmax-vmin)
                if bg_on_stat:
                    face_colors = cmap(stat_map_faces) * face_colors
                else:
                    face_colors = cmap(stat_map_faces)

        p3dcollec.set_facecolors(face_colors)

    return fig, ax
