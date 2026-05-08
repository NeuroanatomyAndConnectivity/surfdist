import matplotlib

matplotlib.use('Agg')

import matplotlib.figure  # noqa: E402
import numpy as np  # noqa: E402

from surfdist.viz import viz  # noqa: E402


def test_viz_returns_figure(grid_surface):
    verts, tris, _ = grid_surface
    stat = np.linspace(0, 1, len(verts))
    fig, ax = viz(verts, tris, stat_map=stat)
    assert isinstance(fig, matplotlib.figure.Figure)


def test_viz_accepts_string_cmap(grid_surface):
    """Regression test for the matplotlib.colormaps[name] migration."""
    verts, tris, _ = grid_surface
    stat = np.linspace(-1, 1, len(verts))
    fig, _ = viz(verts, tris, stat_map=stat, cmap='coolwarm')
    assert isinstance(fig, matplotlib.figure.Figure)


def test_viz_with_bg_map(grid_surface):
    verts, tris, _ = grid_surface
    stat = np.linspace(0, 1, len(verts))
    bg = np.linspace(-1, 1, len(verts))
    fig, _ = viz(verts, tris, stat_map=stat, bg_map=bg, bg_on_stat=True)
    assert isinstance(fig, matplotlib.figure.Figure)
