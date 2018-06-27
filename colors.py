import numpy as np
import colorsys

def _get_colors(num_colors):
    """Return a list of visually distinct colors by uniformly spreading points
    in hue space and selecting normally distributed lightness and saturation.

    arguments:
    num_colors -- the number of distinct colors you wish to generate
    """
    if num_colors == 0:
        return (255, 0, 0)
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors
