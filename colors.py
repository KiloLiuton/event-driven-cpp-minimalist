import itertools as it
import numpy as np
import colorsys


def _get_colors(num_colors, var=10):
    """Return a list of visually distinct colors by uniformly spreading points
    in hue space and selecting normally distributed lightness and saturation.

    arguments:
    num_colors -- the number of distinct colors you wish to generate
    var = 10   -- random variance in lightness and saturation (max 10, min -50)
    """
    if num_colors <= 0:
        return (255, 0, 0)
    colors = []
    if var > 10:
        var = 10
    elif var < -50:
        var = -50
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * var)/100.
        saturation = (90 + np.random.rand() * var)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors


class markerStyles:
    def __init__(self):
        self.thick_markers = [
                'o', 'v', '^', '<', '>', '8', 's', 'p', 'P', '*', 'h', 'H',
                'X', 'D', 'd']
        self.thin_markers = ['.', ',', '1', '2', '3', '4', '+', 'x', '|', '_']
        self.thick_gen = it.cycle(self.thick_markers)
        self.thin_gen = it.cycle(self.thin_markers)
        self.random_thick = it.cycle(np.random.choice(
            self.thick_markers,
            len(self.thick_markers),
            replace=False))
        self.random_thin = it.cycle(np.random.choice(
            self.thin_markers,
            len(self.thin_markers),
            replace=False))
        self.random_all = it.cycle(np.random.choice(
            self.thick_markers + self.thin_markers,
            len(self.thick_markers) + len(self.thin_markers),
            replace=False))

    def next_thick(self):
        '''Return the next thick marker in sequence'''
        return next(self.thick_gen)

    def next_thin(self):
        '''Return the next thick marker in sequence'''
        return next(self.thin_gen)

    def rand_thick(self):
        return next(self.random_thick)

    def rand_thin(self):
        return next(self.random_thin)

    def rand(self):
        return next(self.random_all)
