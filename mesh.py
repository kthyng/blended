"""
Code for editing masks, from Rob's pygridgen
"""

import matplotlib.pyplot as plt
import numpy as np


class edit_mask_mesh(object):

    def _on_key(self, event):
        if event.key == 'e':
            self._clicking = not self._clicking
            plt.title('Editing %s -- click "e" to toggle' % self._clicking)
            plt.draw()

    def _on_click(self, event):
        x, y = event.xdata, event.ydata
        if event.button==1 and event.inaxes is not None and self._clicking == True:
            d = (x-self._xc)**2 + (y-self._yc)**2
            if isinstance(self.xv, np.ma.MaskedArray):
                idx = np.argwhere(d[~self._xc.mask] == d.min())
            else:
                idx = np.argwhere(d.flatten() == d.min())
            self._mask[idx] = float(not self._mask[idx])
            i, j = np.argwhere(d == d.min())[0]
            self.mask[i, j] = float(not self.mask[i, j])
            self._pc.set_array(self._mask)
            self._pc.changed()
            plt.draw()

    def __init__(self, xv, yv, mask, **kwargs):
        if xv.shape != yv.shape:
            raise ValueError('xv and yv must have the same shape')
        for dx, dq in zip(xv.shape, mask.shape):
            if dx != dq+1:
                raise ValueError('xv and yv must be cell verticies '
                                 '(i.e., one cell bigger in each dimension)')

        self.xv = xv
        self.yv = yv

        self.mask = mask

        land_color = kwargs.pop('land_color', (0.6, 1.0, 0.6))
        sea_color = kwargs.pop('sea_color', (0.6, 0.6, 1.0))

        cm = plt.matplotlib.colors.ListedColormap([land_color, sea_color],
                                                 name='land/sea')
        self._pc = plt.pcolor(xv, yv, mask, cmap=cm, vmin=0, vmax=1, **kwargs)
        self._xc = 0.25*(xv[1:,1:]+xv[1:,:-1]+xv[:-1,1:]+xv[:-1,:-1])
        self._yc = 0.25*(yv[1:,1:]+yv[1:,:-1]+yv[:-1,1:]+yv[:-1,:-1])

        if isinstance(self.xv, np.ma.MaskedArray):
            self._mask = mask[~self._xc.mask]
        else:
            self._mask = mask.flatten()

        plt.connect('button_press_event', self._on_click)
        plt.connect('key_press_event', self._on_key)
        self._clicking = False
        plt.title('Editing %s -- click "e" to toggle' % self._clicking)
        plt.draw()
