import pandas
import numpy as np
from math import floor, ceil
from scipy import fftpack, optimize
from scipy.special import erf
from skimage import filters
from collections.abc import Iterable
from numbers import Number
from functools import cached_property
from pandas import DataFrame

np.seterr(invalid='ignore')


class Crop:
    """ Special crop object which never takes data from outside the array, and returns the used extent too,
        together with an image showing how much of each pixel is within the extent
    """
    def __getitem__(self, n):
        array, *n = n
        if isinstance(n, slice):
            n = (n,)
        if isinstance(n, type(Ellipsis)):
            n = (None,) * array.ndim
        if isinstance(n, Number):
            n = (slice(n),)
        n = list(n)
        ell = [i for i, e in enumerate(n) if isinstance(e, type(Ellipsis))]
        if len(ell) > 1:
            raise IndexError('an index can only have a single ellipsis (...)')
        if len(ell):
            if len(n) > array.ndim:
                n.remove(Ellipsis)
            else:
                n[ell[0]] = None
                while len(n) < array.ndim:
                    n.insert(ell[0], None)
        while len(n) < array.ndim:
            n.append(None)

        for i, (e, s) in enumerate(zip(n, array.shape)):
            if e is None:
                e = slice(None)
            elif isinstance(e, Number):
                e = slice(e, e)
            start, stop, step = floor(e.start or 0), ceil(e.stop or -1), round(e.step or 1)
            if start < 0:
                start = 0
            elif start >= s:
                start = s
            if stop >= s:
                stop = s
            elif stop < 0:
                stop = 0
            n[i] = slice(start, stop, step)
        return np.vstack([(i.start, i.stop) for i in n]), array[tuple(n)]


crop = Crop()


def gauss2(I, p, s, bg, shape=None):
    """ Create an image of a Gaussian peak
        I:  peak intensity
        p:  peak location (x, y) or grid with peak at origin
        s:  peak width (x & y)
        bg: tilt (offset, x, y) or just offset only
        S:  image size (x, y)

        returns: 2D array

        wp@tl20220328
    """
    if shape is None:
        xv, yv = p
    else:
        xv, yv = np.meshgrid(*[np.arange(s) - q for s, q in zip(shape, p)])
    im = I * np.exp(-(xv**2/2/s**2 + yv**2/2/s**2)) / (2*np.pi*s**2)
    try:
        return im + bg[1]*xv + bg[2]*yv + bg[0]
    except Exception:
        return im + bg


def gauss3(I, p, s, sz, bg, shape=None):
    """ Create an image of a Gaussian peak
        I:  peak intensity
        p:  peak location (z, y, x) or grid with peak at origin
        s:  peak width (x & y)
        sz: peak width (z)
        bg: tilt (offset, z, y, x) or just offset only
        S:  image size (z, y, x)

        returns: 3D array

        wp@tl20200713
    """
    if shape is None:
        zv, yv, xv = p
    else:
        zv, yv, xv = np.meshgrid(*[np.arange(s) - q for s, q in zip(shape, p)], indexing='ij')
    im = I * np.exp(-zv ** 2 / (2 * sz ** 2) - yv ** 2 / (2 * s ** 2) - xv ** 2 / (2 * s ** 2)) / \
         ((2 * np.pi) ** 1.5 * s * s * sz)
    try:
        return im + bg[3]*xv + bg[2]*yv + bg[1]*zv + bg[0]
    except Exception:
        return im + bg


def invert(f, y, bracket=None, tol=1e-15):
    """ invert a real valued scalar injective function to find x @ f(y) """
    g = lambda x: (f(x) - y) ** 2
    r = optimize.minimize_scalar(g, bracket=bracket, tol=tol)
    if not r.success or r.fun > tol:
        raise UserWarning(f'Could not invert function {f} at {y}')
    return r.x


class GaussianMaskFit:
    def __init__(self, image, offset, coordinates, sigma, window_size, sigma_z=None, window_size_z=None,
                 optimize_location=True, optimize_width=False, subtract_background=None,
                 location_confidence=None, width_confidence=0.05, max_iterations=20,
                 max_distance=None, max_distance_z=None, reset_if_not_fit=True, correct_fit_window=True,
                 refine_fit=False, max_sigma=None, max_sigma_z=None, fit_weight_sigma=2):
        self.image = image.astype('float')
        self.ndim = image.ndim
        if self.ndim not in (2, 3):
            raise NotImplementedError
        self.offset = np.asarray(offset, 'float')
        self.coordinates = np.asarray(coordinates, 'float')
        self.initial_coordinates = self.coordinates.copy()
        self.ideal_sigma = sigma
        self.initial_sigma = sigma
        self.window_size = window_size
        self.initial_window_size = window_size
        self.window_factor = window_size / sigma
        self.max_distance, self.max_distance_z = max_distance, max_distance_z
        self.correct_fit_window = correct_fit_window
        self.refine_fit = refine_fit
        self.max_sigma = max_sigma
        self.max_sigma_z = max_sigma_z
        self.fit_weight_sigma = fit_weight_sigma
        if self.ndim == 3:
            self.ideal_sigma_z = sigma_z
            self.initial_sigma_z = sigma_z
            self.window_size_z = window_size_z
            self.initial_window_size_z = window_size_z
            self.window_factor_z = window_size_z / sigma_z
        if subtract_background is None:
            subtract_background = 'tilt'
        if location_confidence is None:
            location_confidence = 0.05 if self.ndim == 3 else 0.01
        self.iterations = 0
        self.fit = False
        self.status = 'max iterations' if optimize_location else 'location fixed'
        if optimize_location:
            for self.iterations in range(max_iterations):
                S = self.crop
                if S.size == 0:
                    self.status = 'window outside'
                    break
                S, self.background = self.get_background(S, subtract_background)
                self.SN = S * self.N
                new_coordinates = self.get_new_coordinates()
                if optimize_width:
                    new_sigma = max(np.sqrt(self.get_new_sigma() * self.ideal_sigma), self.initial_sigma)
                    new_sigma_z = max(np.sqrt(self.get_new_sigma_z() * self.ideal_sigma_z), self.initial_sigma_z) \
                        if self.ndim == 3 else None
                else:
                    new_sigma = self.ideal_sigma
                    new_sigma_z = self.ideal_sigma_z if self.ndim == 3 else None

                if self.too_far_from_initial_coordinates(2 * new_coordinates - self.coordinates):
                    self.status = 'too far'
                    break

                self.fit = abs(new_coordinates - self.coordinates).max() < location_confidence
                if optimize_width:
                    self.fit = self.fit and abs(new_sigma - self.ideal_sigma) < width_confidence
                    if self.ndim == 3:
                        self.fit = self.fit and abs(new_sigma_z - self.ideal_sigma_z) < width_confidence

                if self.fit:
                    self.coordinates = 2 * new_coordinates - self.coordinates
                    if optimize_width:
                        self.set_sigma_and_window(new_sigma, new_sigma_z)
                    self.status = 'success'
                    break
                elif not np.all(np.isfinite(self.coordinates)):
                    self.status = 'inf coordinates'
                    break
                else:
                    self.coordinates = new_coordinates
                if optimize_width:
                    self.set_sigma_and_window(new_sigma, new_sigma_z)

                self.sigma = new_sigma
                self.sigma_z = new_sigma_z

        if not self.fit and reset_if_not_fit:
            self.coordinates = self.initial_coordinates
            self.ideal_sigma = self.initial_sigma
            self.window_size = self.initial_window_size
            if self.ndim == 3:
                self.ideal_sigma_z = self.initial_sigma_z
                self.window_size_z = self.initial_window_size_z

        self.S = self.crop
        S, self.background = self.get_background(self.S, subtract_background)
        self.SN = S * self.N
        self.sigma = self.get_new_sigma() if self.fit else self.ideal_sigma
        if self.ndim == 3:
            self.sigma_z = self.get_new_sigma_z() if self.fit else self.ideal_sigma_z

        if self.refine_fit:
            if self.ndim == 3:
                iz, iy, ix = self.grid
                W = gauss3(1, (iz - self.coordinates[0], iy - self.coordinates[1], ix - self.coordinates[2]),
                           self.fit_weight_sigma * self.ideal_sigma, self.fit_weight_sigma * self.ideal_sigma_z, 0)

                def fun(p):
                    intensity, coordinates, sigma, sigma_z, background = np.exp(p[0]), p[1:4], p[4], p[5], p[6:]
                    T = gauss3(intensity, (iz - coordinates[0], iy - coordinates[1], ix - coordinates[2]),
                               sigma, sigma_z, background)
                    return np.sum(W * (self.S - T) ** 2)

                guess = (np.log(self.get_corrected_integrated_intensity(self.get_integrated_intensity())),
                         *self.coordinates, self.sigma, self.sigma_z, *self.background)
                bounds = ((0, None), *[(w[0], w[1]) for w in self.window], (1, self.max_sigma), (1, self.max_sigma_z),
                          (None, None), (None, None), (None, None), (None, None))
            else:
                ix, iy = self.grid
                W = gauss2(1, (ix - self.coordinates[0], iy - self.coordinates[1]),
                           self.fit_weight_sigma * self.ideal_sigma, 0)

                def fun(p):
                    intensity, coordinates, sigma, background = np.exp(p[0]), p[1:3], p[3], p[4:]
                    T = gauss2(intensity, (ix - coordinates[0], iy - coordinates[1]), self.ideal_sigma, background)
                    return np.sum(W * (self.S - T) ** 2)

                guess = (np.log(self.get_corrected_integrated_intensity(self.get_integrated_intensity())),
                         *self.coordinates, self.sigma, *self.background)
                bounds = ((0, None), *[(w[0], w[1]) for w in self.window], (1, self.max_sigma),
                          (None, None), (None, None), (None, None))
            fit_result = optimize.minimize(fun, guess, bounds=bounds, tol=1e-40)
            self.corrected_integrated_intensity = np.exp(fit_result.x[0])
            self.coordinates = fit_result.x[1:4]
            self.sigma = fit_result.x[4]
            self.sigma_z = fit_result.x[5]
            self.background = fit_result.x[6:]
            self.integrated_intensity = self.corrected_integrated_intensity * \
                                        self.ideal_sigma ** 2 * self.ideal_sigma_z / self.sigma ** 2 / self.sigma_z
        else:
            self.integrated_intensity = self.get_integrated_intensity()
            self.corrected_integrated_intensity = self.get_corrected_integrated_intensity(self.integrated_intensity)

        if self.ndim == 3:
            iz, iy, ix = self.grid
            background = self.S - (gauss3(self.corrected_integrated_intensity, (iz - self.coordinates[0],
                iy - self.coordinates[1], ix - self.coordinates[2]), self.sigma, self.sigma_z, 0))
            fun = lambda q: gauss3(q[0], (iz - q[1], iy - q[2], ix - q[3]), q[4], q[5], q[6:])
            peak_fun = lambda q: gauss3(q[0], (iz - q[1], iy - q[2], ix - q[3]), q[4], q[5], 0)
            background_fun = lambda q: (q[3] * (ix - self.coordinates[2]) + q[2] * (iy - self.coordinates[1]) +
                                        q[1] * (iz - self.coordinates[0]) + q[0])
            self.X2, dq, self.R2 = fminerr(fun, [self.corrected_integrated_intensity, *self.coordinates,
                                                 self.sigma, self.sigma_z, *self.background], self.S)
            self.R2_peak = fminerr(peak_fun, [self.corrected_integrated_intensity, *self.coordinates,
                                              self.sigma, self.sigma_z], S)[2]
            self.dpeak_intensity = self.peak_intensity * np.sqrt((dq[0] / self.corrected_integrated_intensity) ** 2 +
                                                            4 * (dq[4] / self.sigma) ** 2 + (dq[5] / self.sigma_z) ** 2)
        else:
            ix, iy = self.grid
            background = self.S - gauss2(self.corrected_integrated_intensity,
                                         (ix - self.coordinates[0], iy - self.coordinates[1]), self.sigma, 0)
            fun = lambda q: gauss2(q[0], (ix - q[1], iy - q[2]), q[3], q[4:])
            peak_fun = lambda q: gauss2(q[0], (ix - q[1], iy - q[2]), q[3], 0)
            background_fun = lambda q: q[0] + q[1] * (ix - self.coordinates[0]) + q[2] * (iy - self.coordinates[1])
            self.X2, dq, self.R2 = fminerr(fun, [self.corrected_integrated_intensity, *self.coordinates,
                                                 self.sigma, *self.background], self.S)
            self.R2_peak = fminerr(peak_fun, [self.corrected_integrated_intensity, *self.coordinates,
                                              self.sigma], S)[2]
            self.dpeak_intensity = self.peak_intensity * np.sqrt((dq[0] / self.corrected_integrated_intensity) ** 2 +
                                                                 4 * (dq[4] / self.sigma) ** 2)
        self.R2_background = fminerr(background_fun, self.background, background)[2]
        self.dcorrected_integrated_intensity = dq[0]
        self.dcoordinates = dq[1:self.ndim + 1]
        self.dsigma = dq[self.ndim + 1]
        if self.ndim == 3:
            self.dsigma_z = dq[5]
        self.dbackground = dq[2 * self.ndim:]

    def too_far_from_initial_coordinates(self, coordinates):
        return self.max_distance is not None and \
               ((self.ndim == 2 and
                 np.sum((coordinates - self.initial_coordinates) ** 2 / self.max_distance ** 2) > 1) or
                (self.ndim == 3 and self.max_distance_z is not None and
                 np.sum((coordinates - self.initial_coordinates) ** 2 /
                        np.array((self.max_distance_z, self.max_distance, self.max_distance)) ** 2) > 1) or
                (self.ndim == 3 and not self.max_distance_z and
                 np.sum((coordinates[1:] - self.initial_coordinates[1:]) ** 2 / self.max_distance ** 2) > 1))

    def get_fit_image(self):
        if self.ndim == 3:
            return gauss3(self.corrected_integrated_intensity, self.coordinates, self.sigma, self.sigma_z,
                          self.background, self.image.shape)
        else:
            return gauss2(self.corrected_integrated_intensity, self.coordinates, self.sigma,
                          self.background, self.image.shape)

    def get_residue(self):
        return self.image - self.get_fit_image()

    def set_sigma_and_window(self, new_sigma, new_sigma_z):
        self.ideal_sigma = new_sigma
        self.window_size = self.window_factor * self.ideal_sigma
        if self.ndim == 3:
            self.ideal_sigma_z = new_sigma_z
            self.window_size_z = self.window_factor_z * self.ideal_sigma_z

    @property
    def crop(self):
        if self.ndim == 3:
            self.window, S = crop[(self.image, *[slice(c - w / 2, c + w / 2)
                for c, w in zip(self.coordinates, (self.window_size_z, self.window_size, self.window_size))])]
            self.grid = np.meshgrid(*[range(*w) for w in self.window], indexing='ij')
            iz, iy, ix = self.grid
            self.N = gauss3(1, (iz - self.coordinates[0], iy - self.coordinates[1], ix - self.coordinates[2]),
                            self.ideal_sigma, self.ideal_sigma_z, 0)
        else:
            self.window, S = crop[(self.image, *[slice(c - self.window_size / 2, c + self.window_size / 2)
                                                               for c in self.coordinates[::-1]])]
            self.grid = np.meshgrid(*[range(*w) for w in self.window[::-1]], indexing='xy')
            ix, iy = self.grid
            self.N = gauss2(1, (ix - self.coordinates[0], iy - self.coordinates[1]), self.ideal_sigma, 0)
        return S

    def get_background(self, S, mode):
        if self.ndim == 3:
            iz, iy, ix = self.grid
            if mode == 'tilt':
                x = np.arange(*self.window[2]) - self.coordinates[2]
                bgx = np.polyfit(x, np.r_[S[:, 0], S[:, -1]].mean(0), 1)[0]
                y = np.arange(*self.window[1]) - self.coordinates[1]
                bgy = np.polyfit(y, np.r_[S[:, :, 0], S[:, :, -1]].mean(0), 1)[0]
                z = np.arange(*self.window[0]) - self.coordinates[0]
                bgz = np.polyfit(z, np.c_[S[:, 0], S[:, -1], S[:, 1:-1, 0], S[:, 1:-1, -1]].mean(1), 1)[0]
                S = S - (bgz * (iz - self.coordinates[0]) + bgy * (iy - self.coordinates[1]) +
                         bgx * (ix - self.coordinates[2]))
                # note: the corners are used twice, not sure why
                bg = np.mean(np.hstack((S[:, 0], S[:, -1], S[:, :, 0], S[:, :, -1])))
                return S - bg, np.r_[bg, bgz, bgy, bgx]
            elif mode == 'constant':
                bg = np.mean(np.hstack((S[:, 0], S[:, -1], S[:, :, 0], S[:, :, -1])))
                return S - bg, np.r_[bg, 0., 0., 0.]
            else:
                return S, np.r_[0., 0., 0., 0.]
        else:
            ix, iy = self.grid
            if mode == 'tilt':
                xy = np.arange(*self.window[1]) - self.coordinates[0]
                bgx = np.polyfit(np.tile(xy, 2), np.r_[S[0], S[-1]], 1)[0]
                xy = np.arange(*self.window[0]) - self.coordinates[1]
                bgy = np.polyfit(np.tile(xy, 2), np.r_[S[:, 0], S[:, -1]], 1)[0]
                S = S - (bgx * (ix - self.coordinates[0]) + bgy * (iy - self.coordinates[1]))
                # note: the corners are used twice, not sure why
                bg = np.mean(np.hstack((S[0], S[-1], S[:, 0], S[:, -1])))
                return S - bg, np.r_[bg, bgx, bgy]
            elif mode == 'constant':
                bg = np.mean(np.hstack((S[0], S[-1], S[:, 0], S[:, -1])))
                return S - bg, np.r_[bg, 0., 0.]
            else:
                return S, np.r_[0., 0., 0.]

    def get_new_coordinates(self):
        # halfway towards the peak
        new_coordinates = np.array([np.sum(i * self.SN) for i in self.grid]) / np.sum(self.SN)
        return new_coordinates if np.all(np.isfinite(new_coordinates)) else self.coordinates

    @staticmethod
    def sigma_to_moment(sigma, x_a, x_b):
        """ The theoretical result of int(x^2*SN, x=a..b) / int(SN, x=a..b) for a gaussian with certain sigma """
        if isinstance(sigma, Number) and sigma == 0:
            return 0
        sigma2 = sigma ** 2
        f = erf(x_a / sigma / np.sqrt(2)) - erf(x_b / sigma / np.sqrt(2))
        e_a = np.exp(-x_a ** 2 / sigma2 / 2)
        e_b = np.exp(-x_b ** 2 / sigma2 / 2)
        return sigma2 + np.sqrt(2 / np.pi) * sigma * (x_b * e_b - x_a * e_a) / f

    @staticmethod
    def moment_to_sigma(m, a, b):
        """ Calculate the sigma of SN from the result of int(x^2*SN, x=a..b) / int(SN, x=a..b) """
        if isinstance(m, Number):
            lim = (a ** 2 + a * b + b ** 2) / 3
            if m == lim:
                return np.inf
            elif m > lim:
                return np.nan
            elif m <= 0:
                return 0
            else:
                return abs(invert(lambda sigma: __class__.sigma_to_moment(sigma, a, b), m, (np.sqrt(m), lim)))
        else:
            return [__class__.moment_to_sigma(j, a, b) for j in m]

    @staticmethod
    def correct_sigma(sigma2, sigma2_n):
        """ Calculate the sigma of the S from the sigma's of SN and N """
        return np.clip(1 / np.sqrt(1/sigma2 - 1/sigma2_n) if sigma2_n > sigma2 > 0 else 1, 1, None)

    @staticmethod
    def correct_intensity_fw(I, *args):
        """ Correct I_i of SN for the part of the gaussian that was outside the window
            args: (sigma_x, x_a, x_b), (sigma_y, y_a, y_b), ...
        """
        return I * np.prod([2 / (erf(b/sigma/np.sqrt(2)) - erf(a/sigma/np.sqrt(2))) for sigma, a, b in args])

    def get_new_sigma(self):
        # moment analysis to get peak width, only works well when winSize(Z) > 3 * width
        sigma2 = [np.sum((i - c) ** 2 * self.SN) / np.sum(self.SN)
                  for i, c in zip(self.grid[-2:], self.coordinates[-2:])]
        if self.correct_fit_window:
            sigma = [self.moment_to_sigma(s, w[0] - c, w[1] - c)
                     for s, w, c in zip(sigma2, self.window[-2:], self.coordinates[-2:])]
            sigma2 = sigma[0] * sigma[1]
        else:
            sigma2 = np.sqrt(sigma2[0] * sigma2[1])
        return self.correct_sigma(sigma2, self.ideal_sigma ** 2)

    def get_new_sigma_z(self):
        w, i, c = self.window[0], self.grid[0], self.coordinates[0]
        sigma2 = np.sum(((i - c) ** 2) * self.SN) / np.sum(self.SN)
        if self.correct_fit_window:
            sigma2 = self.moment_to_sigma(sigma2, w[0] - c, w[1] - c) ** 2
        return self.correct_sigma(sigma2, self.ideal_sigma_z ** 2)

    def get_integrated_intensity(self):
        return np.sum(self.SN) / np.sum(self.N ** 2)

    def get_corrected_integrated_intensity(self, intensity):
        # correct the intensity when sigma is different from ideal_sigma and
        # when part of the spot falls outside the fit window
        sigma_sn = 1 / np.sqrt(1 / self.sigma ** 2 + 1 / self.ideal_sigma ** 2)
        if self.ndim == 3:
            sigma_sn_z = 1 / np.sqrt(1 / self.sigma_z ** 2 + 1 / self.ideal_sigma_z ** 2)
            if self.correct_fit_window:
                intensity = __class__.correct_intensity_fw(intensity, *[(s, w[0] - c, w[1] - c)
                    for s, w, c in zip((sigma_sn_z, sigma_sn, sigma_sn), self.window, self.coordinates)])
            return intensity * ((self.sigma ** 2 + self.ideal_sigma ** 2)
                                * np.sqrt(self.sigma_z ** 2 + self.ideal_sigma_z ** 2)) / \
                (2 * np.sqrt(2) * self.ideal_sigma ** 2 * self.ideal_sigma_z)
        else:
            if self.correct_fit_window:
                intensity = __class__.correct_intensity_fw(intensity, *[(s, w[0] - c, w[1] - c)
                    for s, w, c in zip((sigma_sn, sigma_sn), self.window, self.coordinates)])
            return intensity * (self.sigma ** 2 + self.ideal_sigma ** 2) / (2 * self.ideal_sigma ** 2)

    @cached_property
    def peak_intensity(self):
        if self.ndim == 3:
            return self.corrected_integrated_intensity / ((2 * np.pi) ** (3 / 2) * self.sigma ** 2 * self.sigma_z)
        else:
            return self.corrected_integrated_intensity / (2 * np.pi * self.sigma ** 2)

    @cached_property
    def data(self):
        if self.ndim == 3:
            return DataFrame(((self.integrated_intensity, self.corrected_integrated_intensity,
                               self.peak_intensity, *(self.coordinates + self.offset),
                               self.sigma, self.sigma_z,
                               *self.background, self.dcorrected_integrated_intensity,
                               self.dpeak_intensity, *self.dcoordinates, self.dsigma,
                               self.dsigma_z, *self.dbackground, self.R2, self.R2_peak,
                               self.R2_background, self.fit,
                               *(self.initial_coordinates + self.offset),
                               self.iterations + 1, self.window_size, self.window_size_z,
                               self.status),),
                             columns=self.columns(self.ndim))
        else:
            return DataFrame(((self.integrated_intensity, self.corrected_integrated_intensity,
                               self.peak_intensity, *(self.coordinates + self.offset), self.sigma,
                               *self.background, self.dcorrected_integrated_intensity,
                               self.dpeak_intensity, *self.dcoordinates, self.dsigma,
                               *self.dbackground, self.R2, self.R2_peak,
                               self.R2_background, self.fit,
                               *(self.initial_coordinates + self.offset),
                               self.iterations + 1, self.window_size, self.status),),
                             columns=self.columns(self.ndim))

    @staticmethod
    def columns(ndim):
        if ndim == 3:
            return ('I', 'Ii', 'Ip', 'z', 'y', 'x', 's', 'sz', 'O', 'Oz', 'Oy', 'Ox', 'dIi', 'dIp', 'dz', 'dy', 'dx',
                    'ds', 'dsz', 'dO', 'dOz', 'dOy', 'dOx', 'R2', 'R2_peak', 'R2_background', 'fit',
                    'z_ini', 'y_ini', 'x_ini', 'iter', 'w', 'wz', 'status')
        else:
            return ('I', 'Ii', 'Ip', 'x', 'y', 's', 'O', 'Ox', 'Oy', 'dIi', 'dIp', 'dx', 'dy', 'ds', 'dO', 'dOx', 'dOy',
                    'R2', 'R2_peak', 'R2_background', 'fit', 'x_ini', 'y_ini', 'iter', 'w', 'status')

    @staticmethod
    def empty_data(ndim):
        return pandas.DataFrame(columns=__class__.columns(ndim))


def fminerr(fun, a, y, dy=None, diffstep=1e-6):
    """ Error estimation of a fit

        Inputs:
        fun: function which was fitted to data
        a:   function parameters
        y:   ydata, this must be all finite
        dy:  errors on ydata, this must be all finite

        Outputs:
        chisq: Chi^2
        da:    error estimates of the function parameters
        R2:    R^2

        Example:
        x = np.array((-3,-1,2,4,5))
        a = np.array((2,-3))
        y = (15,0,5,30,50)
        fun = lambda a: a[0]*x**2+a[1]
        chisq, dp, adjusted R2 = fminerr(fun,p,y)

        adjusted from Matlab version by Thomas Schmidt, Leiden University
        wp@tl2020
    """
    eps = np.spacing(1)
    a = np.array(a).flatten()
    y = np.array(y).flatten()
    if dy is None:
        dy = np.ones(y.shape)
    elif not isinstance(dy, Iterable):
        dy *= np.ones(y.shape)
    else:
        dy = np.array(dy).flatten()
    nData = np.size(y)
    nPar = np.size(a)
    dy = 1 / (dy + eps)
    f0 = np.array(fun(a)).flatten()
    chisq = np.sum(((f0 - y) * dy) ** 2) / (nData - nPar)

    # calculate adjusted R^2
    sstot = np.sum((y - np.mean(y)) ** 2)
    ssres = np.sum((y - f0) ** 2)
    R2 = 1 - (ssres / (nData - nPar)) / (sstot / (nData - 1))

    # calculate derivatives
    deriv = np.zeros((nData, nPar))
    for i in range(nPar):
        ah = a.copy()
        ah[i] = a[i] * (1 + diffstep) + eps
        f = np.array(fun(ah)).flatten()
        deriv[:, i] = (f - f0) / (ah[i] - a[i]) * dy

    hesse = np.matmul(deriv.T, deriv)

    try:
        if np.linalg.matrix_rank(hesse) == np.shape(hesse)[0]:
            da = np.sqrt(chisq * np.diag(np.linalg.inv(hesse)))
        else:
            try:
                da = np.sqrt(chisq * np.diag(np.linalg.pinv(hesse)))
            except Exception:
                da = np.zeros(a.shape)
    except Exception:
        da = np.zeros(a.shape)
    return chisq, da, R2


sHS = fftpack.fftshift  # Swap half-spaces. sHS(matrix[, axes]). axes=all by default
def hS(m,axes=None):
    if axes==None: axes=range(np.ndim(m))
    elif isinstance(axes, int): axes=[axes]
    elif axes==[]: return m
    return hS(m.swapaxes(0, axes[-1])[:m.shape[axes[-1]]/2].swapaxes(0, axes[-1]), axes[:-1])


def sHSM(m,axes=None):
    if axes==None: axes=range(np.ndim(m))
    elif isinstance(axes, int): axes=[axes]
    m=m.swapaxes(0, axes[0]); max=m[1]+m[-1]; m=(m+max/2)%max-max/2; m=m.swapaxes(0, axes[0])
    return sHS(m, axes)


def bpass(im, r1=1., r2=1.7):
    # TODO: deal with exp underflow (the result of exp(-1000) is too small for floats)
    ker1x=np.exp(-(sHS(sHSM(np.r_[:im.shape[1]]))/r1)**2/2); ker1x/=np.sum(ker1x); fker1x=fftpack.fft(ker1x)
    ker1y=np.exp(-(sHS(sHSM(np.r_[:im.shape[0]]))/r1)**2/2); ker1y/=np.sum(ker1y); fker1y=fftpack.fft(ker1y)
    ker2x=np.exp(-(sHS(sHSM(np.r_[:im.shape[1]]))/r2)**2/2); ker2x/=np.sum(ker2x); fker2x=fftpack.fft(ker2x)
    ker2y=np.exp(-(sHS(sHSM(np.r_[:im.shape[0]]))/r2)**2/2); ker2y/=np.sum(ker2y); fker2y=fftpack.fft(ker2y)
    fim=fftpack.fftn(im)
    return fftpack.ifftn((fim*fker1x).T*fker1y-(fim*fker2x).T*fker2y).real.T


def bpass3D(im, r1=1., r2=1.7, rz1=1., rz2=1.7):
    return filters.gaussian(im, np.clip((r1, r1, rz1), 1, None), mode='mirror') - \
           filters.gaussian(im, np.clip((r2, r2, rz2), 1, None), mode='mirror')
