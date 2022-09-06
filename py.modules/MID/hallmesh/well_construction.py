# -*- coding: utf-8 -*-
import numpy as np


class Tube:
    __scope__ = ['id', 'd', 'th', 'th0', 'sigma', 'mu']

    def __init__(self, *args, **kwargs):
        self.id = kwargs.get('id', '')
        self.d = kwargs.get('d', None)
        self.th = kwargs.get('th', None)
        self.th0 = kwargs.get('th0', self.th)
        self.sigma = kwargs.get('sigma', 4.7E+6)
        self.mu = kwargs.get('mu', 65.9)
        # self.sigma = kwargs.get('sigma', 2E+6)
        # self.mu = kwargs.get('mu', 10)

        self.top = None
        self.bottom = None
        self.joints = []
        self.joint_width = None

    def is_valid(self):
        non_valid = [
            self.d < 0,
            self.th0 < 0,
            self.th0 > self.d
        ]

        if any(non_valid):
            return False
        else:
            return True

    def set_sigma_mu(self, sigma=None, mu=None):
        self.sigma = sigma or self.sigma
        self.mu = mu or self.mu

    def __repr__(self):
        return 'Tube(id={}, d={}, th={}, th0={}, sigma={}, mu={})'.format(
            repr(self.id),
            repr(self.d),
            repr(self.th),
            repr(self.th0),
            repr(self.sigma),
            repr(self.mu),
        )

    def __str__(self):
        return 'Tube(id={:>30},\td={:>7},\tth={:>7},\tth0={:>7},\tsigma={:>11},\tmu={:>7})'.format(
            self.id,
            np.round(self.d, decimals=3),
            np.round(self.th, decimals=3),
            np.round(self.th0, decimals=3),
            np.round(self.sigma, decimals=3),
            np.round(self.mu, decimals=3),
        )

    def to_dict(self) -> dict:
        return {key: val for key, val in self.__dict__.items() if key in self.__scope__}

    def from_dict(self, init_dict: dict):
        [self.__setattr__(key, val) for key, val in init_dict.items() if key in self.__scope__]
        return self


pass
