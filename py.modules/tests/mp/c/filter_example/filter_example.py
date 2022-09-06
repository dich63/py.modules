from ...common import depth_units as units
import numpy as np
from ...constants import WIN_SIZE_M
from dataclasses import dataclass
from functools import partial
from scipy.signal import medfilt


@dataclass
class FilterResult:
    result: np.ndarray
    window_size: int


class FilterExample:
    def __init__(self, depth_step: float, depth_units: units.DepthUnits, window_size_m=None):
        window_size_m = window_size_m or WIN_SIZE_M

        window_size_points = int(window_size_m / depth_units.si_format(depth_step))
        window_size_points = np.clip(window_size_points, 1, None)
        window_size_points += 1 - window_size_points % 2
        self.window_size_points = int(window_size_points)

    def _get_trend(self, x: np.ndarray) -> np.ndarray:
        hs = self.window_size_points // 2
        result = np.zeros(x.shape[0] + self.window_size_points - 1)
        result[:hs], result[hs:-hs], result[-hs:] = x[0], x, x[-1]
        filter_ = partial(medfilt, kernel_size=self.window_size_points)
        return filter_(result)[hs:-hs]

    def process(self, x: np.ndarray) -> FilterResult:
        y = x.copy()
        for i, xi in enumerate(x.T):
            y[:, i] += -self._get_trend(y[:, i]) + np.median(y[:, i])
        return FilterResult(y, self.window_size_points)


pass
