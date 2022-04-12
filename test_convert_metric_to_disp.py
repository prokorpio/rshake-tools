import STALTA_trigger as target
import numpy as np
from obspy.core import Trace, Stream


def test_sine_velocity():
    # create input
    t = np.linspace(0,2*np.pi*10,int(100*2*np.pi*10)+1) # ~1min data, 100Hz sampling
    tr = Trace(data=np.sin(t))
    tr.stats.delta = tr.data[1] - tr.data[0] # 0.01 for sampling_rate of 100
    tr.stats.units = "VEL"
    st = Stream(traces=[tr])

    # execute target fn
    res = target.convert_metric_to_disp(st)

    # check result
    correct_res = -np.cos(t) + 1
    error_bound = 0.001 # within 1 milli meter error is good enough
    max_error = max(abs(res[0].data - correct_res))
    assert max_error < error_bound,\
           "Should be -cos(t) + 1 +/- 1e-03, got max error of " + str(max_error)



