import numpy as np

def smooth_gaussian(x,window_len=4,window='hamming'):
	if x.ndim != 1:
		raise ValueError("smooth_gaussian only accepts 1 dimension arrays.")
	if x.size < window_len:
		raise ValueError("Input vector needs to be bigger than window size.")
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError("Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

	s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')

	y=np.convolve(w/w.sum(),s,mode='valid')

	return y[int(window_len/2-1):-int(window_len/2)]