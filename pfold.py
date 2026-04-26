
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/tiger/Desktop/chandra/code")
from astropy.stats import poisson_conf_interval
import GL_algorithm.funcs as funcs
fsize=24
font1 = {'weight': 'normal',
         'size': 24, }
plt.rc('legend',fontsize=14 )

def trans(t, p_test, shift):
    ti = t
    v = 1.0 / p_test
    turns = v * ti
    turns += shift
    # 初始相位
    for i in range(len(turns)):
        turns[i] = turns[i] - int(turns[i])
        # 折叠后的相位
    return turns

def phase_fold(time, epoch_info, p_test, bin=20, net_percent=0.9, shift=0.0, label='test', save=False, sname=None,
               show=True, bkg=False, time_bkg=None):

    # 默认time和epoch_info已经filter好了
    if bkg:
        net_percent = 1 - len(time_bkg) / len(time) * 1 / (4 * 4 - 2 * 2)
    print('counts:', len(time))
    turns = trans(time, p_test, shift)
    loc = np.zeros(bin)
    for index in turns:
        loc[int(index * bin)] += 1
    AM = 1 - min(loc) / max(loc)
    A0 = AM / (2 - AM)
    print('A0={0}'.format(A0))
    x = np.array([(i / bin + 0.5 / bin) for i in range(bin)])
    src_bkg = 1 - net_percent
    bkg_y = len(time) * src_bkg / bin
    b_1sigma = poisson_conf_interval(bkg_y, interval='frequentist-confidence').T
    bkg_y_low = b_1sigma[0];
    bkg_y_high = b_1sigma[1]
    fig = plt.figure(2, (10, 7.5))
    ax1 = fig.add_subplot(111)
    bkg_x = [0, 2]
    plt.fill_between(bkg_x, bkg_y_low, bkg_y_high, facecolor='blue', alpha=0.5)
    x2 = np.concatenate((x, x + 1))
    y2 = np.concatenate((loc, loc))
    T_in_perbin = funcs.get_T_in_mbins(epoch_info, 2 * np.pi / p_test, bin, shift * 2 * np.pi)

    correct_gap = T_in_perbin / (sum(T_in_perbin) / len(T_in_perbin))
    print('correct_gap=', correct_gap)
    y2 /= np.concatenate((correct_gap, correct_gap))
    y2_err = np.array(poisson_conf_interval(y2, interval='frequentist-confidence'))
    y2_err[0] = y2 - y2_err[0]
    y2_err[1] = y2_err[1] - y2

    plt.title("#{0} P={1:.2f} s C={2}".format(label, p_test, str(len(time))), fontsize=fsize)
    plt.xlabel('phase', font1)
    plt.ylabel('counts/bin', font1)
    plt.tick_params(labelsize=fsize)
    plt.ylim(0, (np.max(y2) + np.max(y2) ** 0.5) * 1.05)
    plt.step(np.concatenate(([0], x2)), np.concatenate(([y2[0]], y2)), color='red')
    plt.errorbar(x2 - 0.5 / bin, y2, yerr=y2_err, fmt='.', capsize=1, elinewidth=1, ecolor='red')

    print('net_percent:{0:.5f}'.format(net_percent))
    ax2 = ax1.twinx()
    yhigh = (np.max(y2) + np.max(y2) ** 0.5) * 1.05 / np.mean(y2)
    ax2.set_ylabel('Normalized flux', font1)
    ax2.plot([0, 2], [1.0, 1.0], '--', color='blue')
    ax2.set_ylim([0, yhigh])
    ax2.tick_params(labelsize=fsize)
    plt.tight_layout()

    if save: plt.savefig(fname=sname,bbox_inches='tight',dpi=1200)
    if show: plt.show()
    plt.close()
    return

if __name__=='__main__':
    fp="/Users/tiger/Library/Containers/com.tencent.xinWeChat/Data/Library/Application Support/com.tencent.xinWeChat/2.0b4.0.9/15934c6431be0b4989ca325fd97c5481/Message/MessageTemp/519c6de839abf322809caa191734d909/File/Sim_time.npy"
    time=np.load(fp)
    time=np.sort(time)
    counts=len(time)
    epoch_info=np.array([[time[0],time[-1],0,time[-1]-time[0]]])
    phase_fold(time, epoch_info, p_test=316.55587211141227, bin=20, net_percent=0.9, shift=0.0, label='test', save=False, sname=None,
               show=True, bkg=False, time_bkg=None)