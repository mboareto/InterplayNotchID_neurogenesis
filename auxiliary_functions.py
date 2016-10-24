import numpy as np
import pylab as plt

def solve_eqs(dde, trange, dtmax=0.1, log2=True):
    dde.run()
    r = dde.sample(trange[0], trange[1], dtmax)
    if log2:
        for k in dde.vars:
            r[k] = np.log2(r[k]+1.0)
    return r
    
def oscillation_parameters(s, var, warning=False):
    i_max = np.r_[False, s[var][1:] > s[var][:-1]] & np.r_[s[var][:-1] > s[var][1:], False]
    i_min = np.r_[False, s[var][1:] < s[var][:-1]] & np.r_[s[var][:-1] < s[var][1:], False]
    if (sum(i_min==True) > 1) & (sum(i_max==True) > 1):
        T_max = s['t'][i_max][-1] - s['t'][i_max][-2]
        t_max = s['t'][i_max][-1]
        T_min = s['t'][i_min][-1] - s['t'][i_min][-2]
        v_mean= np.mean(s[var][(s['t'] > s['t'][i_max][-2]) & (s['t'] < s['t'][i_max][-1])])
        v_max = s[var][i_max][-1]
        v_min = s[var][i_min][-1]
    else:
        T_max = 0.0
        T_min = 0.0
        v_mean= np.mean(s[var])
        v_max = v_mean
        v_min = v_mean
        t_max = 0
        
    if warning:
        if T_max > 0.1:
            if abs(1.0 - T_min/T_max) > 0.1:
                print 'Warning: Diference in period T_max and T_min higher than 10%', T_min, T_max
                T_max = T_min = 0.0
    return {'T'    : 0.5*(T_min + T_max), 
            'amp'  : v_max - v_min, 
            'max'  : v_max, 
            'min'  : v_min,
            'mean' : v_mean,
            't_max': t_max}

def bifurcation(dde, par, r_par, var, trange, c=['-b','-b','--b'], label=None, yrange=None, dyrange=2.0, 
                xlim=None, log2=True, ax=None, fs=[6,5]):
    if ax == None:
        fig, ax = plt.subplots(figsize=(fs[0],fs[1]))

    p = {}
    p[par] = dde.params[par]
    out = np.zeros((len(r_par), 3))

    for i in range(len(r_par)):  
        dde.params[par] = r_par[i]
        s = solve_eqs(dde, trange=trange, dtmax=dde.params['dt'])
        out[i,:] = [oscillation_parameters(s, var)[k] for k in ['min','max','mean']]
    dde.params[par] = p[par]

    for i in range(3):
        if par=='gp':
            ax.plot(r_par*np.log(2), out[:,i], c[i], lw=1.5)
        elif log2:  
            ax.plot(np.log2(r_par+1),out[:,i], c[i], lw=1.5)
        else:
            ax.plot(r_par, out[:,i], c[i], lw=1.5)
        if label==None:
            plt.xlabel(par)
            plt.ylabel(var)
        else:
            plt.xlabel(label[par])
            plt.ylabel(label[var])
        if yrange!=None:
            plt.ylim(yrange)
            plt.yticks(np.arange(yrange[0], yrange[1]+1.0, dyrange))
        if xlim!=None:
            plt.xlim(xlim)

def phase_diagram(dde, p_dic, axis, var, trange, a_tr=0.01, keys=['T','amp','max','min','mean'], warning=False):
    p_init = {}
    for k in axis:
        p_init[k] = dde.params[k]
        
    out = {}
    for key in var:
        out[key] = {}
        for k in keys:
            out[key][k] = np.zeros((len(p_dic[axis[0]]), len(p_dic[axis[1]])))
    
    for i in range(len(p_dic[axis[0]])):
        dde.params[axis[0]] = p_dic[axis[0]][i]
        for j in range(len(p_dic[axis[1]])):
            dde.params[axis[1]] = p_dic[axis[1]][j]
            s = solve_eqs(dde, trange, dtmax=dde.params['dt'])
            for key in var:
                x = oscillation_parameters(s, key, warning=warning)
                for k in keys:
                    out[key][k][i,j] = x[k]
    for k in axis:
        dde.params[k] = p_init[k]
    return out   

def plot_pdiagram(m, p_dic, axis, label=None, k='mean', amp_tr=0.1, colorbar=True, clabel=None, 
                  clim=None, nc=2, fs=[6,4], cmap=None, ax=None):
    if ax == None:
        fig, ax = plt.subplots(figsize=(fs[0],fs[1]))
        
    m['T'  ][m['amp'] < amp_tr*m['mean']] = 0.0
    m['amp'][m['amp'] < amp_tr*m['mean']] = 0.0
    plt.pcolor(np.log2(p_dic[axis[0]]), np.log2(p_dic[axis[1]]), np.transpose(m[k]), cmap=cmap)
    if colorbar:
        clb = plt.colorbar(label=k)
        if clabel!=None:
            clb.set_label(clabel, rotation=90)
    if clim!=None:
        plt.clim(clim)
    if label!=None:
        plt.xlabel(label[axis[0]])
        plt.ylabel(label[axis[1]])
    plt.xlim([np.log2(np.min(p_dic[axis[0]])),np.log2(np.max(p_dic[axis[0]]))])
    plt.ylim([np.log2(np.min(p_dic[axis[1]])),np.log2(np.max(p_dic[axis[1]]))])
