#DivO DP module
#Version: 0.1.2
#Autors: Maciej Pietrzak, Michal Seweryn, Grzegorz Rempala
#Maintainer: Maciej Pietrzak <pietrzak.20@osu.edu>
#License: GPL (>=2)


import numpy as np
import math, sys, csv, os

def save_csv(f_out, listToSave):
    with open (f_out, "wt") as f: csv.writer(f).writerows(listToSave)
    return

def saveBootstrap(bs, arrays_for_analysis):
    bs_dir=bs
    if not os.path.exists(bs_dir):os.mkdir(bs_dir)
    n=1
    os.chdir(bs_dir)
    for i in arrays_for_analysis:
        save_csv(str(bs_dir)+"_"+str(n)+".csv", i)
        n+=1
    os.chdir("..")
    return

def draw_arrays(dat, resample, size):
    lista=dat.astype(float)
    lista=lista.reshape(1, len(lista[:,1])*len(lista[1,:]))[0]
    sum_dat=sum(lista)*float(size)
    x= np.random.multinomial(sum_dat, lista/lista.sum(), size=resample)
    arrays_for_analysis=[]
    for z in x:
        z=z.reshape(len(dat[:,1]), len(dat[1,:]))
        arrays_for_analysis.append(z)
    return arrays_for_analysis

def SH(x): 
    x[x==0] = 0.00000000000000000001
    return -sum((np.log(x))*x)

def RE(freq,ENS,alpha):
    freq[freq==0] = 0.00000000000000000001
    if alpha==1: re=SH(freq)
    else: re=np.log(sum(freq**alpha))/(1-alpha)
    if ENS==True: re=math.exp(re)
    return re

def HT(x, ENS, alpha, n):
    x=x.T   
    if alpha==1: alpha=1-0.000000001
    x=x[np.nonzero(x)]
    kisiel=alpha
    ht_out=(math.log(sum((x**alpha) / (1-(1-x)**n)))/(1-alpha))-(math.log(sum(x/(1-((1-x)**n))))/(1-alpha))
    if ENS==True: ht_out=math.exp(ht_out)
    return ht_out

def CVG(items):   
    items=list(items)
    return 1-(float(items.count(1))/sum(items))

def single_profile(repetition, alphas, ENS):
    return [RE(repetition,ENS,a) for a in alphas]

def single_profile_HT(repetition, alphas, ENS, n):
    return [HT(repetition,ENS,a, n) for a in alphas]

def CI_list_row(row, conf):
    row_mean= np.mean(row)
    row.sort()
    if (int(len(row)*conf))>0: row = row[(int(len(row)*conf)):-(int(len(row)*conf))]
    return [row_mean, row_mean-row[0], row[-1]-row_mean]

def confidence_interval_for_list(row_repetitions,conf): 
    conf=(1-conf)/2
    row_tmp=[CI_list_row(i, conf) for i in row_repetitions]
    
    return row_tmp

def DP_single_population(dat_single, alphas, ENS):
    col_repetitions= np.array([single_profile(dat_single[i], alphas, ENS) for i in range(len(dat_single))])

    return confidence_interval_for_list(col_repetitions.T,conf)

def DP_single_population_CVG(dat_single, alphas, cvg_n, ENS):
    col_repetitions= np.array([single_profile(dat_single[i], alphas*cvg_n[i], ENS) for i in range(len(dat_single))])

    return confidence_interval_for_list(col_repetitions.T,conf)

def DP_single_population_HT(dat_single, alphas, ENS):
    ns= np.sum(dat_single, axis=1)
    dat_single=dat_single.T
    ps=dat_single/np.sum(dat_single, axis=0)
    col_repetitions= np.array([single_profile_HT(ps[:,i], alphas, ENS, ns[i]) for i in range(len(ns))])

    return confidence_interval_for_list(col_repetitions.T,conf)

def DP_single_population_CVG_HT(dat_single, alphas, cvg_n, ENS):   
    ns= np.sum(dat_single, axis=1)
    dat_single=dat_single.T
    ps=dat_single/np.sum(dat_single, axis=0)
    col_repetitions= np.array([single_profile_HT(ps[:,i], alphas*cvg_n[i], ENS, ns[i]) for i in range(len(ns))])
    return confidence_interval_for_list(col_repetitions.T,conf)

def DP(dat, alpha_prof, arrays_for_analysis, conf,  cvg, ENS=False):
    if alpha_prof=='Def': alphas= np.arange(0.1, 2.1, 0.1, dtype=None)
    else: alphas=np.load('tmp_alpha.pyc')
    listCvg=[]
    if cvg == 'TRUE':
        for x in arrays_for_analysis:
            listCvg.append([CVG(x[:,i]) for i in range(len(x[1,:]))])
        arrays_for_analysis=[a.astype(float) for a in arrays_for_analysis]
        arrays_for_analysis_p=[a/np.sum(a, axis=0) for a in arrays_for_analysis]
        cos=[]
        cvg_out=[]
        for a in range(len(np.array(arrays_for_analysis_p).T)):
            i=np.array(arrays_for_analysis_p).T[a]
            i= i[~np.all(i == 0, axis=1)]
            cos.append( DP_single_population_CVG(i.T, alphas, np.array(listCvg).T[a], ENS))
            cvg_out.append(np.mean(np.array(listCvg).T[a]))
        np.save('DP_tmp_CVG', cvg_out)
        for i in range(len(cos)):
            cos[i] = np.append(np.array([alphas]).T, cos[i], 1)
            np.save('DP_tmp_output_'+str("%03d" % (i+1,)), cos[i])
 
    else:
        arrays_for_analysis=[a.astype(float) for a in arrays_for_analysis]
        arrays_for_analysis_p=[a/np.sum(a, axis=0) for a in arrays_for_analysis]
        cos=[]
        cvg_out=[]
        for a in range(len(np.array(arrays_for_analysis_p).T)):
            i=np.array(arrays_for_analysis_p).T[a]
            i= i[~np.all(i == 0, axis=1)]
            cos.append( DP_single_population(i.T, alphas, ENS))
        for i in range(len(cos)):
            cos[i] = np.append(np.array([alphas]).T, cos[i], 1)
            np.save('DP_tmp_output_'+str("%03d" % (i+1,)), cos[i])
    return

def DP_HT(dat, alpha_prof, arrays_for_analysis, conf,  cvg, ENS=False): 
    if alpha_prof=='Def': alphas= np.arange(0.1, 2.1, 0.1, dtype=None)
    else: alphas=np.load('tmp_alpha.pyc')
    listCvg=[]
    if cvg == 'TRUE':
        for x in arrays_for_analysis:
            listCvg.append([CVG(x[:,i]) for i in range(len(x[1,:]))])
        arrays_for_analysis=[a.astype(float) for a in arrays_for_analysis]
        array_of_ns=[np.sum(a, axis=0) for a in arrays_for_analysis]
        cos=[]
        cvg_out=[]
        for a in range(len(np.array(arrays_for_analysis).T)):
            i=np.array(arrays_for_analysis).T[a]
            i= i[~np.all(i == 0, axis=1)] 
            cos.append( DP_single_population_CVG_HT(i.T, alphas, np.array(listCvg).T[a], ENS))
            cvg_out.append(np.mean(np.array(listCvg).T[a]))
        np.save('DP_tmp_CVG', cvg_out)
        for i in range(len(cos)):
            cos[i] = np.append(np.array([alphas]).T, cos[i], 1)
            np.save('DP_tmp_output_'+("%03d" % (i+1,)), cos[i])
    else:
        arrays_for_analysis=[a.astype(float) for a in arrays_for_analysis]
        array_of_ns=[np.sum(a, axis=0) for a in arrays_for_analysis]
        cos=[]
        cvg_out=[]
        for a in range(len(np.array(arrays_for_analysis).T)):
            i=np.array(arrays_for_analysis).T[a]
            i= i[~np.all(i == 0, axis=1)]           
            cos.append(DP_single_population_HT(i.T, alphas, ENS))
        for i in range(len(cos)):
            cos[i] = np.append(np.array([alphas]).T, cos[i], 1)
            np.save('DP_tmp_output_'+("%03d" % (i+1,)), cos[i])
    return

test = sys.argv[1]
alpha= sys.argv[2]
resample=int(sys.argv[3])
conf = float(sys.argv[4])
f=(sys.argv[5])
PlugIn=(sys.argv[6])
size=sys.argv[7]
alpha_prof = sys.argv[8]
cvg= sys.argv[9]
bs= sys.argv[10]

if __name__ == '__main__':

    dat=np.load('tmp_in.pyc') 
    if len(dat.shape) ==1: dat=np.array([dat, dat]).T
    if dat.shape[1]==1: dat=np.concatenate([dat, dat], axis=1)
    if PlugIn=='TRUE': arrays_for_analysis=[dat]
    else: arrays_for_analysis=draw_arrays(dat, resample, size)
    if bs=='TRUE':
        bs_dir = "Bootstrap"
        saveBootstrap(bs_dir, arrays_for_analysis)
    if bs!="TRUE" and bs!="FALSE": saveBootstrap(bs, arrays_for_analysis)
    if test == 'DP'  and f == 'RE': DP(dat, alpha_prof, arrays_for_analysis, conf, cvg)                         
    if test == 'ENS' and f == 'RE': DP(dat, alpha_prof, arrays_for_analysis, conf, cvg, ENS=True)              
    if test == 'DP'  and f == 'HT': DP_HT(dat, alpha_prof, arrays_for_analysis, conf, cvg)          
    if test == 'ENS' and f == 'HT': DP_HT(dat, alpha_prof, arrays_for_analysis, conf, cvg, ENS=True) 




