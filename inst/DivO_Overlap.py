#DivO Overlap module
#Version: 0.1
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

def CVG(items):
    items=list(items)
    return 1-(float(items.count(1))/sum(items))

def CVG_IN(x):
    x=x.reshape(1, len(x[:,1])*len(x[1,:]))
    n= x[0]
    n=list(n)
    return 1-(float(n.count(1))/sum(n))

def I_Index_pooled(numbers, alpha):
    numbers= numbers[~np.all(numbers == 0, axis=1)] 
    numbers=numbers.astype(float)
    sum_global=np.sum(numbers) 
    cols= numbers.sum(axis=0)/float(sum_global) 
    rows= numbers.sum(axis=1)/float(sum_global) 
    table1= numbers.T.ravel()/sum_global
    if alpha==1: i_ind=1-((SH(rows)+SH(cols)- SH(table1))/SH(cols))
    else:        
        listForTab2=[]
        for d in range(len(cols)): listForTab2.append(np.multiply(rows, cols[d]))
        table2=np.concatenate(listForTab2)
        Fa=(((1/(alpha-1))*math.log(np.sum(np.multiply((table1**alpha), (table2**(1-alpha))))))) 
        p1=Fa/((1/(alpha-1))*math.log(sum((cols)**(2-alpha)))) 
        i_ind=1-p1 
    return i_ind

def I_Index(numbers, alpha):
    numbers= numbers[~np.all(numbers == 0, axis=1)]   
    numbers=numbers.astype(float)
    sum_global=np.sum(numbers) 
    cols= numbers.sum(axis=0)/float(sum_global) 
    rows= numbers.sum(axis=1)/float(sum_global) 
    table1= numbers.T.ravel()/sum_global
    if alpha==1: i_ind=1-((SH(rows)+SH(cols)- SH(table1))/SH(cols))
    else:        
        table2=np.concatenate((np.multiply(rows, cols[0]), np.multiply(rows, cols[1])))  
        Fa=(((1/(alpha-1))*math.log(np.sum(np.multiply((table1**alpha), (table2**(1-alpha))))))) 
        p1=Fa/((1/(alpha-1))*math.log(sum((cols)**(2-alpha)))) 
        i_ind=1-p1 
    return i_ind

def MH(numbers): 
    numbers= numbers[~np.all(numbers == 0, axis=1)]
    numbers=numbers.astype(float)
    smsq= (np.sum(numbers, axis=1))**2
    sq=numbers**2
    num= 2*sum(numbers[0,:]*numbers[1,:])*sum(numbers[0,:])*sum(numbers[1,:])
    den= np.sum(sq[0,:].T)*smsq[1]+np.sum(sq[1,:].T)*smsq[0]
    return num/den
    
def PG(alpha, beta, x): 
    x= x.astype(float)
    x= x[~np.all(x == 0, axis=1)] 
    x= x.T
    x= np.array(map(lambda a: a/np.sum(a, axis=0), x))
    x= x.T
    return  (sum((x[:,0]**alpha)*(x[:,1]**beta))+sum((x[:,1]**alpha)*(x[:,0]**beta)))/(sum((x[:,0]**(beta+alpha)))+ sum((x[:,1]**(alpha+beta))))

def PG_ht(alpha, beta, x): 
    x= x.astype(float)
    nlist=np.sum(x, axis=0)
    dpx=sum(x[:,0][x[:,0]==1])
    dpy=sum(x[:,1][x[:,1]==1])
    if dpx == nlist[0]: dpx= nlist[0]-1
    if dpy == nlist[0]: dpy= nlist[0]-1
    dpCx= 1-np.array([dpx,dpy])/nlist
    paxy=(x/nlist)*dpCx 
    abc= 1-((1-paxy)**(nlist)) 
    abc[abc==0]=1
    numerator1= sum(((paxy[:,0]**alpha)*(paxy[:,1]**beta))/(abc[:,0]*abc[:,1]))
    denominator=sum((paxy[:,0]**(2*alpha))/abc[:,0])+sum((paxy[:,1]**(2*beta))/abc[:,1])
    return (2*numerator1)/denominator
   
def LI(numbers): 
    numbers=numbers.astype(float)
    return  float(2*len(numbers[~(numbers==0).any(1)]))/len(numbers[np.nonzero(numbers)])

def JI(x): 
    x[x!=0] = 1
    a=x.T
    val_11=((a[:,0]==1) & (a[:,1]==1)).sum()
    val_01=((a[:,0]==0) & (a[:,1]==1)).sum()
    val_10=((a[:,0]==1) & (a[:,1]==0)).sum()
    return float(val_11)/(val_11+val_01+val_10)

def RDS(x, alpha):
    x= x.astype(float)
    freqs=np.asarray(map(lambda a: a/np.sum(a, axis=0), x))
    if alpha !=1: rd_out= 0.5*(((1/(alpha-1))*math.log(sum((freqs[0,:]**alpha)*(freqs[1,:]**(1-alpha))))) + ((1/(alpha-1))*math.log(sum((freqs[1,:]**alpha)*(freqs[0,:]**(1-alpha))))))
    else: rd_out=0 
    return rd_out

def RD(x, alpha):
    x= x.astype(float)
    #x=x.T
    #x=x[np.all(x!=0, axis=1)] #remove zeros for alpha >=1
    #x=x.T
    freqs=np.asarray(map(lambda a: a/np.sum(a, axis=0), x))
    if alpha !=1: rd_out=((1/(alpha-1))*math.log(sum((freqs[0,:]**alpha)*(freqs[1,:]**(1-alpha)))))
    else: rd_out= sum(freqs[0,:] * np.log(freqs[0,:]/freqs[1,:]))
    return rd_out

def confidence_interval(all_Resamples, conf): 
    conf=(1-conf)/2
    np.save('tmp_out_mean',  np.mean(all_Resamples, axis=0))
    sums= [np.linalg.norm(i) for i in all_Resamples]
    sumsprim=sums[:]
    sumsprim.sort()
    if (int(len(sumsprim)*conf))>0: sumsprim = sumsprim[(int(len(sumsprim)*conf)):-(int(len(sumsprim)*conf))]
    for position, item in enumerate(sums):
        if item == sumsprim[0]: np.save('tmp_out_min', all_Resamples[position])
        if item == sumsprim[-1]: np.save('tmp_out_max', all_Resamples[position])           
    return

def IN(dat, f, alpha, beta, arrays_for_analysis, conf, cvg): 
    all_Resamples=[]
    all_CVG=[]
    for r in range(len(arrays_for_analysis)):
        i_index_array=[]
        list_of_columns=[]
        numbers=arrays_for_analysis[r]
        if cvg=='TRUE':
            cvg_nums=map(lambda i: CVG(numbers[:,i]), range(len(numbers[1,:])))
            all_CVG.append(cvg_nums)
            for i in range(len(numbers[1,:])):               
                for j in range(len(numbers[1,:])):
                    if [j,i] not in list_of_columns and i != j:
                        if f == 'PG':
                            alpha = cvg_nums[i]
                            beta = cvg_nums[j]
                        if f =='PG_HT': 
                            alpha = cvg_nums[i]
                            beta = cvg_nums[j]
                        else: alpha = (cvg_nums[i]+cvg_nums[j])/2
                        if f=='INP': i_index_array.append(I_Index(np.array([numbers[:,i],numbers[:,j]]).T, alpha))
                        if f=='PG': i_index_array.append(PG(cvg_nums[i], cvg_nums[j], np.array([numbers[:,i],numbers[:,j]]).T))
                        if f=='PG_HT': i_index_array.append(PG_ht(alpha, beta, np.array([numbers[:,i],numbers[:,j]]).T))
                        if f=='MH': i_index_array.append(MH(np.array([numbers[:,i],numbers[:,j]])))
                        if f=='LI': i_index_array.append(LI(np.array([numbers[:,i],numbers[:,j]])))
                        if f=='RD': i_index_array.append(RD(np.array([numbers[:,i],numbers[:,j]]), alpha))
                        if f=='RDS': i_index_array.append(RDS(np.array([numbers[:,i],numbers[:,j]]), alpha))
                        if f=='JI': i_index_array.append(JI(np.array([numbers[:,i],numbers[:,j]])))
                        list_of_columns.append([i, j])
                    else: i_index_array.append(0)
        else:                
            for i in range(len(numbers[1,:])): 
                for j in range(len(numbers[1,:])):
                    if [j,i] not in list_of_columns and i != j:
                        if f=='INP': i_index_array.append(I_Index(np.array([numbers[:,i],numbers[:,j]]).T, alpha))
                        if f=='PG': i_index_array.append(PG(alpha, beta, np.array([numbers[:,i],numbers[:,j]]).T))
                        if f=='PG_HT': i_index_array.append(PG_ht(alpha, beta, np.array([numbers[:,i],numbers[:,j]]).T))
                        if f=='MH': i_index_array.append(MH(np.array([numbers[:,i],numbers[:,j]])))
                        if f=='LI': i_index_array.append(LI(np.array([numbers[:,i],numbers[:,j]]).T))
                        if f=='RD': i_index_array.append(RD(np.array([numbers[:,i],numbers[:,j]]), alpha))
                        if f=='RDS': i_index_array.append(RDS(np.array([numbers[:,i],numbers[:,j]]), alpha))
                        if f=='JI': i_index_array.append(JI(np.array([numbers[:,i],numbers[:,j]])))
                        list_of_columns.append([i, j])
                    else: i_index_array.append(0)               
        all_Resamples.append(np.reshape(np.array(i_index_array),(len(numbers[1,:]), len(numbers[1,:]))))

    if cvg=='TRUE': np.save('tmp_out_CVG',np.mean(all_CVG,axis=0))
    confidence_interval(all_Resamples, conf)
    return

test = sys.argv[1]
alpha= float(sys.argv[2])
resample=int(sys.argv[3])
conf = float(sys.argv[4])
faaa=(sys.argv[5]) 
PlugIn=(sys.argv[6])
size=sys.argv[7]
beta=sys.argv[8]
cvg= sys.argv[9]
bs= sys.argv[10]

if __name__ == '__main__':

    dat=np.load('tmp_in.pyc')
    if PlugIn=='TRUE': arrays_for_analysis=[dat]
    else: arrays_for_analysis=draw_arrays(dat, resample, size)
    if bs=='TRUE':
        bs_dir = "Bootstrap"
        saveBootstrap(bs_dir, arrays_for_analysis)
    if bs!="TRUE" and bs!="FALSE": saveBootstrap(bs, arrays_for_analysis)                                                                               
    if test=='IN': 
        conf=(1-conf)/2
        list_IN=[]
        for r in range(len(arrays_for_analysis)):
            if cvg == "TRUE": alpha= CVG_IN(arrays_for_analysis[r])
            list_IN.append(I_Index_pooled(np.array(arrays_for_analysis[r]), alpha))
        IN_mean= np.mean(list_IN)
        list_IN.sort()
        if (int(len(list_IN)*conf))>0: list_IN = list_IN[(int(len(list_IN)*conf)):-(int(len(list_IN)*conf))]
        np.save('tmp_IN',[[IN_mean, list_IN[0], list_IN[-1]]])
    if test=='INP': IN(dat, 'INP', float(alpha), 'beta', arrays_for_analysis, conf, cvg)
    if test=='MH': IN(dat, 'MH', 'alpha', 'beta', arrays_for_analysis, conf, cvg)  
    if test=='PG': IN(dat, 'PG', float(alpha), float(beta), arrays_for_analysis, conf, cvg)  
    if test=='PG_HT': IN(dat, 'PG_HT', float(alpha), float(beta), arrays_for_analysis, conf, cvg)
    if test=='RD': IN(dat, 'RD', alpha, 'beta', arrays_for_analysis, conf, cvg)                 
    if test=='RDS': IN(dat, 'RDS', alpha, 'beta', arrays_for_analysis, conf, cvg)
    if test=='LI': IN(dat, 'LI', 'alpha', 'beta', arrays_for_analysis, conf,cvg)  
    if test=='JI': IN(dat, 'JI', 'alpha', 'beta', arrays_for_analysis, conf, cvg)                



