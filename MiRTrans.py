import os
import sys
from collections import defaultdict
from sklearn import linear_model
import numpy as np
import scipy.stats
#import pdb
#pdb.set_trace()
import warnings
warnings.filterwarnings('ignore')
class parameter_struc(object):
    def __int__(self):
        self.Seq_Pair='N'
        self.Degrad_Pair='N'
        self.MirExp='N'
        self.MrnaExp='N'
def input_parameter(argv,parameter_struc):
    deter=1
    f=open(argv,"r")
    for line in f:
        Par=line.split('=')
        if len(Par)==2:
           if Par[0]=='Pred_seq':
              parameter_struc.Seq_Pair=Par[1].strip('\n')
           elif Par[0]=='Pred_deg':
              parameter_struc.Degrad_Pair=Par[1].strip('\n')
           elif Par[0]=='Exp_miRNA':
              parameter_struc.MirExp=Par[1].strip('\n')
           elif Par[0]=='Exp_mRNA':
              parameter_struc.MrnaExp=Par[1].strip('\n')
    f.close()
    if os.path.isfile(parameter_struc.Seq_Pair)=='N':
       deter=0
       print('The file with sequence based microRNA target prediction does not exist')
    if os.path.isfile(parameter_struc.Degrad_Pair)=='N':
       deter=0
       print('The file with degradome sequencing does not exist')
    if os.path.isfile(parameter_struc.MirExp)=='N':
       deter=0
       print('The file with microRNA expression does not exist')
    if os.path.isfile(parameter_struc.MrnaExp)=='N':
       deter=0
       print('The file with mRNA expression does not exist')
    return deter
class Candidate_pair(object):
    def __ini__(self,mirna,mrna):
        self.mirna=mirna
        self.mrna=mrna
def input_target_pair(in_path):
    pair_info=[]
    f=open(in_path,"r")
    line_index=0
    pair_dict=defaultdict(list)
    for line in f:
        if line_index>0:
           pair_full=line.strip('\n')
           pair_div=pair_full.split('\t')
           pair_dict[pair_div[0]].append(pair_div[1])
        line_index=line_index+1
    f.close()
    return pair_dict
def input_deg(in_path):
    pair_info=[]
    f=open(in_path,"r")
    line_index=0
    pair_dict=defaultdict(list)
    for line in f:
        if line_index>0:
           pair_full=line.strip('\n')
           pair_div=pair_full.split('\t')
           pair_dict[(pair_div[0],pair_div[1])].append(pair_div[2])
        line_index=line_index+1
    f.close()
    return pair_dict
def input_expression(in_path):
    dict_rna={}
    f=open(in_path,"r")
    line_index=0
    for line in f:
        if line_index>0:
           exp_full=line.strip('\n')
           exp_div=exp_full.split('\t')
           dict_rna[exp_div[0]]=list(map(float,exp_div[1:]))
        line_index=line_index+1
    f.close()
    return dict_rna
def cal_dep(Seq_dict,Deg_dict,miR_exp,mRNA_exp):
    f=open('MicroTrans_results.txt',"w")
    f.write('microRNA')
    f.write('\t')
    f.write('mRNA')
    f.write('\t')
    f.write('p-value')
    f.write('\n')
    #perform lasso regression
    lasso_reg={}
    for (Key_mirna,Value_mrna) in Seq_dict.items():
        print('processing microRNA '+Key_mirna,end="\r" )
        TotalmRNA=[]
        for i in range(len(Value_mrna)):
           TotalmRNA.append(mRNA_exp[Value_mrna[i]])
        clf = linear_model.LassoLarsIC(criterion='bic')
        clf.fit(np.transpose(np.asarray(TotalmRNA)),np.asarray(miR_exp[Key_mirna]))
        if len(np.nonzero(clf.coef_))==0:
           continue
        stdev=bootstrap(np.asarray(miR_exp[Key_mirna]),np.asarray(TotalmRNA),len(Value_mrna))
        for j in range(len(clf.coef_)):
            if clf.coef_[j]!=0:
               lasso_reg[(Key_mirna,Value_mrna[j])]=1-round(scipy.stats.norm(0, 1).cdf(clf.coef_[j]/stdev[j]),3)
    lasso_reg_set=set(lasso_reg)
    deg_set=set(Deg_dict)
    sharedKey={}
    for inter_key in lasso_reg_set.intersection(deg_set):
        sharedKey[(inter_key[0],inter_key[1])]=1
        Pvalue=scipy.stats.combine_pvalues([float(lasso_reg[inter_key]),float(deg_set[inter_key])], method='fisher')
        output(inter_key[0],inter_key[1],Pvalue,f)
    for uniq_key in lasso_reg.keys():
        if uniq_key not in sharedKey.keys():
           output(uniq_key[0],uniq_key[1],lasso_reg[uniq_key],f)
    for uniq_key in Deg_dict.keys():
        if uniq_key not in sharedKey.keys():
           output(uniq_key[0],uniq_key[1],Deg_dict[uniq_key][0],f)
    f.close()
    print('Succesfully finished, the results are in MicroTrans_results.txt')
    return None  
def output(mirna_name,mrna_name,p_value,f):
    f.write(mirna_name)
    f.write('\t')
    f.write(mrna_name)
    f.write('\t')
    f.write(str(round(float(p_value),4)))
    f.write('\n')
    return None
def bootstrap(miR_exp,mRNA_exp,num_mrna):
    clf = linear_model.LassoLarsIC(criterion='bic')
    boot_coll=[]
    std=[]
    for i in range(1000):
        rand_mrna=np.transpose(np.random.permutation(np.transpose(mRNA_exp)))
        clf.fit(np.transpose(rand_mrna),np.transpose(miR_exp))
        boot_coll.append(clf.coef_)
    stdev=np.std(boot_coll,axis=0)
    return stdev
def main():
    Par=parameter_struc()
    deter=input_parameter(sys.argv[1],Par)
    if deter==1:
        print('Analysis Begin...')
        #input sequence based prediction
        Seq_dict=input_target_pair(Par.Seq_Pair)
        print('Finish reading sequence based prediction results')
        #input degradome sequencing
        Deg_dict=input_deg(Par.Degrad_Pair)
        print('Finish reading degradome sequencing results')
        #input miRNA expression
        miR_exp=input_expression(Par.MirExp)
        print('Finish reading microRNA expression')
        #input mRNA expression
        mRNA_exp=input_expression(Par.MrnaExp)
        print('Finish reading mRNA expression')
        cal_dep(Seq_dict,Deg_dict,miR_exp,mRNA_exp)
    return None
if __name__=='__main__':
    main()
