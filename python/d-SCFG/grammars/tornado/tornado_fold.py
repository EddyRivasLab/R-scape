import os
from pathlib import Path
import subprocess

import re

import math
import numpy as np

file_path = Path("/home/user/documents/my_file.txt")
directory = file_path.parent


def tornado_fold(outdir, method, param_file, grm_file, postgrm_file, test_file):

    test_path = Path(test_file)
    test_dir  = str(test_path.parent)
    test_name = str(test_path.stem)
    
    param_path = Path(param_file)
    param_name = str(param_path.stem)

    test_name_noss = re.search(r'(\S+).noss', test_name)
    if (test_name_noss):
        test_sto  = test_dir+"/"+str(test_name_noss.group(1))+".sto"
    else:  
        test_sto  = test_dir+"/"+test_name+".sto"
 
    sen = 0
    ppv = 0
    f1  = 0
    if (os.path.exists(test_sto) == False):
        return sen, ppv, f1

    res_file   = str(outdir)+"/"+param_name+"."+test_name+"."+method+".sto"
    out_file   = res_file+".tornado"
    stats_file = res_file+".tornado.stats"
    #print(f"results file: {stats_file}")
   
    # grm-fold
    cmd = ""
    cmd += "../../bin/grm-fold --"+method
    cmd += " --lprob "
    cmd += " --gpostfile "+postgrm_file
    cmd += " "+grm_file
    cmd += " "+test_file
    cmd += " "+res_file
    cmd += " "+param_file
    cmd += " > "+out_file
    #print('mea_fold:', cmd)
    os.system(cmd)  

    # compare to given structure
    cmd  = ""
    cmd += "../../lib/hmmer/easel/miniapps/esl-compstruct"
    cmd += " "+test_sto
    cmd += " "+res_file
    cmd += " > "+stats_file
    #print('compstruct:', cmd)
    os.system(cmd)

    sen, ppv, f1 = grmfold_stats_parse(stats_file)
        
    return sen, ppv, f1


def grmfold_stats_parse(stats_file):
    sen = 0
    ppv = 0
    f1  = 0

    # Overall prediction accuracy (1415 sequences, 104851 positions)
    # 11326/29469 trusted pairs predicted (38.43% sensitivity)
    # 11326/47045 predicted pairs correct (24.07% PPV)

    sen_pattern = r'trusted pairs predicted \((\d+.\d+)\% sensitivity\)'
    ppv_pattern = r'predicted pairs correct \((\d+.\d+)\% PPV\)'
    
    fp = open(stats_file, "r")
    lines = fp.readlines()
    for line in lines:
        sen_match = re.search(sen_pattern, line)
        if sen_match:
            sen = float(sen_match.group(1))
        ppv_match = re.search(ppv_pattern, line)
        if ppv_match:
            ppv = float(ppv_match.group(1))
            
    if sen + ppv > 0:
        f1 = 2.0 * sen * ppv / (sen + ppv)
        
    return sen, ppv, f1
