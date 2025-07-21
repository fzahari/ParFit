#!/usr/bin/env python3

import subprocess as sub
import threading
import os
import shutil

 
class RunCmd(threading.Thread):
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
        self.p = sub.Popen(self.cmd)
        self.p.wait()

    def Run(self):
        self.start()
        self.status=True
        self.join(self.timeout)
        if self.is_alive():
            self.p.terminate()
            self.join()
            self.status=False
        return self.status

def pert_add_param(add_name):
    f=open("../Data/Engine/"+add_name,'r')
    lines=f.readlines()
    f.close()
    f=open("../Data/Engine/"+add_name,'w')
    sline0=lines[0].split()
    print("%s          %s    %s  %s    %s       %.3f %s    %s %s  %s %s"%(
               sline0[0],sline0[1],sline0[2],sline0[3],sline0[4],float(sline0[5])+0.001,sline0[6],sline0[7],sline0[8],sline0[9],sline0[10]), file=f)
    for line in lines[1:]:
        print(line[:-1], file=f)
    f.close()

def run_engine_timeout(engine_path,coengine_name,timeout):
    os.chdir("../Data/Engine")
    status=RunCmd([engine_path+"/engine.x",coengine_name],timeout).Run()
    #os.system(engine_path+"/engine.x "+coengine_name)
    #status=True
    os.chdir("../../ParFit")
    return status

def restore_hang_prm():
    shutil.copy("../Data/Engine/add_MM3_hang.prm_orig","../Data/Engine/add_MM3_hang.prm")

if __name__=="__main__":
    engine_path="/Users/federicozahariev1/Work/Programs/CODES/MM_engine_source/"
    os.environ["Engine_DIR"]=engine_path+"engine_dir"
    coeng_name="coengine_hang"
    timeout=10
    restore_hang_prm()
    while True:
        status=run_engine_timeout(engine_path,coeng_name,timeout)
        print("inside 'while':", status)
        if status: 
            break
        else:
            pert_add_param("add_MM3_hang.prm")
    status=run_engine_timeout(engine_path,coeng_name,timeout)
    print(status)
    restore_hang_prm()
