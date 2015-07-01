#!/usr/bin/env python

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
    f=open("../ENGINE/"+add_name,'r')
    lines=f.readlines()
    f.close()
    f=open("../ENGINE/"+add_name,'w')
    sline0=lines[0].split()
    print >>f,"%s          %s    %s  %s    %s       %.3f %s    %s %s  %s %s"%(
               sline0[0],sline0[1],sline0[2],sline0[3],sline0[4],float(sline0[5])+0.001,sline0[6],sline0[7],sline0[8],sline0[9],sline0[10])
    for line in lines[1:]:
        print >>f,line[:-1]
    f.close()

def run_engine_timeout(engine_path,coengine_name,timeout):
    os.chdir("../ENGINE")
    status=RunCmd([engine_path+"engine.x","../ENGINE/"+coengine_name],timeout).Run()
    os.chdir("../ParFit")
    return status

def restore_hang_prm():
    shutil.copy("../ENGINE/add_MM3_hang.prm_orig","../ENGINE/add_MM3_hang.prm")

if __name__=="__main__":
    engine_path="/home/fzahari/Work/CMI/CODES/MM_engine_source/"
    os.environ["ENGINE_DIR"]=engine_path+"engine_dir"
    coeng_name="coengine_hang"
    timeout=10
    restore_hang_prm()
    while True:
        status=run_engine_timeout(engine_path,coeng_name,timeout)
        print "inside 'while':",status
        if status: 
            break
        else:
            pert_add_param("add_MM3_hang.prm")
    status=run_engine_timeout(engine_path,coeng_name,timeout)
    print status
    restore_hang_prm()
