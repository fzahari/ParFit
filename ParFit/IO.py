#!/usr/bin/env python

from numpy import zeros

def ginp_inp(input_fname):
    f=open(input_fname,'r')
    lines=f.readlines()
    f.close()
    line1=lines[0][:-1].split(',')
    if len(line1)==3:
       gopt_scan_fnameb,tup,bes=line1
       scan_type="diha"
    else:
       gopt_scan_fnameb,tup,bes,scan_type=line1
    gopt_scan_fnameb=gopt_scan_fnameb.strip()
    if scan_type=="diha":
       t1,t2,t3,t4=map(int,tup.split())
       tup=(t1,t2,t3,t4)
    elif scan_type=="bond":
       t1,t2=map(int,tup.split())    
       tup=(t1,t2)
    elif scan_type=="angl":
       t1,t2,t3=map(int,tup.split())    
       tup=(t1,t2,t3)
    b,e,s=map(int,bes.split())
    bes=(b,e,s)
    return scan_type,gopt_scan_fnameb,tup,bes

def par_fit_inp(input_fname):
    f=open(input_fname,'r')
    lines=f.readlines()
    f.close()
    line1=lines[0][:-1].split(',')
    if line1[0].strip()=="mult": 
       n=int(line1[1])
       gopt_type=[]
       gopt_scan_fnameb=[]
       tup=[]
       bes=[]
       if len(line1)==2:
          scan_type="diha"      
       else:
          scan_type=line1[2].strip()
       for i in range(n):
          sline=lines[i+1].split(',')
          if sline[0].strip()=="comp":
             gopt_type.append(sgopt_type.strip())
             gopt_scan_fnameb.append(sgopt_scan_fnameb.strip())
             tup.append(())
             bes.append(())
          else:
             gopt_type.append(sgopt_type.strip())
             gopt_scan_fnameb.append(sgopt_scan_fnameb.strip())
             if scan_type=="diha":
                t1,t2,t3,t4=map(int,stup.split())    
                tup.append((t1,t2,t3,t4))
             elif scan_type=="bond":
                t1,t2=map(int,stup.split())    
                tup.append((t1,t2))
             elif scan_type=="angl":
                t1,t2,t3=map(int,stup.split())    
                tup.append((t1,t2,t3))
             b,e,s=map(int,sbes.split())    
             bes.append((b,e,s))
    else:
       n=0
       if line1[0].strip()=="comp":
          if len(line1)==2:
             gopt_type,gopt_scan_fnameb=line1
             scan_type="diha"
          else:
             gopt_type,gopt_scan_fnameb,scan_type=line1
          gopt_type=[gopt_type.strip(),]
          gopt_scan_fnameb=[gopt_scan_fnameb.strip(),]
          tup=[(),]
          bes=[(),]
       else:
          if len(sline)==4:
             gopt_type,gopt_scan_fnameb,tup,bes=line1
             scan_type="diha"
          else:
             gopt_type,gopt_scan_fnameb,tup,bes,scan_type=line1
          gopt_type=[gopt_type.strip(),]
          gopt_scan_fnameb=[gopt_scan_fnameb.strip(),]
          if scan_type=="diha":
             t1,t2,t3,t4=map(int,tup.split())
             tup=[(t1,t2,t3,t4),]
          elif scan_type=="bond":
             t1,t2=map(int,tup.split())    
             tup=[(t1,t2),]
          elif scan_type=="angl":
             t1,t2,t3=map(int,tup.split())    
             tup=[(t1,t2,t3),]
          b,e,s=map(int,bes.split())
          bes=[(b,e,s),]
    if gopt_type[0]=="ginp":
       return scan_type,gopt_type,gopt_scan_fnameb,tup,bes,None,None,None,None,None,None,None,None,None
    engine_path=lines[n+1][:-1]
    l2s=lines[n+2].split()
    mm=l2s[0].lower()
    if len(l2s)==2:
       mode=l2s[1].lower()
    else:
       mode="opt"
    np=0
    #plist=[]
    nc=0
    opt_lin={}
    alg=lines[n+3].strip()
    for line in lines[n+4:-1]:
	t=line.split()
        n=len(t)-1
        for i in range(n):
            t[i+1]=t[i+1].lower()
            if t[i+1][0]=='p': 
               #plab=t[i+1][1:]
               #if len(plab)>0:               
               #   if plab in plist:
               #      continue
               #   else:
               #      plist.append(plab)
               np+=1
            if t[i+1]=='c': 
               nc+=1
        if n==3:
           opt_lin.update({int(t[0])-1:(t[1],t[2],t[3])})
        elif n==2:
           opt_lin.update({int(t[0])-1:(t[1],t[2])})
        elif n==1:
           opt_lin.update({int(t[0])-1:(t[1],)})
    last_line=lines[-1].split()
    if len(last_line)==1:
       step_int=10
       csv=lines[-1].strip()
    elif len(last_line)==2:
       step_int,csv=lines[-1].split()
    else:
       print "Wrong csv line in the input file!"

    return scan_type,gopt_type,gopt_scan_fnameb,tup,bes,engine_path,mm,mode,alg,opt_lin,np,nc,step_int,csv

def read_add(mm,opt_lin,np,nc,fl,scan_type):

    p=zeros(np)
    c=zeros(nc)
     
    if mm=="mm3":
        if fl==1:
           add_name="../Data/Engine/add_MM3.prm"
        else:
           add_name="../Data/Engine/add_MM3_PF.prm"
    elif mm=="mmff94":
        if fl==1:
           add_name="../Data/Engine/add_MMFF94.prm"
        else:
           add_name="../Data/Engine/add_MMFF94_PF.prm"
    else:
        print "read_add: Wrong MM-type!"

    f=open(add_name,'r')
    lines=f.readlines()
    f.close()
        
    i=j=0
    ol_keys=opt_lin.keys()
    ol_keys.sort()
    pdict={}
    ol_templ={}
    scan_type=scan_type.strip()
    if scan_type=="diha":
       n=3
    elif scan_type=="bond":
       n=2
    elif scan_type=="angl":
       n=2
    for ol_sk in ol_keys:
        s=""
        t=opt_lin[ol_sk]
        for k in range(n):
            if t[k][0]=='p':
                plab=t[k][1:]
                if len(plab)>0:
                   if pdict.has_key(plab):
                      ii=pdict[plab]
                   else:
                      pdict.update({plab:i})
                      ii=i
                      i+=1
                else:
                   ii=i
                   i+=1
                s+="p["+str(ii)+"],"
                #s+="p["+str(i)+"],"
                #i+=1
            elif t[k]=='c':
                s+="c["+str(j)+"],"
                j+=1
            else:
                print "Woops!" 
        v=lines[ol_sk].split()
        ol_vars=s[:-1]
        if scan_type=="diha":
           exec ol_vars+"="+v[5]+","+v[7]+","+v[9]
           v[1:5]=map(int,v[1:5])
           s1='"%s        %3i  %3i  %3i  %3i'%(v[0],v[1],v[2],v[3],v[4])
           s2='      %6.3f +1   %6.3f -2 %6.3f +3      "%('+ol_vars+')'
        elif scan_type=="bond":
           exec ol_vars+"="+v[3]+","+v[4]
           v[1:3]=map(int,v[1:3])
           s1='"%s        %3i  %3i'%(v[0],v[1],v[2])
           s2='      %6.3f    %6.3f        "%('+ol_vars+')'
        elif scan_type=="angl":
           exec ol_vars+"="+v[4]+","+v[5]
           v[1:4]=map(int,v[1:4])
           s1='"%s        %3i  %3i  %3i'%(v[0],v[1],v[2],v[3])
           s2='      %6.3f    %6.3f        "%('+ol_vars+')'
        ol_templ.update({ol_sk:s1+s2})

    return p,c,ol_templ,lines

def write_add(sdir,p,c,mm,ol_templ,lines,fl,step,step_int):
    if mm=="mm3":
       if fl==1:
          add_name="../Data/Engine/add_MM3.prm"
       else:
          add_name="../Data/Engine/add_MM3_PF.prm"
    elif mm=="mmff94":
       if fl==1:
          add_name="../Data/Engine/add_MMFF94.prm"
       else:
          add_name="../Data/Engine/add_MMFF94_PF.prm"
    else:
       print "write_add: Wrong MM-type!"

    f=open(add_name,'w')

    ol_keys=ol_templ.keys()
    ol_keys.sort()
    for ol_sk in ol_keys:
        exec 'lines['+str(ol_sk)+']='+ol_templ[ol_sk]

    for line in lines:
       print >>f,line[:-1]

    f.close()
   
    if step==None: return
    
    if step=="ga" or (step-1)%step_int==0:
       if mm=="mm3":
          add_name_arch=sdir[0]+"/add_MM3_"+str(step)+".prm"
       elif mm=="mmff94":
          add_name_arch=sdir[0]+"/add_MMFF94_"+str(step)+".prm"

       f=open(add_name_arch,'w')
       for line in lines:
          print >>f,line[:-1]
       f.close()

    return

if __name__=="__main__":
    gmol=DihGOpt_Molecule(mm="mmff94",name="opmmm")    
    print gmol.name
    #gmol.dt=(13,4,0,1)
    gmol.dt=(1,0,4,13)
    print gmol.dt
    fn_b="opmmm-mp2-opt"
    gmol.read_gopt_log(fn_b)
    print "tl:",gmol.tl
    print gmol.e
    print gmol.sl
    print gmol.rl
    print gmol.cl
    print gmol.da
    gmol.dih_rot(pi/2.0)
    print gmol.rl
    print gmol.da
    print gmol.write_inp_pcm("test")
    gmol=DihGOpt_Molecule(mm="mmff94")
    gmol.read_ginp("tempC2")    
    gmol.dt=(13,4,0,1)
    print gmol._ginp_templ
    print gmol._rl
    gmol.dih_rot(pi/2.0)
    print gmol._rl
