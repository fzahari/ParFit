#!/usr/bin/env python3

from numpy import zeros

def ginp_inp(input_fname):
    f=open(input_fname,'r')
    lines=f.readlines()
    f.close()
    line1=lines[0][:-1].split(',')
    if line1[0].lower()=="nwchem":
       prog_type="nwchem"
       line1=line1[1:]
    else:
       prog_type="gamess"   
    if len(line1)==3:
       gopt_scan_fnameb,tup,bes=line1
       scan_type="diha"
    else:
       gopt_scan_fnameb,tup,bes,scan_type=line1
    gopt_scan_fnameb=gopt_scan_fnameb.strip()
    scan_type=scan_type.strip()
    if scan_type=="diha":
       t1,t2,t3,t4=map(int,tup.split())
       tup=(t1,t2,t3,t4)
    elif scan_type=="bond":
       t1,t2=map(int,tup.split())    
       tup=(t1,t2)
    elif scan_type=="angl":
       t1,t2,t3=map(int,tup.split())    
       tup=(t1,t2,t3)
    b,e,s=map(float,bes.split())
    bes=(b,e,s)
    return prog_type,scan_type,gopt_scan_fnameb,tup,bes

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
             sgopt_type,sgopt_scan_fnameb=sline
             gopt_type.append(sgopt_type.strip())
             gopt_scan_fnameb.append(sgopt_scan_fnameb.strip())
             tup.append(())
             bes.append(())
          else:
             sgopt_type,sgopt_scan_fnameb,stup,sbes=sline
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
             b,e,s=map(float,sbes.split())    
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
          scan_type=scan_type.strip()
          tup=[(),]
          bes=[(),]
       else:
          if len(line1)==4:
             gopt_type,gopt_scan_fnameb,tup,bes=line1
             scan_type="diha"
          else:
             gopt_type,gopt_scan_fnameb,tup,bes,scan_type=line1
          gopt_type=[gopt_type.strip(),]
          gopt_scan_fnameb=[gopt_scan_fnameb.strip(),]
          scan_type=scan_type.strip()
          if scan_type=="diha":
             t1,t2,t3,t4=map(int,tup.split())
             tup=[(t1,t2,t3,t4),]
          elif scan_type=="bond":
             t1,t2=map(int,tup.split())    
             tup=[(t1,t2),]
          elif scan_type=="angl":
             t1,t2,t3=map(int,tup.split())    
             tup=[(t1,t2,t3),]
          b,e,s=map(float,bes.split())
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
    #alg=lines[n+3].strip()
    ls=lines[n+3].split()
    lls=len(ls)
    if lls==1:
       alg=ls[0]
       ref_p=-1
    elif lls==2:
       alg=ls[0]
       ref_p=ls[1]
       if ref_p=="min": ref_p=-1
    else:
       print("Warning: more variables on the 'alg' line than allowed!")
    for line in lines[n+4:-1]:
        t=line.split()
        n=len(t)-1
        for i in range(n):
            sflag=False
            t[i+1]=t[i+1].lower()
            if t[i+1][0]=='-':
                sflag=True
            elif t[i+1][0]=='+':
                sflag=True
            if sflag:
               if t[i+1][1]=='p': 
                  np+=1
            else:
               if t[i+1][0]=='p': 
            #if t[i+1][0]=='p': 
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
           #opt_lin.update({int(t[0])-1:(t[1],t[2],t[3])})
           opt_lin.update({int(t[0])-1:[t[1],t[2],t[3]]})
        elif n==2:
           #opt_lin.update({int(t[0])-1:(t[1],t[2])})
           opt_lin.update({int(t[0])-1:[t[1],t[2]]})
        elif n==1:
           #opt_lin.update({int(t[0])-1:(t[1],)})
           opt_lin.update({int(t[0])-1:[t[1],]})
    last_line=lines[-1].split()
    if len(last_line)==1:
       step_int=10
       csv=lines[-1].strip()
    elif len(last_line)==2:
       csv,step_int=lines[-1].split()
    else:
       print("Wrong csv line in the input file!")

    return scan_type,gopt_type,gopt_scan_fnameb,tup,bes,engine_path,mm,mode,alg,ref_p,opt_lin,np,nc,step_int,csv

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
        print("read_add: Wrong MM-type!")

    f=open(add_name,'r')
    lines=f.readlines()
    f.close()
        
    i=j=0
    ol_keys=opt_lin.keys()
    ol_keys = sorted(ol_keys)
    pdict={}
    ol_templ={}
    scan_type=scan_type.strip()
    if scan_type=="diha":
       n=3
    elif scan_type=="bond":
       n=2
    elif scan_type=="angl":
       n=2
    constr_count=0
    for ol_sk in ol_keys:
        s=""
        s2=""
        t=opt_lin[ol_sk]
        minuses=[]
        for k in range(n):
            mflag=False
            if t[k][0]=='-':
                t[k]=t[k][1:]
                minuses.append(k)
                mflag=True
            elif t[k][0]=='+':
                t[k]=t[k][1:]
            if t[k][0]=='p':
                plab=t[k][1:]
                if len(plab)>0:
                   if plab in pdict:
                      ii=pdict[plab]
                      constr_count+=1
                   else:
                      pdict.update({plab:i})
                      ii=i
                      i+=1
                else:
                   ii=i
                   i+=1
                s+="p["+str(ii)+"],"
                if mflag:
                   s2+="-p["+str(ii)+"],"
                else:
                   s2+="p["+str(ii)+"],"
                #s+="p["+str(i)+"],"
                #i+=1
            elif t[k]=='c':
                s+="c["+str(j)+"],"
                s2+="c["+str(j)+"],"
                j+=1
            else:
                print("Woops!") 
        v=lines[ol_sk].split()
        ol_vars=s[:-1]
        ol_vars2=s2[:-1]
        if scan_type=="diha":
           v[5:10:2]=map(float,v[5:10:2])
           #v[4:9:2]=map(float,v[4:9:2])
           for k in minuses:
           #    print v[5+2*k]
              v[5+2*k]=-v[5+2*k]  
           #   v[5+2*k]=str(v[5+2*k])  
           v[5:10:2]=list(map(str,v[5:10:2]))
           exec(ol_vars+"="+v[5]+","+v[7]+","+v[9])
           #print "minuses",ol_vars
           v[1:5]=list(map(int,v[1:5]))
           s1='"%s        %3i  %3i  %3i  %3i'%(v[0],v[1],v[2],v[3],v[4])
           s2='      %6.3f +1   %6.3f -2 %6.3f +3      "%('+ol_vars2+')'
        elif scan_type=="bond":
           #v[3:5]=map(float,v[3:5])
           #for k in minuses:
           #   v[3+k]=-float(v[3+k])  
           #   v[3+k]=str(v[3+k])  
           #v[3:5]=map(str,v[3:5])
           exec(ol_vars+"="+v[3]+","+v[4])
           v[1:3]=list(map(int,v[1:3]))
           s1='"%s        %3i  %3i'%(v[0],v[1],v[2])
           s2='      %6.3f    %6.3f        "%('+ol_vars2+')'
        elif scan_type=="angl":
           #v[4:6]=map(float,v[4:6])
           #for k in minuses:
           #   v[4+k]=-float(v[4+k])  
           #   v[4+k]=str(v[4+k])  
           #v[4:6]=map(str,v[4:6])
           exec(ol_vars+"="+v[4]+","+v[5])
           v[1:4]=list(map(int,v[1:4]))
           s1='"%s        %3i  %3i  %3i'%(v[0],v[1],v[2],v[3])
           s2='      %6.3f    %6.3f        "%('+ol_vars2+')'
        ol_templ.update({ol_sk:s1+s2})
    
    if constr_count>0:
       p=p[:-constr_count] 
   
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
       print("write_add: Wrong MM-type!")

    f=open(add_name,'w')

    ol_keys=ol_templ.keys()
    ol_keys = sorted(ol_keys)
    for ol_sk in ol_keys:
        exec('lines['+str(ol_sk)+']='+ol_templ[ol_sk])

    for line in lines:
       print(line[:-1], file=f)

    f.close()
   
    if step is None: return
    
    if step=="ga" or (step-1)%step_int==0:
       if mm=="mm3":
          add_name_arch=sdir[0]+"/add_MM3_"+str(step)+".prm"
       elif mm=="mmff94":
          add_name_arch=sdir[0]+"/add_MMFF94_"+str(step)+".prm"

       f=open(add_name_arch,'w')
       for line in lines:
          print(line[:-1], file=f)
       f.close()

    return

if __name__=="__main__":
    gmol=DihGOpt_Molecule(mm="mmff94",name="opmmm")    
    print(gmol.name)
    #gmol.dt=(13,4,0,1)
    gmol.dt=(1,0,4,13)
    print(gmol.dt)
    fn_b="opmmm-mp2-opt"
    gmol.read_gopt_log(fn_b)
    print("tl:", gmol.tl)
    print(gmol.e)
    print(gmol.sl)
    print(gmol.rl)
    print(gmol.cl)
    print(gmol.da)
    gmol.dih_rot(pi/2.0)
    print(gmol.rl)
    print(gmol.da)
    print(gmol.write_inp_pcm("test"))
    gmol=DihGOpt_Molecule(mm="mmff94")
    gmol.read_ginp("tempC2")    
    gmol.dt=(13,4,0,1)
    print(gmol._ginp_templ)
    print(gmol._rl)
    gmol.dih_rot(pi/2.0)
    print(gmol._rl)
