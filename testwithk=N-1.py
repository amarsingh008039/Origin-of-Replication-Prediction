import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
from numpy import *
char=[]
l=0;lob=0;inc=0
def check():
    a=0;t=0;n=0;g=0;c=0
    global l
    f=open("D:/DNA/DNA_Seq/DNA_Seq/NC_008207.fna",'r')
    l=len(f.read())
    f.seek(0)
    for i in range(l):
                #print(f.read(1),end='')    
                r=f.read(1)
                if r is 'A':
                    a=a+1
                if r is 'G':
                    g=g+1
                if r is 'T':
                    t=t+1            
                if r is 'C':
                    c=c+1
                if r is '\n':
                    n=n+1
                    continue
                char.append(r)
    l=l-n
    print("Following file consists :-")
    print("% of A ",str(round(((a/l)*100),2)))
    print("% of T ",str(round(((t/l)*100),2)))                
    print("% of G ",str(round(((g/l)*100),2)))                
    print("% of C ",str(round(((c/l)*100),2)))                
    f.close()
    l=len(char)
    
GC=[]
gc_skew=[]
def gc_window():
    i=0
    for i in range (0,len(char)-lob,i+inc):
        nC=0;
        nG=0;
        for j in range(i,i+lob):
            if char[j] is 'C':
                    nC=nC+1
            if char[j] is 'G':
                    nG=nG+1
        GC.append((nG,nC))
    i=len(GC)
    for j in range(i):
        g=GC[j][0]
        c=GC[j][1]
        skew=(c-g)/(c+g)
        gc_skew.append(skew)


def func(k,r):
    t=0
    for j in range(len(r)-k):
        temp=r[j]*r[j+k]
        t=t+temp
    ck=t/(len(r)-k)
    return ck

G_seq=[]
cg=[]
def corr_G():
    i=0
    for i in range(0,len(char)-lob,i+inc):
        g=[]
        for j in range(i,i+lob):
            r=char[j]
            if r is 'G':
                g.append(1)
            elif r in ['A','T','C']:
                g.append(-1)
            else:
                g.append(0)
        G_seq.append(g)
    
    for i in range(len(G_seq)):
        r=G_seq[i]
        c=0
        for k in range(len(r)-1):
            ck=abs(func(k,r))
            c=c+ck
        temp=c/(len(r)-1)
        #print (temp)
        cg.append(temp)
    

def _plot():
    check()
    global lob
    global inc
    lob=l//100
    inc=lob//5
    
    print(lob)
    print(inc)
    print("-------GC SKEW-------")
    gc_window()
    l1=len(GC)
    x=list(range(0,l1))
    plt.plot(x,gc_skew)
    #plt1.xticks(list(range(0,2000,100)))
    plt.xlabel('Window Number')
    plt.ylabel('G-C Skew')
    plt.show()
    print("-------Correlation Method-------")
    corr_G()
    l1=len(G_seq)
    x=list(range(0,l1))
    plt1.plot(x,cg)
    plt1.xlabel('Window Number')
    plt1.ylabel('Correlation')
  #  plt1.xticks(list(range(0,l,20)))
  #  plt1.yticks([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09])
    plt1.show()
    

if __name__=='__main__':
    _plot()
    
