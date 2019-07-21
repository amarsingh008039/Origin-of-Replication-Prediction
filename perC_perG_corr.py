import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt2

from numpy import *
char=[]
l=0;lob=0;inc=0;gc_mean=0;co_mean=0
def check():
    a=0;t=0;n=0;g=0;c=0
    global l
    f=open("D:/DNA/DNA_Seq/DNA_Seq/NC_000911.fna",'r')
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
        for k in range(20):
            ck=abs(func(k,r))
            c=c+ck
        temp=c/(20)
        cg.append(temp)
perC=[]
perG=[]
def count_perC_perG():
    i=0
    global perC;global perG
    for i in range(0,len(char)-lob,i+inc):
        g=0;c=0
        for j in range(i,i+lob):
            r=char[j]
            if r is 'G':
                g=g+1
            if r in 'C':
                c=c+1
        p=(g/lob)*100
        q=(c/lob)*100

        perG.append(p)
        perC.append(q) 
    
    
corrC=[]
corrG=[]
kc=0
perc_mean=0;co_mean=0;perg_mean=0
def correlation_per(l1):
    global perc_mean;global co_mean;global perg_mean
    global corr;global kc
    global perC;global perG;global corrC;global corrG
    l1=len(perC)
    perc_mean=sum(perC)/float(l1)
    co_mean=sum(cg)/float(l1)
    perg_mean=sum(perG)/float(l1)
    if(l1>200):
        kc=450
    else:
        kc=100
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((cg[j]-co_mean)*(perC[j+i]-perc_mean))
            t=t+temp
        c=t/(l1-i)
        corrC.append(c)
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((cg[j]-co_mean)*(perG[j+i]-perg_mean))
            t=t+temp
        c=t/(l1-i)
        corrG.append(c)

corrCG=[]
def corr_perC_perG(l1):
    global perc_mean;global co_mean;global perg_mean
    global kc
    global perC;global perG
    
    perc_mean=sum(perC)/float(l1)
    perg_mean=sum(perG)/float(l1)
    if(l1>200):
        kc=450
    else:
        kc=100    
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((perC[j]-co_mean)*(perG[j+i]-perc_mean))
            t=t+temp
        c=t/(l1-i)
        corrCG.append(c)

def _plot():
    check()
    global lob
    global inc
    lob=(l//100)
    inc=(lob//5)
    '''print("Genome Length: -  ",l)
    print("Window Size: - ",lob)
    print("Increment : - ",inc)
    print("-------GC SKEW calculation-------")
    gc_window()
    l1=len(GC)'''
    '''x=list(range(0,l1))
    plt.plot(x,gc_skew)
    #plt1.xticks(list(range(0,2000,100)))
    plt.xlabel('Window Number')
    plt.ylabel('G-C Skew')
    plt.show()'''
    print("-------Correlation Measure Calculation-------")
    corr_G()
    #l2=len(G_seq)
    '''x=list(range(0,l1))
    plt1.plot(x,cg)
    plt1.xlabel('Window Number')
    plt1.ylabel('Correlation')
  #  plt1.xticks(list(range(0,l,20)))
  #  plt1.yticks([0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09])
    plt1.show()
    '''
    print("-------Comparison Correlation Calculation-------")
    global corr;global kc
    global perC;global perG
    count_perC_perG()
    l1=len(perC)
    correlation_per(l1)
    x=range(kc)
    plt1.plot(x,corrC)
    plt1.xlabel('k value')
    plt1.ylabel(' Correlation %C')
    plt1.show()
    plt2.plot(x,corrG)
    plt2.xlabel('k value')
    plt2.ylabel('Correlation %G')
    plt2.show()
    corr_perC_perG(l1)
    plt.plot(x,corrCG)
    plt.xlabel('k value')
    plt.ylabel('Correlation %G%C')
    plt.show()
    
if __name__=='__main__':
    _plot()
    
