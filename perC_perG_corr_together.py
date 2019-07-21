import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt2
from numpy import *

char=[]
l=0;lob=0;inc=0;gc_mean=0;co_mean=0;n=[]
files=[]
    
def check(f2):
    global char
    char=[]
    global l;
    n1=0
    f=open(f2,'r')
    l=len(f.read())
    f.seek(0)
    for i in range(l):  
                r=f.read(1)
                if r is '\n':
                    n1=n1+1
                    continue
                char.append(r)
    l=l-n1      
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

def corr_G1():
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
all_corrC=[]
corrG=[]
all_corrG=[]
kc=0
perc_mean=0;co_mean=0;perg_mean=0
def correlation_per(l1):
    global perc_mean;global co_mean;global perg_mean;global kc
    global perC;global perG;global corrC;global corrG
    corrG=[];corrC=[]
    l1=len(perC)
    perc_mean=sum(perC)/float(l1)
    co_mean=sum(cg)/float(l1)
    perg_mean=sum(perG)/float(l1)
    
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((cg[j]-co_mean)*(perC[j+i]-perc_mean))
            t=t+temp
        c=t/(l1-i)
        corrC.append(c)
    all_corrC.append(corrC)
    
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((cg[j]-co_mean)*(perG[j+i]-perg_mean))
            t=t+temp
        c=t/(l1-i)
        corrG.append(c)
    all_corrG.append(corrG)
    

corrCG=[]
all_corrCG=[]
def corr_perC_perG(l1):
    global perc_mean;global co_mean;global perg_mean
    global kc
    global perC;global perG
    corrCG=[]
    perc_mean=sum(perC)/float(l1)
    perg_mean=sum(perG)/float(l1)
      
    for i in range(kc):
        t=0
        for j in range(l1-i):
            temp= ((perC[j]-co_mean)*(perG[j+i]-perc_mean))
            t=t+temp
        c=t/(l1-i)
        corrCG.append(c)
    all_corrCG.append(corrCG)
    

def _plot(a,c):
    global kc
    x=range(kc)
    y=range(100)
    plt1.plot(y,a[0],'#ff0000')
    plt1.plot(x,a[1],'#000000')
    plt1.plot(x,a[2],'#00FF00')
    plt1.plot(x,a[3],'#00008B')
    plt1.plot(x,a[4],'#FFFF00')
    plt1.plot(x,a[5],'#FF00FF')
    plt1.plot(x,a[6],'#FF8C00')
    plt1.plot(x,a[7],'#808000')
    plt1.plot(x,a[8],'#00CED1')
    plt1.plot(x,a[9],'#B8860B')
    plt1.plot(x,a[10],'#cd5c5c')
    plt1.plot(x,a[11],'#8b4513')
    plt1.plot(x,a[12],'#708090')
    plt1.plot(x,a[13],'#ee82ee')
    plt1.plot(x,a[14],'#fa8072')
    plt1.plot(x,a[15],'#bdb76b')
    plt1.plot(x,a[16],'#a020f0')
    plt1.plot(x,a[17],'#bc8f8f')
    
    plt1.xlabel('k value')
    plt1.ylabel(c)
    plt1.show()
    
    
def calc_corrs():
    global lob
    global inc
    global n;global kc
    global perC;global perG
    global G_seq;global cg
    for i in range(18):
        f=input("Enter file->")
        files.append(f)
        n1=int(input("Enter the number to divide N:- "))
        n.append(n1)
    for j in range(18):
        f=files[j]
        check(f)
        n1=n[j]
        if (n1==40):
            kc=100
        else:
            kc=200
        
        lob=(l//n1)
        inc=(lob//5)
        perG=[]
        perC=[]
        G_seq=[]
        cg=[]
        print("File ",j+1)
        print("Genome Length: -  ",l)
        print("Window Size: - ",lob)
        print("Increment : - ",inc)
        count_perC_perG()
        if(j==1):
            print("-------Correlation Measure Calculation-------")
            corr_G1()
            l1=len(perC)
            print("------Correlation with %C and with %G------")
            correlation_per(l1)
            print("------Correlation between %C and  %G------")
            corr_perC_perG(l1)
            continue
        print("-------Correlation Measure Calculation-------")
        corr_G()
        l1=len(perC)
        print("------Correlation with %C and with %G------")
        correlation_per(l1)
        print("------Correlation between %C and  %G------")
        corr_perC_perG(l1)
        
        
    
if __name__=='__main__':
    calc_corrs()
    _plot(all_corrC,'Correlation %C')
    _plot(all_corrG,'Correlation %G')
    _plot(all_corrCG,'Correlation %C_%G')
    
