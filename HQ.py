import numpy as np
import matplotlib.pyplot as plt
class HQ:
    def __init__(self):
        self.k=3
        self.A=[1,0,0.7]
        self.B=[0,1,-0.7]
        self.C=[0,0,0]
        self.G=[5,5,5]
       
    # fonction 
    def f(self,x=None,y=None):
        s=0
        for i in range (self.k):
            s=s+abs(self.A[i]*x+self.B[i]*y+self.C[i])**self.G[i]
        return s-1  #pour le contour recherché
    
    def f1(self,i):
        return -self.C[i]/self.A[i]-1/self.A[i]
        
    def f2(self,i):
        return -self.C[i]/self.A[i]+1/self.A[i]
        
    def f3(self,x,i):
        return (-self.A[i]/self.B[i])*x-self.C[i]/self.B[i]-1/self.B[i]
        
    def f4(self,x,i):
        return (-self.A[i]/self.B[i])*x-self.C[i]/self.B[i]+1/self.B[i]
    
        
    def Droite_englobantes(self,f,xmin=-1.5,xmax=1.5,ymin=-1.5,ymax=1.5,nx=50,ny=50,figsize=(5,6)):
        plt.figure(figsize=figsize)
        tab_x = np.linspace(xmin,xmax,nx)
        tab_y = np.linspace(ymin,ymax,nx)
        #--------------------------------
        x1d = np.linspace(xmin,xmax,nx)
        y1d = np.linspace(ymin,ymax,ny)
        x2d, y2d = np.meshgrid(x1d, y1d)
        
        for i in range(self.k):
            if self.B[i]==0 and self.A[i]!=0: 
                try:
                    plt.plot([self.f1(i) ,self.f1(i)],[xmin, xmax], 'r-', lw=2) # Red straight line
                    plt.plot([self.f2(i) ,self.f2(i)],[ymin, ymax], 'k-', lw=2) # Black straight line
                except: pass
                
        for i in range(self.k):
            if self.B[i]!=0:
                try:
                    plt.plot([xmin, xmax],[self.f3(xmin,i),self.f3(xmax,i)], 'b-', lw=2) # Blue straight line
                    plt.plot([xmin, xmax],[self.f4(ymin,i),self.f4(ymax,i)], 'g-', lw=2) # Green straight line
                except: pass
        
        plt.contour(x2d,y2d,f(x2d,y2d), cmap = 'rainbow')
        plt.title('les doites englobantes')
        plt.xlim([xmin, xmax])
        plt.ylim([ymin, ymax])


    
        

