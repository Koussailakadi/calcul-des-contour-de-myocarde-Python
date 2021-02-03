import csv
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from HQ import HQ

class FitData(HQ):
    def __init__(self):
        super().__init__()
        self.x=[]
        self.y=[] 
        self.a=None
        self.b=None
    
    def app(self):
        root= tk.Tk()
        canvas1 = tk.Canvas(root, width = 300, height = 300, bg = 'lightsteelblue2', relief = 'raised')
        canvas1.pack()
        browseButton_CSV = tk.Button(text=" Import CSV File ", command=self.getCSV, bg='green',fg='white', font=('helvetica', 12, 'bold'))
        canvas1.create_window(150, 150, window=browseButton_CSV)
        root.mainloop()
        
    def getCSV(self):
        x,y=[],[]
        import_file_path = filedialog.askopenfilename()
        file=open(import_file_path)
        data=csv.reader(file)
        for i,row in enumerate(data):
            if i==0: x=row
            else:    y=row
        self.x=np.asarray([float(i) for i in x])
        self.y=np.asarray([float(i) for i in y])
        
        
    def plotData(self):
        plt.figure(figsize=(5,5))
        plt.scatter(self.x,self.y,marker='+',c='b')
        plt.title('data scatter')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()  
        
    #---------------------------------------------
    def w(self,x=None,y=None):
        s=(self.a*x+self.b*y)**4+(x+y)**4-1
        return s
    
    def J(self,a,b):
        s=0
        for i in range(len(self.x)):
               s=s+((a*self.x[i]+b*self.y[i])**4+(self.x[i]+self.y[i])**4-1)**2
        return s
    
    
    # methode de gradient: 
    def Gradient_descent(self,a=0.5,b=-0.5,itermax=1000):
        a_ini=a
        b_ini=b
        alpha=0.001
        iter=0
        eps=10e-6
        sa=sb=1
        while  (abs(sa)>eps or abs(sb)>eps) and iter<itermax :#self.J(a,b)>0.3:
            sa,sb=0,0
            for i in range(len(self.x)):
                sa+=8*self.x[i]*((a*self.x[i]+b*self.y[i])**3)*((a*self.x[i]+b*self.y[i])**4+(self.x[i]+self.y[i])**4-1)
                sb+=8*self.y[i]*((a*self.x[i]+b*self.y[i])**3)*((a*self.x[i]+b*self.y[i])**4+(self.x[i]+self.y[i])**4-1)
            a=a-alpha*sa
            b=b-alpha*sb
            iter+=1
        self.a=a
        self.b=b
        self.A=[a,1]
        self.B=[b,1]
        self.C=[0,0]
        self.k=2
        mes='\n---------------\nGradient descent:\n---------------\n'
        mes+='a_ini={} , b_ini= {}\n'.format(a_ini,b_ini)
        mes+='a_final={} , b_final= {}\n'.format(a,b)
        mes+='Nombre iter={}\n'.format(iter)
        mes+='valeur final J ={}'.format(self.J(a,b))
        print(mes)
        return mes
    
    def Newton(self,a=0.5,b=-0.5,itermax=1000):
        a_ini=a
        b_ini=b
        dX=1 #longueur du déplacement
        eps = 10e-6
        iter=0
        d1a,d1b,d2a,d2b,d2ab=0,0,0,0,0
        X = np.array([[a],[b]])
        while dX>eps and  iter<itermax:
            for i in range(len(self.x)):
                a=X[0][0]
                b=X[1][0]
                d1a+=8*self.x[i]*((a*self.x[i]+b*self.y[i])**3)*((a*self.x[i]+b*self.y[i])**4+(self.x[i]+self.y[i])**4-1)
                d1b+=8*self.y[i]*((a*self.x[i]+b*self.y[i])**3)*((a*self.x[i]+b*self.y[i])**4+(self.x[i]+self.y[i])**4-1)
                d2a+=32*self.x[i]**2*(a*self.x[i]+b*self.y[i])**6 + 24*self.x[i]**2*(a*self.x[i]+b*self.y[i])**2*((a*self.x[i]+b*self.y[i])**4 + (self.x[i]+self.y[i]**4)-1)
                d2b+=32*self.y[i]**2*(a*self.x[i]+b*self.y[i])**6 + 24*self.x[i]**2*(a*self.x[i]+b*self.y[i])**2*((a*self.x[i]+b*self.y[i])**4 + (self.x[i]+self.y[i]**4)-1)
                d2ab+= 32*self.x[i]*self.y[i]**2*(a*self.x[i]+b*self.y[i])**6 + 24*self.x[i]**2*(a*self.x[i]+b*self.y[i])**2*((a*self.x[i]+b*self.y[i])**4 + (self.x[i]+self.y[i]**4)-1)
            Grad = np.array([[d1a],[d1b]])
            Hess = np.array([[d2a,d2ab],[d2ab,d2b]])
            deltaX =np.dot(np.linalg.inv(Hess),-Grad)
            X += deltaX
            dX = np.linalg.norm(deltaX)
    
            iter+=1
        self.a=X[0][0]
        self.b=X[1][0]
        self.A=[self.a,1]
        self.B=[self.b,1]
        self.C=[0,0]
        self.k=2

        mes='\n---------------\nNewton:\n---------------\n'
        mes+='a_ini={} , b_ini= {}\n'.format(a_ini,b_ini)
        mes+='a_final={} , b_final= {}\n'.format(self.a,self.b)
        mes+='Nombre iter={}\n'.format(iter)
        mes+='valeur final J ={}'.format(self.J(self.a,self.b))
        print(mes)
        return mes

    def graph(self):
        self.Droite_englobantes(f=self.w,xmin=-1.5,xmax=1.5,ymin=-1.5,ymax=1.5,nx=50,ny=50,figsize=(5,5))
        plt.scatter(self.x,self.y,marker='+',c='b')
        plt.title('graph')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()

