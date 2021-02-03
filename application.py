from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from FitData import *

class application(FitData):
    def __init__(self):
        super().__init__()
        self.root=Tk()

    def quit_app(self):
        self.root.destroy()


    def app(self):
        self.root.title('Hyperquadrique et mesure de contours')
        self.root.geometry("1080x800")
        self.root.minsize(1080,800)
        self.root.maxsize(1080,800)
        self.root.iconbitmap("logo.ico")
        self.root.config(background='#34495E')

        frame1 = Frame(self.root, relief=GROOVE, width=100, height=50, bd=1)
        frame1.place(x=2, y=2)

        img = PhotoImage(file="img.png")
        img_zone=Canvas(frame1,bg='#B3B6B7',height=50, width=100)
        CAN_Zone_Image = img_zone.create_image(0,0, image=img,anchor="nw")
        img_zone.grid()

        frame2= Frame(self.root,bg='#34495E',width=100, height=50)#'#5D6D7E'
        frame2.place(x=440,y=2)
        label_title = Label(frame2,text='Projet d\'Optimisation', font=('courrier', 20), bg='#34495E',fg='#F8F9F9').pack()

        frame3= Frame(self.root,bg='#B3B6B7')
        frame3.place(x=500,y=70)
        team= Button(frame3, text='Réalisé par', font=('courrier', 15),bg='#B3B6B7', fg='black', command=self.team)
        team.pack(fill=X)

        frame5 = Frame(self.root, bg='#B3B6B7')
        frame5.place(x=505, y=130)
        Instruction = Button(frame5, text='Instruction', font=('courrier', 15),bg='#B3B6B7', fg='black',
                      command=self.instraction)
        Instruction.pack(fill=X)

        frame7 = Frame(self.root, bg='#B3B6B7')
        frame7.place(x=440, y=200)
        Methode_label=Label(frame7,text='Choisir une methode',font=('courrier', 20), bg='#34495E',fg='#F8F9F9').pack()

        frame8 = Frame(self.root, bg='#B3B6B7')
        frame8.place(x=300, y=250)
        Newton = Button(frame8, text='Methode Newton', font=('courrier', 15), bg='#B3B6B7', fg='black',
                      command=self.PlotNewton)
        Newton.pack(fill=X)

        frame9 = Frame(self.root, bg='#B3B6B7')
        frame9.place(x=700, y=250)
        Newton = Button(frame9, text='Gradient descent', font=('courrier', 15), bg='#B3B6B7', fg='black',
                        command=self.PlotGradient)
        Newton.pack(fill=X)

        frame14=Frame(self.root,bg='#B3B6B7')

        frame14.place(x=900,y=50)
        bouton_quit = Button(frame14, text='Quit Application', font=('courrier', 15), bg='#EC7063' , fg='black',
                        command=self.quit_app)
        bouton_quit.pack(fill=X)



        self.root.mainloop()


    def team(self):
        res='--------------------------------------------------------------------------------\n'
        res+='But de projet: calcul des contours du myocarde en utilisant les hyperquadriques\n'
        res+='--------------------------------------------------------------------------------\n'
        res+='réalisé par:\n'
        res+='-Koussaila KADI\n'
        res+='-Nathaniel DAHON\n'
        res+='-mehdi ZEMOUM\n'
        res+='-reponsable TP: Mr. BETOUX\n'
        res+='Date:2021\n'
        res+='Sorbonne université\nMerci!'
        frame4 = Frame(self.root, bg='#F8F9F9')
        frame4.place(x=650, y=70)
        team_Label = Label(frame4, text=res).pack()
        bouton=Button(frame4,text='QUIT', font=('courrier', 10),bg='#B3B6B7', fg='black', command=frame4.destroy).pack()
        del(res)

    def instraction(self):
        res = '------------------------------------------------------------------\n'
        res += ' calcul des contours du myocarde en utilisant les hyperquadriques\n'
        res += '-----------------------------------------------------------------\n'
        res += 'comment bien utiliser l\'application?\n'
        res += 'choisissez une methode, une fenêtr va apparaitre pour selectionner le fichier\n'
        res += 'CSV qui est dans le même dossier que l\'application,ce fichier contient les \n'
        res += 'données réel de notre contour de myocarde, une fois c\'est fait, il y\'aura les\n'
        res += 'de chaque méthode. à la fin vous cliquer sur quitter application pour sortir\n'

        frame6 = Frame(self.root, bg='#F8F9F9')
        frame6.place(x=650, y=100)
        team_Label = Label(frame6, text=res).pack()
        bouton = Button(frame6, text='QUIT', font=('courrier', 10), bg='#B3B6B7', fg='black',
                        command=frame6.destroy).pack()
        del(res)

    def PlotNewton(self,xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, nx = 50, ny = 50):
        res='Methode de Newton:'
        frame6 = Frame(self.root, bg='#F8F9F9')
        frame6.place(x=150, y=300)
        team_Label = Label(frame6, text=res).pack()

        data = FitData()
        data.getCSV()
        mes=data.Newton()

        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        #------------------------------------------
        tab_x = np.linspace(xmin, xmax, nx)
        tab_y = np.linspace(ymin, ymax, nx)
        # --------------------------------
        x1d = np.linspace(xmin, xmax, nx)
        y1d = np.linspace(ymin, ymax, ny)
        x2d, y2d = np.meshgrid(x1d, y1d)

        for i in range(data.k):
            if data.B[i] == 0 and data.A[i] != 0:
                try:
                    ax.plot([data.f1(i), data.f1(i)], [xmin, xmax], 'r-', lw=2)  # Red straight line
                    ax.plot([data.f2(i), data.f2(i)], [ymin, ymax], 'k-', lw=2)  # Black straight line
                except:
                    pass

        for i in range(data.k):
            if data.B[i] != 0:
                try:
                    plt.plot([xmin, xmax], [data.f3(xmin, i), data.f3(xmax, i)], 'b-', lw=2)  # Blue straight line
                    plt.plot([xmin, xmax], [data.f4(ymin, i), data.f4(ymax, i)], 'g-', lw=2)  # Green straight line
                except:
                    pass

        ax.contour(x2d, y2d, data.w(x2d, y2d), cmap='rainbow')
        #------------------------------------------
        ax.scatter(data.x,data.y)
        ax.set_xlim([-2,2])
        ax.set_ylim([-2, 2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid()
        canvas=FigureCanvasTkAgg(fig,frame6)
        canvas.get_tk_widget().pack()
        bouton = Button(frame6, text='QUIT', font=('courrier', 10), bg='#B3B6B7', fg='black',
                        command=frame6.destroy).pack()

        frame10 = Frame(self.root, bg='#F8F9F9')
        frame10.place(x=500, y=680)
        team_Label = Label(frame6, text=mes).pack()
        del(data)
        del(mes)

    def PlotGradient(self,xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5, nx = 50, ny = 50):
        res = 'Gradient Descent:'
        frame6 = Frame(self.root, bg='#F8F9F9')
        frame6.place(x=650, y=300)
        team_Label = Label(frame6, text=res).pack()

        data = FitData()
        data.getCSV()
        mes = data.Gradient_descent()

        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        # ------------------------------------------
        tab_x = np.linspace(xmin, xmax, nx)
        tab_y = np.linspace(ymin, ymax, nx)
        # --------------------------------
        x1d = np.linspace(xmin, xmax, nx)
        y1d = np.linspace(ymin, ymax, ny)
        x2d, y2d = np.meshgrid(x1d, y1d)

        for i in range(data.k):
            if data.B[i] == 0 and data.A[i] != 0:
                try:
                    ax.plot([data.f1(i), data.f1(i)], [xmin, xmax], 'r-', lw=2)  # Red straight line
                    ax.plot([data.f2(i), data.f2(i)], [ymin, ymax], 'k-', lw=2)  # Black straight line
                except:
                    pass

        for i in range(data.k):
            if data.B[i] != 0:
                try:
                    plt.plot([xmin, xmax], [data.f3(xmin, i), data.f3(xmax, i)], 'b-', lw=2)  # Blue straight line
                    plt.plot([xmin, xmax], [data.f4(ymin, i), data.f4(ymax, i)], 'g-', lw=2)  # Green straight line
                except:
                    pass

        ax.contour(x2d, y2d, data.w(x2d, y2d), cmap='rainbow')
        # ------------------------------------------
        ax.scatter(data.x, data.y)
        ax.set_xlim([-2, 2])
        ax.set_ylim([-2, 2])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.grid()
        canvas = FigureCanvasTkAgg(fig, frame6)
        canvas.get_tk_widget().pack()
        bouton = Button(frame6, text='QUIT', font=('courrier', 10), bg='#B3B6B7', fg='black',
                        command=frame6.destroy).pack()

        frame11 = Frame(self.root, bg='#F8F9F9')
        frame11.place(x=540, y=680)
        team_Label = Label(frame6, text=mes).pack()
        del(data)
        del(mes)





