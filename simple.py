import queue
import random
import sys
from turtle import color
from numpy import arange, sin, pi, array
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.uic import loadUi
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_template import FigureCanvas
from matplotlib.figure import Figure
import ressources_rc
import sys
from PyQt5 import QtWidgets, uic, QtGui, QtCore
from PyQt5.QtWidgets import QDialog, QApplication, QMainWindow, QMessageBox, QLabel, QGraphicsItem, QGraphicsView, \
    QGraphicsScene, QSizePolicy, QVBoxLayout, QGraphicsDropShadowEffect
import  matplotlib
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import copy
import sys
import random
import math
import copy
import time
import numpy as np
from decimal import *

#-----------------------------------------------------------------------
def random_energie( tab ):
     for i in range(len(tab)):
         tab[i] = round(random.uniform(0.1 , 10.0) , 2)
     return tab
def randomnotsame(tab , nbr_rand):
	j=0
	lst = []
	while j < nbr_rand :
		rnd =random.randrange(0,len(tab))
		if rnd not in lst:
			lst.append(rnd)
			tab[rnd] = 1
			j= j + 1
	return tab
#------------------------------ initialiser individus --------------
def initialiser_individue(nbr_BSs, width, height , nbsmin, rayon):
    actarray = nbr_BSs * [0]
    energy = nbr_BSs * [0.0]
    xarray = no_overlapping(nbr_BSs , width , height ,  rayon)
    actarray = randomnotsame(actarray ,nbsmin)
    energy = random_energie(energy)
    data = []
    for i in range(nbr_BSs):
        data.append((xarray[i][0], xarray[i][1] , actarray[i] , energy[i] ))
    return data

#-----------------------------------------------
#---------------------------------------------------------------------------------
def calculer_Scov(rayon):
	#area = math.pow(3 , 1/3)  * ((rayon ** 2)  / 2)
    area = rayon*rayon * math.pi
    return area
#--------------------------------------------------------------------------
def nbr_site(Scov , Stot):
	nbr = int ( Stot / Scov) +1
	return int(nbr * 1.8)+1
#-------------------------------------------------------------------------
def surface(width, height):
	return width*height
#--------------------------------------------------------------------------
def fitness_function(individu ,width, height,rayon,  nbr_BSs , nbsmin , T_p ,alpha, beta, lamda):
    Acov = verif_overlapp(individu ,int(width), int(height), rayon)
    alpha = 0.7
    beta = 0.2
    lamda = 0.1
    UE = 1000
    fitness_value = alpha * (Acov) + beta * (nbsmin/(nbr_BSs)) + lamda * T_p \
                   # * (UE**2 / (T_p * nbsmin ** 2))
    return fitness_value


def onlyactif(individu):
    result=[]
    for i in range(len(individu)):
        if individu[i][2] == 1:
            result.append(individu[i])

    return  result

class wolf:
    def __init__(self,nbr_BSs, width, height, rayon,alpha, beta, lamda ):
        self.nbsmin =random.randint(int(nbr_BSs*0.2), nbr_BSs)
        # self.nbsmin =nbr_BSs
        self.individu =initialiser_individue(nbr_BSs, width, height, self.nbsmin , rayon)
        self.puisst = 0.0
        for i in range(len(self.individu)):
            if(self.individu[i][2]==1):
                self.puisst += self.individu[i][3]
        self.fitness = fitness_function(self.individu , width , height ,rayon, nbr_BSs, self.nbsmin , self.puisst,alpha, beta, lamda) # curr fitness
    def __repr__(self):
        return '{ ' + str(self.individu) + ',' + 'fitness = ' + str(self.fitness) +  'puissance totale ='\
               +  str(self.puisst)  +'} \n'

# grey wolf optimization (GWO)
def gwo(max_iter, nbr_individu, nbr_BSs, width, height, rayon, alpha, beta, lamda):

    population = [wolf(nbr_BSs, width, height, rayon,alpha, beta, lamda) for i in range(nbr_individu)]
    population = sorted(population, key=lambda x: x.fitness, reverse=True)
    alpha_wolf, beta_wolf, delta_wolf = copy.copy(population[: 3])
    #**************************************************************
    covrage = [0] * max_iter
    alphalist = [0] * max_iter
    Iter = 0
    temps = 0.0
    while Iter < max_iter:
        start= time.time()
        a = 2 - Iter * (2 / max_iter)
        for i in range(nbr_individu):
            Xnew = copy.copy(population[i])
            for j in range(nbr_BSs):
                lst = list(Xnew.individu[j])
                if lst[2] == 1:
                    r1x = random.random()
                    r2x = random.random()
                    r1y = random.random()
                    r2y = random.random()
                    A1, A2, A3 = round(a * (2 * r1x - 1), 4), round(a * (2 * r1x - 1), 4), round(a * (2 * r1x - 1), 4)
                    A1y, A2y, A3y = round(a * (2 * r1y - 1), 4), round(a * (2 * r1y - 1), 4), round(a * (2 * r1y - 1),
                                                                                                    4)
                    C1, C2, C3 = round(2 * r2x, 4), round(2 * r2x, 4), round(2 * r2x, 4)
                    C1y, C2y, C3y = round(2 * r2y, 4), round(2 * r2y, 4), round(2 * r2y, 4)

                    X1 = round(alpha_wolf.individu[j][0] - A1 * abs(
                        C1 * alpha_wolf.individu[j][0] - population[i].individu[j][0]), 2)
                    X2 = round(beta_wolf.individu[j][0] - A2 * abs(
                        C2 * beta_wolf.individu[j][0] - population[i].individu[j][0]), 2)
                    X3 = round(delta_wolf.individu[j][0] - A3 * abs(
                        C3 * delta_wolf.individu[j][0] - population[i].individu[j][0]), 2)
                    lst[0] = round((X1 + X2 + X3) / 3, 2)

                    Y1 = round(alpha_wolf.individu[j][1] - A1y * abs(
                        C1y * alpha_wolf.individu[j][1] - population[i].individu[j][1]), 2)
                    Y2 = round(beta_wolf.individu[j][1] - A2y * abs(
                        C2y * beta_wolf.individu[j][1] - population[i].individu[j][1]), 2)
                    Y3 = round(delta_wolf.individu[j][1] - A3y * abs(
                        C3y * delta_wolf.individu[j][1] - population[i].individu[j][1]), 2)
                    lst[1] = round((Y1 + Y2 + Y3) / 3, 2)
                    lst[0] = boundry(lst[0], width, rayon)
                    lst[1] = boundry(lst[1], height, rayon)
                    Xnew.individu[j] = tuple(lst)
                    population[i] = copy.copy(Xnew)
#********************************************************
            # Xnew.fitness = fitness_function(Xnew.individu, width, height, rayon, nbr_BSs, Xnew.nbsmin,Xnew.puisst)
            # if Xnew.fitness > population[i].fitness:
            #     population[i] = copy.copy(Xnew)
        for j in range(nbr_individu):
            population[j].fitness = fitness_function(population[j].individu, width, height, rayon, nbr_BSs,population[j].nbsmin, population[j].puisst,alpha, beta, lamda)
        for k in range(nbr_individu):
            fit = population[k].fitness
            if fit > alpha_wolf.fitness:
                delta_wolf = copy.copy(beta_wolf)
                delta_wolf.individu = copy.copy(beta_wolf.individu)
                beta_wolf = copy.copy(alpha_wolf)
                beta_wolf.individu = copy.copy(alpha_wolf.individu)
                alpha_wolf = copy.copy(population[k])
                alpha_wolf.individu = copy.copy(population[k].individu)
            if fit < alpha_wolf.fitness and fit > beta_wolf.fitness:
                delta_wolf = copy.copy(beta_wolf)
                delta_wolf.individu = copy.copy(beta_wolf.individu)
                beta_wolf = copy.copy(population[k])
                beta_wolf.individu = copy.copy(population[k].individu)
            if fit < alpha_wolf.fitness and fit < beta_wolf.fitness and fit > delta_wolf.fitness:
                delta_wolf = copy.copy(population[k])
                delta_wolf.individu = copy.copy(population[k].individu)


        #-*****************************************

        print( ["At iteration " + str(Iter) + " the best fitness is " ,(alpha_wolf.fitness)] )
#//-------------------------------------
        end = time.time()
        duree = end-start
        temps += duree
        alphalist  = copy.copy(alpha_wolf)
        covrage[Iter] = (temps, alpha_wolf.fitness)
        Iter = Iter + 1
    #------------------************************
    print("end while")
    return alphalist,covrage
#----------------------------------------------------------------------------------------------------------------------
def verif_overlapp( indv , w, h, rayon):

    ctrl = 10
    mat = np.zeros( (w*ctrl, h*ctrl) )
    #mat = np.array([[0 for y in range(w*ctrl)] for x in range(h*ctrl)] , dtype=np.byte)

    for k in range(len(indv)) :
        xi = (int)(indv[k][0] - rayon) * ctrl
        yi = (int)(indv[k][1] - rayon) * ctrl
        xm = (int)(indv[k][0]  + rayon) * ctrl
        ym = (int)(indv[k][1]  + rayon) * ctrl

        if (indv[k][2] == 1):
            for i in range(xi, xm):
                for j in range(yi, ym):
                    if (distance(indv[k][0] * ctrl, indv[k][1] * ctrl, i, j) <= rayon * ctrl):
                        #if (i >= 0 and i < w * ctrl and j >= 0 and j < h * ctrl):
                            mat[i][j] = 1

    cpt =0
    for i in range(w * ctrl):
        for j in range(h * ctrl):
            if(mat[i][j] ==1):
                cpt +=1
    #print("------------------- end while verif-------------" , cpt , w * ctrl* h * ctrl )
    return round((cpt) / (w *ctrl* h*ctrl), 3)
def nombre_site_actif(x):
    cpt =0
    for i in range(len(x.individu)):
        if x.individu[i][2]==1:
            cpt += 1
    return cpt
def distance(x1,y1,x2,y2):
    return  math.sqrt((x1-x2)**2 + (y1-y2)**2)
def no_overlapping(taille , w , h , rayon):
    circles= []
    overlap = False
    while len(circles) < taille:
        x = round(random.uniform(rayon, w - (rayon)), 2)
        y = round(random.uniform(rayon, h - (rayon)), 2)
        circles.append((x, y))
        overlap= False
        #for i in range(len(circles)):
        #     other = circles[i]
        #     if(distance(other[0], other[1], x, y) < rayon):
        #         overlap = True
        #         break
        #
        # if not overlap:
        #     circles.append((x,y))

    #print("------------------- end while overlapping-------------")
    return  circles
def boundry(x, h , r):
    if x < r:
        x=r
    if x > h-r :
        x=h-r
    return x
#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
counter =0

#--------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------
def strToFloat(x):
    x = x.replace(',', '.')
    return float(x)

# ------------------------------ loading plan ---------------------
class loading(QMainWindow):
    def __init__(self):
        super(loading, self).__init__()
        uic.loadUi("loading.ui", self)
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        # self.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        self.movie = QMovie('loading.gif')
        self.label.setMovie(self.movie)
    def start(self):
        self.show()
        self.movie.start()
    def stop(self):
        self.movie.stop()
        self.close()
#----------------------------------------first window -----------------------------------------
class splash(QMainWindow):
    def __init__(self):
        super(splash, self).__init__()
        uic.loadUi("splash.ui",self)
        self.setWindowFlag(QtCore.Qt.FramelessWindowHint)
        self.setAttribute(QtCore.Qt.WA_TranslucentBackground)
        ## QTIMER ==> START
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.progress)
        # TIMER IN MILLISECONDS
        self.timer.start(35)
        self.show()
        ## ==> END ##
    ## ==> APP FUNCTIONS
    ########################################################################
    def progress(self):
        global counter
        # SET VALUE TO PROGRESS BAR
        self.progressBar.setValue(counter)
        # CLOSE SPLASH SCREE AND OPEN APP
        if counter > 100:
            # STOP TIMER
            self.timer.stop()
            # SHOW MAIN WINDOW
            window2.show()
           # CLOSE SPLASH SCREEN
            self.close()
        # INCREASE COUNTER
        counter += 1


#---------------------- 2nd window
class secondWindow(QDialog):
    def __init__(self):
        super(secondWindow, self).__init__()
        uic.loadUi("secondwindow.ui",self)
        #----- CONDITION SUR LES CHAMPS
        self.calculer.clicked.connect(self.calculerfct)
        self.longitude.setValidator(QtGui.QDoubleValidator())
        self.latitude.setValidator(QDoubleValidator())
        self.rayon.setValidator(QDoubleValidator())
        self.nbr_population.setValidator(QIntValidator())
        self.nbr_generation.setValidator(QIntValidator())
        self.alpha.setValidator(QDoubleValidator())
        self.beta.setValidator(QDoubleValidator())
        self.gamma.setValidator(QDoubleValidator())

        self.max_iter =0
        self.population = 0
        self.long =0.0
        self.larg =0.0
        self.radius = 0.0
        self.nbrsite = 0
        self.tempfit = []
        self.individu = []
        self.fitness = 0.0
        self.puissance = 0.0
        self.active = 0





    def calculerfct(self):
        print("clicked")



        # ----------- RECUPERER LES VALEURS

        if len(self.longitude.text()) > 0 and len(self.latitude.text()) > 0 and len(self.rayon.text()) > 0 \
                and len(self.alpha.text()) >0 and len(self.beta.text()) >0  and len(self.gamma.text()) >0\
                and len(self.nbr_generation.text()) >0 and len(self.nbr_population.text()) > 0 :
            self.max_iter = int(self.nbr_generation.text())
            self.population = int(self.nbr_population.text())
            self.long = strToFloat(self.longitude.text())
            self.larg = strToFloat(self.latitude.text())
            self.radius = strToFloat(self.rayon.text())
            alpha = round(strToFloat(self.alpha.text()),2)
            beta = round(strToFloat(self.beta.text()),2)
            lamda = round(strToFloat(self.gamma.text()),2)
            print(alpha,beta,lamda)
            somme = math.fsum([Decimal(alpha) , Decimal(beta) , Decimal(lamda)])
            print(somme)
            if somme!= 1.0:
                print(' somme!= 1.0')
                dlg1 = QMessageBox(self)
                dlg1.setStyleSheet("color : rgb(225 , 225 , 225)")
                dlg1.setWindowTitle("warning !!!")
                dlg1.setText(" alpha + beta + lamda= 1!")
                dlg1.exec()
            if self.population < 4  :
                print(' self.population < 4  ')
                dlg1 = QMessageBox(self)
                dlg1.setStyleSheet("color : rgb(225 , 225 , 225)")
                dlg1.setWindowTitle("warning !!!")
                dlg1.setText(" nombre de population doit etre > 3!")
                dlg1.exec()

            if self.population > 3 and somme == 1.0 :
                # # appl gwo
                self.nbrsite = nbr_site(calculer_Scov(self.radius), surface(self.long, self.larg))
                width = int(math.ceil(self.long))
                height = int(math.ceil( self.larg))
                best_sol, self.tempfit = gwo(self.max_iter,self.population, self.nbrsite,width , height, self.radius, alpha,beta,lamda)
                print("gwo")
                self.individu = best_sol.individu
                self.fitness = best_sol.fitness
                self.puissance = best_sol.puisst
                self.active = onlyactif(self.individu)
                self.close()
                window3.show()
                window3.draw()




        else :
            print("rempli tous les champs")
            dlg = QMessageBox(self)
            dlg.setStyleSheet("color : rgb(225 , 225 , 225)")
            dlg.setWindowTitle("warning !!!")
            dlg.setText("rempli tous les champs!")
            dlg.exec()


#---------------------- third window
class thirdWindow(QDialog):
    def __init__(self):
        super(thirdWindow, self).__init__()
        uic.loadUi("thirdwindow.ui",self)
        self.affichage.setStyleSheet("color : rgb(225 , 225 , 225)")
        # self.affichage.setText("Nombre d'ameleoration = " , window2.nbr_generation , " \n\n"
        #                        "Fitness = " , window2.fitness , " \n\n"
        #                        "Couverture = " , round(window2.fitness , 2)  , "%\n\n"
        #                        "energie de transmission = " , window2.puissance , "  watt\n\n"
        #                        "Nombre de gNode = " , window2.nbrsite , " \n\n"
        #                        "Nombre de gNode actifs =" , window2.active , " \n\n"
        #                        "Temps d'execution = " , window2.tempfit[len(window2.tempfit)-1][0] , "  \n\n")
        self.resultats.clicked.connect(self.Resultat)

        #------------------ affichage -----------------------
        self.planification.setStyleSheet("color : rgb(225 , 225 , 225)")


        #---------------------- GWo call--------------
        #---------------------------------------------------

    def draw(self):
        scene = QGraphicsScene()
        self.ellipse_color = QtGui.QBrush(QtGui.QColor(222, 74, 50), 1)
        self.rectangle_color = QtGui.QBrush(QtGui.QColor("#ffffff") )
        self.pen = QPen(Qt.black , 0.2)
        self.planification.setScene(scene)
        #-------------- couverture %
        width = int(math.ceil(window2.long))
        height = int(math.ceil(window2.larg))
        couverture = verif_overlapp(window2.individu,width,height,window2.radius) *100
        temps_exec = window2.tempfit[len(window2.tempfit)-1][0]

        self.progressBar.setValue(int(couverture))
        n = 20
        rayon = window2.radius *2
        lst = window2.individu
        rect = scene.addRect(0, 0, window2.long * n, window2.larg * n , self.pen , self.rectangle_color )
        cpt=0
        for i in range(len(lst)):
            if(lst[i][2] == 1):
                cpt += 1
                x=lst[i][0] -rayon/2
                y =lst[i][1] - rayon/2
                if x <0 :
                    x = 0
                if y <  0 :
                    y=0
                ellipse = scene.addEllipse(x*n,y*n , rayon *n , rayon*n , self.pen, self.ellipse_color )


            # compteur = 0
            # while compteur< 100:
            #     compteur += 1
            #     self.progressBar.setValue(compteur)

        #----------------- pourcetage-------------
        #-------------- affichage
        #
        self.affichage.setText(
            'Population size  = ' + str(window2.population) + '\n\n'+ \
            'The number of iterations   = ' + str(window2.max_iter) +'\n\n'+\
            'Fitness = ' + str(window2.fitness) +  '\n''\n'+ \
            'Total transmission power = ' + str(window2.puissance) + '\n''\n' + \
            'Coverage = ' +  str(couverture)+ '%' + '\n''\n'+\
            'The total number of gNodeBs = '+  str(window2.nbrsite)+'\n''\n'+ \
            'The number of active gNodeBs = ' + str(cpt)+'\n''\n'+\
            'Execution time = ' +       str(temps_exec)  )


    def Resultat(self):
        self.close()
        window4.displayinfo()
        window4.courbe()
        window4.show()

#----------------------------------------- 4th window
class forthWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(forthWindow, self).__init__()
        self.ui = uic.loadUi("FORTHWINDOW.ui",self)

    def displayinfo(self):
        lst = window2.individu
        tab = [0.0] * len(lst)
        for i in range(len(lst)):
            tab[i]= {"N° gNodeB": i ,"X": lst[i][0] , "Y":lst[i][1] , "Etat":lst[i][2] , "Energie":lst[i][3]}


        print(tab)
        row = 0
        self.positions.setRowCount(len(tab))
        for i in tab:
            self.positions.setItem(row ,0, QtWidgets.QTableWidgetItem(str(i["N° gNodeB"])))
            self.positions.setItem(row ,1, QtWidgets.QTableWidgetItem(str(i["X"])))
            self.positions.setItem(row ,2, QtWidgets.QTableWidgetItem(str(i["Y"])))
            self.positions.setItem(row ,3, QtWidgets.QTableWidgetItem(str(i["Etat"])))
            self.positions.setItem(row ,4, QtWidgets.QTableWidgetItem(str(i["Energie"])))
            row = row+1
        """" --------------------------------------"""


    def courbe(self):

        lis = window2.tempfit
        print(lis)
        temp =[0.0]* len(lis)
        fitness = [0.0]* len(lis)

        for i in range(len(lis)):
            temp[i] = int(lis[i][0])
            fitness[i] = lis[i][1]

        fig, ax = plt.subplots()
        ax.plot(temp, fitness)
        ax.grid()
        matplotlib.pyplot.ylabel('fitness',fontsize=8)
        matplotlib.pyplot.xlabel('time (s)',fontsize=8)
        matplotlib.pyplot.xticks(fontsize=6)
        matplotlib.pyplot.yticks(fontsize=6)
        fig.savefig("test.png")
        fig.set_size_inches(4,3.3)
        canvas = FigureCanvasQTAgg(fig)
        canvas.setParent(self.courbefitness)

#-------------------------------
#------------------ main ------------------

app = QtWidgets.QApplication(sys.argv)
load = loading()
sp = splash()
window2 = secondWindow()
window3 = thirdWindow()
window4 = forthWindow()
try:
    sys.exit(app.exec_())
except:
    print("Exiting")
