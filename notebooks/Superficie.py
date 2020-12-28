import sympy as sy

class surface:
    def __init__(self,sup,u,v):
        
        self.sup=sup
        self.u=u
        self.v=v
        self.sdu=[sy.diff(sup[i],u) for i in range(len(sup))]
        self.sdv=[sy.diff(sup[i],v) for i in range(len(sup))]
        self.sduXsdv = [sy.Matrix(self.sdu).cross(sy.Matrix(self.sdv))[i] for i in range(len(self.sup))] 
        norm = sy.Matrix(self.sduXsdv).norm()
        self.n = [sy.simplify(self.sduXsdv[i]/norm) for i in range(len(self.sduXsdv))]
        print(self.n)
        self.sdudu = [sy.diff(self.sdu[i],u) for i in range(len(self.sdu))]
        self.sdudv =  [sy.diff(self.sdu[i],v) for i in range(len(self.sdu))] 
        self.sdvdv = [sy.diff(self.sdv[i],v) for i in range(len(self.sdv))]

        self.firstFF = self.firstFF()
        self.secondFF = self.secondFF()
        self.Weingarten = self.Weingarten()
        self.curvatureGauss = self.curvatureGauss()
        self.curvatureMedia = self.curvatureMedia()
        self.christoffel1 = self.christoffel1()
        self.christoffel2 = self.christoffel2()
        

    def firstFF(self):
        e = sy.Matrix(self.sdu).dot(self.sdu)
        f = sy.Matrix(self.sdu).dot(self.sdv)
        g = sy.Matrix(self.sdv).dot(self.sdv)

        firstFF=sy.simplify(sy.Matrix([[e,f],[f,g]]))
        return firstFF
    def secondFF(self):
        E = sy.Matrix([self.sdu,self.sdv,self.sdudu]).det()
        F = sy.Matrix([self.sdu,self.sdv,self.sdudv]).det()
        G = sy.Matrix([self.sdu,self.sdv,self.sdvdv]).det()
        secondFF=sy.simplify(sy.Matrix([[E,F],[F,G]])/sy.Matrix(self.sduXsdv).norm())
        return secondFF
    def Weingarten(self):
        invfirstFF = self.firstFF.inv()
        
        #return sy.simplify(sy.Matrix(self.firstFF().inv()*(self.secondFF)))
        return sy.simplify(invfirstFF*self.secondFF)

    def curvatureGauss(self):
        return sy.simplify(self.Weingarten.det())
    def curvatureMedia(self):
        return sy.simplify(1/2*sy.trace(self.Weingarten))
    def christoffel1(self):
        E = self.firstFF[0,0]
        F = self.firstFF[0,1]
        G = self.firstFF[1,1]

        Edu = sy.diff(E,self.u)
        Edv = sy.diff(E,self.v)
        Fdu = sy.diff(F,self.u)
        Fdv = sy.diff(F,self.v)
        Gdu = sy.diff(G,self.u)
        Gdv = sy.diff(G,self.v)

        return sy.simplify(1/2*sy.Matrix([[Edu,Edv,2*Fdv-Gdu],[2*Fdu-Edv,Gdu,Gdv]]))
    def christoffel2(self):
        return sy.simplify(self.firstFF.inv()*self.christoffel1)
    def mainCurvaturesAndDirections(self):
       
        self.curvatureMain = [sy.simplify(self.curvatureMedia-sy.sqrt(self.curvatureMedia**2-self.curvatureGauss)),sy.simplify(self.curvatureMedia+sy.sqrt(self.curvatureMedia**2-self.curvatureGauss))]
        self.directionsMain = self.Weingarten.eigenvects()
        
        return self.curvatureMain, self.directionsMain