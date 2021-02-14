import numpy as np 
from numpy import sin, cos
from math import pi 
from manimlib.imports import * 

g = 5
dt = 0.1
class MorinProblem(Scene):
    def construct(self):
        axes = NumberPlane().set_opacity(0.2)
        self.add(axes) 

        Conditions = TexMobject(r"""
        m_1 = m_2 \\ 
        \theta_0 = \frac{\pi}{6}\\
        g = 5
        """).to_edge(UP+LEFT)
        self.add(Conditions) 
        def y(t_1):
            r = np.zeros((4), float)
            k1 = np.zeros((4), float) 
            k2 = np.zeros((4), float)
            k3 = np.zeros((4), float)
            k4 = np.zeros((4), float)
            t = 0.
            r[0] = 1. 
            r[1] = 0. 
            r[2] = pi/6. 
            r[3] = 0. 
            def f(r): 
                h = np.zeros((4), float) 
                h[0] = r[1] 
                h[1] = (r[0]*r[3]**2. - g*(1. - cos(r[2])))/2.
                h[2] = r[3]
                h[3] = -2.*r[1]*r[3]/r[0]-g*sin(r[2])/r[0]
                return h
        
    
        
            def rk4():
                k1 = dt*f(r)
                k2 = dt*f(r+k1/2.)
                k3 = dt*f(r+k2/2.)
                k4 = dt*f(r+k3)
                return (k1+2.*(k2+k3)+k4)/6.

        
            while t < t_1 : 
                n = int(t_1/dt)
                t = t + dt 
                r = r + rk4()
            return r 
        def eta(t_1):
            return y(t_1)[0]
        def v(t_1):
            return y(t_1)[1]
        def o(t_1):
            return y(t_1)[2]
        def omega(t_1):
            return y(t_1)[3]
    

        t1 = np.arange(0, 20, dt)
        n=len(t1)
        f1 = np.vectorize(eta)
        f2 = np.vectorize(o)
        l = np.array([6]*(n))
        x1 = f1(t1)*sin(f2(t1))
        y1 = f1(t1)*cos(f2(t1))
        y2 = l - f1(t1) 

        Pulley_1 = Dot(radius = 0.25).move_to(2*RIGHT+2*UP).set_color(RED)
        Pulley_2 = Dot(radius = 0.25).move_to(2*LEFT+2*UP).set_color(RED)
        Body_1 = Dot(radius =0.1).move_to((x1[0]+2)*RIGHT + (y1[0]-2)*DOWN).set_color(BLUE)
        Body_2 = Dot(radius =0.1).move_to(2*LEFT + (y2[0]-2)*DOWN).set_color(BLUE)
    
        Line1 = self.getline(Pulley_1,Pulley_2) 

        Line2 = self.getline(Pulley_1,Body_1)
        Line2.add_updater(
                        lambda mob: mob.become(
                       self.getline(Pulley_1,Body_1)
                        ))

        Line3 = self.getline(Pulley_2,Body_2)
        Line3.add_updater(
                        lambda mob: mob.become(
                        self.getline(Pulley_2,Body_2)
                        ))

        self.add(Line1,Line2,Line3,Pulley_1,Pulley_2,Body_1,Body_2)
        
        traj = VGroup()

        for i in range(len(x1)-1):
            a = 30
            self.remove(traj)
            traj = VGroup()
            
            if i >= 30:
                for n in range(i,i-30,-1):
                    line1 = Line((2+x1[n])*RIGHT + (y1[n]-2)*DOWN , (2+x1[n-1])*RIGHT + (y1[n-1]-2)*DOWN).set_stroke(YELLOW,width=2).set_opacity(0.03*a)
                    traj.add(line1)
                    a -= 1
            
            else:
                for n in range(i,0,-1):
                    line1 = Line((2+x1[n])*RIGHT + (y1[n]-2)*DOWN,(2+x1[n-1])*RIGHT + (y1[n-1]-2)*DOWN).set_stroke(YELLOW,width=2).set_opacity(0.05*a)
                    traj.add(line1)
                    a -= 1
            
            self.add(traj)
            self.remove(Body_1,Body_2)
            self.add(Body_1,Body_2)
            self.play(
                Body_1.move_to,(2+ x1[i+1])*RIGHT + (y1[i+1]-2)*DOWN,
                Body_2.move_to,2*LEFT + (y2[i+1]-2)*DOWN,
                run_time=0.1,rate_func=linear)
        

     

    
    def getline(self,Point1,Point2):
            start_point = Point1.get_center()
            end_point = Point2.get_center() 
            line = Line(start_point,end_point).set_stroke(width=2) 
            return line



        

        
    

    


            

        





