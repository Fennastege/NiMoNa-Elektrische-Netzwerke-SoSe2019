import matplotlib.pyplot as plt

def f(p1=0, p2=0, e1=0, e2=0, g1=0, g2=0):
    return (p1*(e1-g1*p2), -p2*(e2-g2*p1))

def rk4(e1=0, e2=0, g1=0, g2=0, h=0, p1s=0, p2s=0, w=0, ts=0):
    t=[ts]
    p1=[p1s]
    p2=[p2s]
    if (w<0 or p2s<0 or p1s<0):
        print("w, p1s und p2s müssen positiv sein.")
    else:
        i=0
        while i < (w/h):
            t.append(h+i*h)
            k1, l1 =f(p1[i], p2[i], e1, e2, g1, g2)
            k2, l2 =f(p1[i]+h/2*k1, p2[i]+h/2*l1, e1, e2, g1, g2)
            k3, l3 =f(p1[i]+h/2*k2, p2[i]+h/2*l2, e1, e2, g1, g2)
            k4, l4 =f(p1[i]+h*k3, p2[i]+h*l3, e1, e2, g1, g2)
            p1.append(p1[i]+h/6*(k1+2*k2+2*k3+k4))
            p2.append(p2[i]+h/6*(l1+2*l2+2*l3+l4))
            i=i+1
        plt.ylabel('p(t)')
        plt.xlabel('t')
        plt.title("rk4")
        plt.plot(t, p1,'g', label='p1')
        plt.plot(t, p2,'r', label='p2')
        plt.legend()
        plt.grid(True)
        plt.show()
        plt.ylabel('p2')
        plt.xlabel('p1')
        plt.title("rk4")
        plt.plot(p1, p2,'g')
        plt.grid(True)
        plt.show()

rk4(2.0, 0.8, 0.02, 0.0002, 0.025, 1000, 20, 50, 0)