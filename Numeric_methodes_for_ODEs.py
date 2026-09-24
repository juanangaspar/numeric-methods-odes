"""

@author: Juan Antonio Gaspar Pascual


The goal of this project is to implement some numerical methods that will allow us to
approximate the solution of a differential equation. To compare them, we will use the error
of each one with respect to the exact solution, which is known in this case, and in the last
experiment we will compare how the orders of the methods change depending on the meshes.

"""
from pylab import *





#1. The differential equation we will be trying to approximate is the following:

def fun(t,y):
    return -y + exp(-t)*cos(t) #that is, the ODE given by y'= fun(t,y)


#Since this is a "toy" problem, we know the exact solution of the ODE. It is:
def exacta(t):
    return exp(-t)*sin(t)



#2. Now let's move on to the numerical methods we will use:

def adams_bashforth2(a,b,fun, N,y0):
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a)/float(N) 
    y[0] = y0
    f[0] = fun(t[0],y[0])    
    y[1] = y[0] + h*f[0] # First iteration with Euler's method
    t[1] = t[0]+h
    f[1] = fun(t[1], y[1])
    for k in range(1,N):
        y[k+1] = y[k]+0.5*h*(3.0*f[k] - f[k-1])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
        
    return (t,y)



def adams_bashforth3(a,b,fun, N,y0):
     
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a)/float(N) 
    y[0] = y0
    f[0] = fun(t[0],y[0])
    for k in range(2):
        y[k+1] = y[k] + h*fun(t[k]+0.5*h, y[k]+0.5*h*f[k]) # First two iterations with the midpoint method.
        t[k+1] = t[k]+h
        f[k+1] = fun(t[k+1], y[k+1])
    for k in range(2,N):
        y[k+1] = y[k]+h/12*(23*f[k] - 16*f[k-1] + 5*f[k-2])
        t[k+1] = t[k] + h
        f[k+1] = fun(t[k+1], y[k+1])
         
    return (t,y)

def adams_moulton3(a,b,fun,N,y0): 
     
    y = zeros(N+1)
    t = zeros(N+1)
    f = zeros(N+1)
    t[0] = a
    h = (b-a)/float(N) 
    y[0] = y0
    f[0] = fun(t[0],y[0])
    for k in range(2): # We start with RK4
        K1 = f[k]
        K2 = fun(t[k]+0.5*h, y[k]+0.5*h*K1)
        K3 = fun(t[k]+0.5*h, y[k]+0.5*h*K2)
        K4 = fun(t[k]+h, y[k]+h*K3)
        y[k+1] = y[k] + h/6*(K1+2*K2+2*K3+K4) 
        t[k+1] = t[k]+h
        f[k+1] = fun(t[k+1], y[k+1])
    maxiter = 0
    for k in range(2,N):
        C = y[k] + h/24 * (19*f[k] -5*f[k-1] + f[k-2])
        t[k+1] = t[k]+h
        
        z_0 = y[k] #seed of the fixed-point algorithm
        
        for l in range(200):
            z_sig = h*9/24*fun(t[k+1],z_0) + C #To solve the equation required at each step, we have used this fixed-point algorithm
            if abs(z_sig-z_0) <1e-12: #this is the error-based stopping criterion
                iteraciones = l+1 
                if iteraciones > maxiter:
                    maxiter = iteraciones
                break
            z_0 = z_sig
        else:
            print('El algoritmo del punto fijo no converge')
        
        y[k+1] = z_sig
        f[k+1] = fun(t[k+1],y[k+1])
         
    return (t,y)





#3. Experiments

    #3.1 The first test will be an experiment with a fixed mesh N = 40,
    #that is, the interval where we want to approximate the ODE will be divided into 40 equal parts
def experimento_mallado_fijo40():
    #data
    a = 0
    b = 5
    y0 = 0 #this will be the initial value we take
    N = 40 
    
    
    #Numerical approximations
    (t_ab2,y_ab2) = adams_bashforth2(a, b, fun, N, y0)
    (t_ab3,y_ab3) = adams_bashforth3(a, b, fun, N, y0)
    (t_am3, y_am3) = adams_moulton3(a, b, fun, N, y0)
    
    y_exacta = exacta(t_ab2)
    
    
    #We compute and display the maximum errors of each method:
    error_ab2 = max(abs(y_ab2-y_exacta))
    error_ab3 = max(abs(y_ab3-y_exacta))
    error_am3 = max(abs(y_am3-y_exacta))
    
    print('Los errores para N=40 son:')
    print('Error en AB2:', error_ab2)
    print('Error en AB3:', error_ab3)
    print('Error en AM3:', error_am3)
    
    
    #Now we will plot the approximations on the same graph, together with the exact solution
    figure('Comparación métodos para N = 40')
    
    plot(t_ab2,y_exacta,'k-', label = 'Solución exacta')
    plot(t_ab2,y_ab2, 'o--', label = 'AB2')
    plot(t_ab3,y_ab3, 's--', label = 'AB3')
    plot(t_am3,y_am3, '^--', label = 'AM3')
    
    xlabel('t')
    ylabel('y(t)')
    legend()
    title('Comparación métodos para N = 40')
    grid(True)
    show()

experimento_mallado_fijo40()

print()
print()
print()
print()

    #3.2 The next experiment will consist of comparing the errors of these methods
    #by changing the mesh size. We will also make a plot for each
    #method showing it visually:

def experimento_varios_mallados(metodo, nombre_metodo):
    #data
    a = 0
    b = 5
    y0 = 0 #this will be the initial value we take
    mallados = [10,20,40,80,160,320]
    
    print('Errores para', nombre_metodo)
    
    figure('Comparación variando N - ' + nombre_metodo)
    for N in mallados:
        (t_metodo, y_metodo) = metodo(a,b,fun,N,y0)
        y_exacta = exacta(t_metodo)
        
        error = max(abs(y_exacta-y_metodo))
        
        print('Para N =', N, 'el error es', error)
    
        #Now we will build the plot where this can be seen visually:
        plot(t_metodo,y_metodo, 'o--', label = 'N='+str(N))
    
    #we plot the exact solution too
    t_exacta = linspace(a,b,1000)
    y_exacta = exacta(t_exacta)
    plot(t_exacta,y_exacta,'k-', label = 'exacta')
    
    xlabel('t')
    ylabel('y(t)')
    title('Aproximaciones de ' + nombre_metodo +' para varios mallados')
    legend()
    grid(True)
    show()

    #now we call each method:
experimento_varios_mallados(adams_bashforth2, "AB2")
experimento_varios_mallados(adams_bashforth3, "AB3")
experimento_varios_mallados(adams_moulton3, "AM3")

print()
print()
print()
print()


    #3.3 Now for the last experiment, which will show us the orders of convergence
    #computed empirically and approximately as we vary N
def experimento_ordenes_de_convergencia(metodo,nombre_metodo):
    #data
    a = 0
    b = 5
    y0 = 0 #this will be the initial value we take
    mallados = [10,20,40,80,160,320]
    
    print('Los órdenes de convergencia aproximado del método', nombre_metodo,'son:')
    print('El orden del primer mallado N = 10 no lo podemos calcular dado que necesitamos otro error para compararlo')
    
    errores=[]
    for N in mallados:
        (t_metodo,y_metodo) = metodo(a,b,fun,N,y0)
        y_exacta = exacta(t_metodo)
        
        
        error = max(abs(y_exacta-y_metodo))
        errores.append(error)
        
        
    for i in range(1,len(mallados)):
        orden = log(errores[i-1]/errores[i])/log(2)
        print('Para N=', mallados[i],'el orden aproximado es:',orden)


experimento_ordenes_de_convergencia(adams_bashforth2, "AB2")
print()
experimento_ordenes_de_convergencia(adams_bashforth3, "AB3")
print()
experimento_ordenes_de_convergencia(adams_moulton3, "AM3") 