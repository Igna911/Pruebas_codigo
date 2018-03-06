####Optimizador 2 etapas####

from math import e,exp,radians, cos, pi, sin, degrees, atan, log10

'''
Este programa servirá para establecer la trayectoria a describir por el
avión.  Contiene la maniobra de giro, seguida por otra de ascenso.  Al final se
procederá a exportar las características y cómo varían a lo largo de la
trayectoria para poder obtener sus gráficas.
'''

'''Condiciones g0itatorias y constantes atmosféricas'''
G = 6.673e-11  # Constante de g0itación.
MT = 5.972e24  # Masa terrestre.
MU = G * MT
RT = 6378136.3  # Radio terrestre.
R_AIR = 287  # Constante de los gases ideales.
g0 = 9.80665  # Aceleración gravitatoria inicial
RHO_SL = 101325 / (R_AIR * 288.15)  # Densidad del aire a nivel del mar.
GAMMA = 1.4  # Coeficiente de dilatación adiabática.

##Condiciones de viscosidad
beta_visc = 0.000001458
S_visc = 110.4

'''
-----------------------ATMÓSFERA ESTÁNDAR INTERNACIONAL-----------------------
Esta subrutina nos permitirá obtener los valores de presión, temperatura y
densidad del aire en funcióin de la altura. Están sacados de la ISA.  Más
adelante, en los siguientes bucles, se llamará a las funciones de P, T y rho
que variarán con respecto a h.
'''

#Estas son las alturas estipuladas según la normativa.

H_ISA1 = 11000
H_ISA2 = 20000
H_ISA3 = 32000
H_ISA4 = 47000
H_ISA5 = 51000
H_ISA6 = 71000
H_ISA7 = 84852

  

#Ahora se programan las variables termodinámicas, en función de la altura,
# y se relacionarán con los valores de T y alfa para cada altura estipulada.

def temperature(alt):
    '''Cálculo de la temperatura en función de la altura dada por el modelo ISA
    '''
    if alt < H_ISA1:
        h_0 = 0
        t_0 = 288.15
        alfa_isa = -.0065
    elif alt < H_ISA2:
        h_0 = H_ISA1
        t_0 = 216.65
        alfa_isa = 0
    elif alt < H_ISA3:
        h_0 = H_ISA2
        t_0 = 216.65
        alfa_isa = .001
    elif alt < H_ISA4:
        h_0 = H_ISA3
        t_0 = 228.65
        alfa_isa = .0028
    elif alt < H_ISA5:
        h_0 = H_ISA4
        t_0 = 270.65
        alfa_isa = 0
    elif alt < H_ISA6:
        h_0 = H_ISA5
        t_0 = 270.65
        alfa_isa = -.0028
    elif alt < H_ISA7:
        h_0 = H_ISA6
        t_0 = 214.65
        alfa_isa = -.002
    else:
        h_0 = H_ISA7
        t_0 = 214.65 - .002 * (H_ISA7 - H_ISA6)
        alfa_isa = 0
    return t_0 + alfa_isa * (alt - h_0)

def density(alt):
    '''Cálculo de la densidad en función de la altura dada por el modelo ISA
    '''
    rho0 = RHO_SL
    t_isa = temperature(alt)
    h_0 = 0
    t_0 = 288.15
    alfa_isa = -.0065
    if alt < H_ISA1:
    
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA1) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA1
    if alt < H_ISA2:
        
        return rho0 * exp(-g0 * (alt - h_0) / (R_AIR * t_isa))
    rho0 = rho0 * exp(-g0 * (H_ISA2 - h_0) / (R_AIR * temperature(H_ISA2)))
    h_0 = H_ISA2
    t_0 = 216.65
    alfa_isa = .001
    if alt < H_ISA3:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA3) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA3
    t_0 = 228.65
    alfa_isa = .0028
    if alt < H_ISA4:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA4) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA4
    if alt < H_ISA5:
        
        return rho0 * exp(-g0 * (alt - h_0)/(R_AIR * t_isa))
    rho0 = rho0 * exp(-g0 * (H_ISA5 - h_0)/(R_AIR * temperature(H_ISA5)))
    h_0 = H_ISA5
    t_0 = 270.65
    alfa_isa = -.0028
    if alt < H_ISA6:
        
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA6) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA6
    t_0 = 214.65
    alfa_isa = -.002
    if alt < H_ISA7:
       
        return rho0 * (t_isa / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    rho0 = rho0 * (temperature(H_ISA7) / t_0)**(-g0 / (R_AIR * alfa_isa) - 1)
    h_0 = H_ISA7
    return rho0 * exp(-g0 * (alt - h_0) / (R_AIR * t_isa))

def pressure(alt):
    '''Cálculo de la presión en función de la altura dada por el modelo ISA
    '''
    return density(alt) * R_AIR * temperature(alt)

def viscosity(alt):	
    '''Cálculo de la viscosidad en función de la altura dada por el modelo ISA
    '''							
    if alt < H_ISA1:
	    h_0 = 0
	    t_0 = 288.15
	    alfa_isa = -0.0065
	    t = t_0 + alfa_isa * (alt- h_0)
     
    elif alt< H_ISA2:
      
	    h_0 = H_ISA1
	    t_0 = 216.65   
	    alfa_isa = 0
	    t = t_0 + alfa_isa * (alt- h_0)  
		 
    elif alt< H_ISA3:
   
	    h_0 = H_ISA2
	    t_0 = 216.65
	    alfa_isa = 0.001
		   
	    t = t_0 + alfa_isa * (alt- h_0) 
                    
    elif alt< H_ISA4:
          
        h_0 = H_ISA3
        t_0 = 228.65                          
        alfa_isa = 0.0028
     
        t = t_0 + alfa_isa * (alt- h_0)                           
                          
    elif alt< H_ISA5:
        
        h_0 = H_ISA4
        t_0 = 270.65                                                               
        alfa_isa = 0                                 
        t = t_0 + alfa_isa * (alt- h_0)                                  
    
    elif alt< H_ISA6:
       
        h_0 = H_ISA5
        t_0 = 270.65                                     
        alfa_isa = -0.0028
    
        t = t_0 + alfa_isa * (alt- h_0)                                         
    
    elif alt< H_ISA7:
    
        h_0 = H_ISA6
        t_0 = 214.65                                               
        alfa_isa = -0.002
                                               
        t = t_0 + alfa_isa * (alt- h_0)                                         
        
    elif alt> H_ISA7: 
                                                                                             
	    h_0 = H_ISA7
	    t_0 = 214.65 - 0.002 * (H_ISA7 - H_ISA6)													                                                    
	    alfa_isa = 0
	    t = t_0 + alfa_isa * (alt- h_0)                                                    
		
    return (beta_visc * t ** (3 / 2)) / (t + S_visc)   

'''
-------------------CARACTERÍSTICAS GEOMÉTRICAS DEL MISIL-------------------
'''

diametro_m=0.5           #diametro misil
longitud_cono = 0.9      #longitud cono misil
longitud_misil = 3       #longitud total misil


Sw_aleta = 0.07875                           #Superficie de una aleta del AIM (tomado como ref.)
Craiz_aleta = 0.24                           #cuerda raiz de aleta
Cmedia_aleta = 0.18                          #cuerda media de aleta
espesor_aleta = 0.0065                       #espesor de aleta
tao_aleta = espesor_aleta/Cmedia_aleta       #tao aleta
num_aletas = 4                               #numero aletas
 
Swtotal_aletas=Sw_aleta*num_aletas           #Superficie total de aletas
Sref_aletas = Swtotal_aletas/2               #Superficie de referencia aletas


Sup_cono=pi*(diametro_m/2)*(longitud_cono**2+(diametro_m/2)**2)**(1/2)
Sup_total=2*pi*(diametro_m/2)*(longitud_misil-longitud_cono)
Sref_misil=pi*(diametro_m/2)**2
Sgases=pi*((diametro_m)*0.9/2)**2 #Area salida de los gases (consideramos el área de salida de la tobera)
Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 

angulo_cono = atan(0.5*diametro_m/longitud_cono)*(180/pi)  #angulo de la generatriz del cono

'''CDO del misil'''

def Cdll(Ml):

    def coef_resistencia_base_misil(Ml):
    
    	if Ml < 0.8:
    		x0 = 0
    		x1 = 0
    		x2 = 0
    		x3 = 0
    		x4 = 0
    	elif Ml < 1:
    		x0 = -1.548523
    		x1 = 6.05972764
    		x2 = -7.30548391
    		x3 = 2.96129532
    		x4 = 0
    	elif Ml < 1.1:
    		x0 = 5.79090984*10**3
    		x1 = -2.19843314*10**4
    		x2 = 3.12774812*10**4
    		x3 = -1.97644892*10**4
    		x4 = 4.68059822*10**3
    	elif Ml < 1.5:
    		x0 = -4.11856506
    		x1 = 1.42267421*10**1
    		x2 = -1.69678524*10**1
    		x3 = 8.771665
    		x4 = -1.67398037
    	elif Ml < 2.2:
    		x0 = 3.0748*10**-1
    		x1 = -1.3258*10**-1
    		x2 = 2.8812*10**-2
    		x3 = 0
    		x4 = 0
    	elif Ml <=3.5:                            
    		x0 = 1.8481*10**-1
    		x1 = -2.2895*10**-2
    		x2 = 5.1876*10**-3
    		x3 = -4.0742*10**-4
    		x4 = 0
    	elif Ml >3.5:                            
    		x0 = 0.15
    		x1 = 0
    		x2 = 0
    		x3 = 0
    		x4 = 0
    		    		
    	return x4*Ml**4 + x3*Ml**3 + x2*Ml**2 + x1*Ml + x0
     
    CD_base_misil = coef_resistencia_base_misil(Ml)*Ratio_areas
    
  
    #CALCULO DEL COEFICIENTE DE FRICCION
    ##COEFICIENTE DE FRICCIÓN DEL CONO
    ###CALCULO DEL REYNOLDS1 
    
    Re_cono=rho*vl*longitud_cono/Mu_Visc #REYNOLDS 2
    
    def cfcono_misil(Re_cono):
    		
        #LAMINAR
    	if Re_cono < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    	    cfi_cono=0.664*Re_cono**(-1/2)
    		
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    	    cf_cono=2*cfi_cono
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    	    cfm_cono=cf_cono*(1/(1+0.17*Ml**2))**0.1295
    
    			#CALCULO COEFICIENTE DE FRICCION DEL CONO
    	    
    		
    		#TURBULENTO
    	else:	
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cono=0.288*((log10(Re_cono))**(-2.45))
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf_cono=cfi_cono*1.597*((log10(Re_cono))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfm_cono=cf_cono*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			#CALCULO COEFICIENTE DE FRICCION DEL CONO
    
    		
    		
    	return cfm_cono*Sup_cono/Sref_misil
    
    
    ##COEFICIENTE DE FRICCIÓN DEL CILINDRO
    
    
    Re_cil=rho*vl*(longitud_misil-longitud_cono)/Mu_Visc  #REYNOLDS 2
    
    
    def cfcil(Re_cil):
    	
    		
    		#LAMINAR
    	if Re_cil < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cil=0.664*Re_cil**(-1/2)
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    		cf_cil=2*cfi_cil
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    		cfm_cil=cf_cil*(1/(1+0.17*Ml**2))**0.1295
    
    			
    		
    		
    		
                #TURBULENTO
    	else:			
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfi_cil=0.288*((log10(Re_cil))**(-2.45))
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf_cil=cfi_cil*1.597*((log10(Re_cil))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfm_cil=cf_cil*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			
    
    		
    		
    	return cfm_cil*Sup_total/Sref_misil
    
    
    #CALCULO DEL COEFICIENTE DE FRICCION TOTAL REFERIDO A LA SUPERFICIE TRANSVERSAL
    CDFriccion_cono = cfcono_misil(Re_cono)
    CDFriccion_cil = cfcil(Re_cil)
    CDFriccion=CDFriccion_cono+CDFriccion_cil
    
    #CALCULO DEL COEFICIENTE DE ONDA
    def cd_onda(Ml,angulo_cono):
    	if Ml >= 1 :
    		
    		return (0.083+0.096/(Ml**2))*(angulo_cono/10)**1.69
       #REGIMEN SUBSONICO
    
    	elif Ml <1 :
    		
    		return 	((60/((longitud_cono/diametro_m)**3))+0.0025*(longitud_cono/diametro_m))*CDFriccion
    
    CD_onda = cd_onda(Ml,angulo_cono)
    	
    	
    ####RESISTENCIA ALETAS
    ######################
    
    #COEFICIENTE DE ONDA
    def cd_onda_aletas(Ml,angulo_cono):
    	if Ml >= 1 :
    		
    		return ((4*tao_aleta**2)/((Ml**2-1)**0.5))*Swtotal_aletas/Sref_misil
       #REGIMEN SUBSONICO
    
    	elif Ml <1 :
    		
    		return 	0
    
    CD_onda_aletas = cd_onda_aletas(Ml,angulo_cono)
        
    Re_aletas=rho*vl*Craiz_aleta/Mu_Visc #REYNOLDS aletas
    
    #coeficiente fricción aletas
    def cf_aletas(Re_aletas):
    		
        #LAMINAR
    	if Re_aletas < 1e6 :
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    	    cfialetas=0.664*Re_aletas**(-1/2)
    		
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL MEDIO
    
    	    cf1aletas=2*cfialetas
    
    			#CALCULO COEFICIENTE DE FRICCION COMPRESIBLE
    
    	    cfmaletas=cf1aletas*(1/(1+0.17*Ml**2))**0.1295
    
    			
    	    
    		
    		#TURBULENTO
    	else:	
    			#CALCULO COEFICIENTE DE FRICCION LOCAL INCOMPRESIBLE
    
    		cfialetas=0.288*(log10(Re_aletas))**(-2.45)
    
    			#CALCULO COEFICIENTE DE FRICCION LOCAL COMPRESIBLE
    
    		cf1aletas=cfialetas*1.597*((log10(Re_aletas))**(-0.15))
    
    			#CALCULO COEFICIENTE DE FRICCION MEDIO
    
    		cfmaletas=cf1aletas*(1/(1+(GAMMA-1)/2*Ml**2)**0.467)
    
    			
    
    		
    		
    	return cfmaletas*Swtotal_aletas/Sref_misil
    
    
    CDFriccion_aletas = cf_aletas(Re_aletas)
    
    return CD_base_misil + CDFriccion + CD_onda + CD_onda_aletas + CDFriccion_aletas



'''Inicialización de variables y diferenciales para la maniobra del misil'''

vl = 288.875 #Inicialización de la velocidad
yl = 18586  #Inicialización de la altitud 
thetal = 88.83*pi/180 #Inicialización del ángulo de asiento
thetalgrados=thetal*(180/pi) #Conversión de radianes a grados del ángulo de asiento



vxl=vl*cos(thetal) #Inicialización de la componente horizontal de velocidad  
vyl=vl*sin(thetal) #Inicialización de la componente verical de velocidad 
tl = 0 #Inicialización temporal 
xl = 0 #Inicialización de la posición en el eje x 
sl = 0 #Inicialización del arco recorrido
dvxl=0 #Inicialización del diferencial de la componente horizontal de velocidad 
dvyl=0 #Inicialización del diferencial de la componente vertical de velocidad 
dsl=0 #Inicialización del diferncial del arco recorrido
dxl=0 #Inicialización del diferencial de la posición
dyl=0 #Inicialización del diferencial de la altitud
dtl=0.001 #Inicialización del diferencial de tiempo
dthetal = 0 #Inicialización del diferencial del ángulo de asiento                                

rho = density(yl)  # Densidad inicial del aire (kg/m3).
p = pressure(yl)  # Presión inicial del aire (Pa).
T = temperature(yl)  # Temperatura inicial del aire (K).
Mu_Visc = viscosity(yl) # Viscosidad 
#A la altura inicial el avión vuela en vuelo estacionario.   
dt = 0.1  # Diferencial de tiempo (s).

Ml=vl/((GAMMA*R_AIR*T)**0.5)
Cdl=Cdll(Ml)
D_misil=0.5*rho*Cdl*Sref_misil*vl**2

'''Caracterisitcas del cohete'''

gasto = 30
gasto_etapa2= gasto 
isp1 = 280
Empuje_misil = gasto * isp1 * 9.81

ratio_estructural = .05 #'Relacion estrutural, sera un 5% del peso del propulsante de etapa
masa_maniobras = 0
masa_util = 30 + masa_maniobras

'''Velocidades implicadas'''
v_orb = 7800
v_rotacional = 400

#ESTE SERÁ EL PARAMETOR A ITERAR#
v_loss = 782
v_loss_empirica = 0
error = 800 
contador = 1

'''Calculos de la optimizacion'''
while error > .1 and contador <= 50: #calcula la diferencia netre la v_loss claculada y la anterior

    '''Inicialización de variables y diferenciales para la maniobra del misil'''
    
    vl = 288.875 #Inicialización de la velocidad
    yl = 18586.2  #Inicialización de la altitud 
    thetal = 88.83*pi/180 #Inicialización del ángulo de asiento
    thetalgrados=thetal*(180/pi) #Conversión de radianes a grados del ángulo de asiento
    
    
    
    vxl=vl*cos(thetal) #Inicialización de la componente horizontal de velocidad  
    vyl=vl*sin(thetal) #Inicialización de la componente verical de velocidad 
    tl = 0 #Inicialización temporal 
    xl = 0 #Inicialización de la posición en el eje x 
    sl = 0 #Inicialización del arco recorrido
    dvxl=0 #Inicialización del diferencial de la componente horizontal de velocidad 
    dvyl=0 #Inicialización del diferencial de la componente vertical de velocidad 
    dsl=0 #Inicialización del diferncial del arco recorrido
    dxl=0 #Inicialización del diferencial de la posición
    dyl=0 #Inicialización del diferencial de la altitud
    dtl=0.01 #Inicialización del diferencial de tiempo
    dthetal = 0 #Inicialización del diferencial del ángulo de asiento                                
    
    rho = density(yl)  # Densidad inicial del aire (kg/m3).
    p = pressure(yl)  # Presión inicial del aire (Pa).
    T = temperature(yl)  # Temperatura inicial del aire (K).
    Mu_Visc = viscosity(yl) # Viscosidad 
    #A la altura inicial el avión vuela en vuelo estacionario.   
    dt = 0.1  # Diferencial de tiempo (s).
    
    Ml=vl/((GAMMA*R_AIR*T)**0.5)
    Cdl=Cdll(Ml)
    D_misil=0.5*rho*Cdl*Sref_misil*vl**2
    
    r = RT + yl  # Distancia al centro de la Tierra (m).
    g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
    
    print(contador, v_loss)
    velocidad_ideal = v_orb - v_rotacional + v_loss - vl
    print(contador, velocidad_ideal)
    f = ((1 + ratio_estructural)/exp(velocidad_ideal / (2 * isp1 * g0))) - ratio_estructural
    
    masa_total = masa_util / (f**2)
    masa_etapa2 = masa_util / f
    masa_propulsante_etapa1 = (masa_total - masa_etapa2) / (1 + ratio_estructural)
    masa_propulsante_etapa2 = (masa_etapa2 - masa_util) / (1 + ratio_estructural)
    masa_estructura1 = masa_propulsante_etapa1 * ratio_estructural
    masa_estructura2 = masa_propulsante_etapa2 * ratio_estructural
    
    t_combustion= masa_propulsante_etapa1/gasto
    t_combustion_2_etapa=masa_propulsante_etapa2/gasto_etapa2 
    retardo_encendido = 0
    t_fin_combustion2= t_combustion + t_combustion_2_etapa + retardo_encendido
    
    masa_misil = masa_total #esta sera la masa que ira variando
    v_loss_empirica = 0
    while tl<= t_combustion and thetal>0:
            
        '''Combustion de la primera etapa
        '''

        
        tl=tl+dtl #Evolución temporal (s)
        xl=xl+dxl #Posición horizontal (m)
        yl=yl+dyl #Altitud (m) 
        sl=sl+dsl #Arco recorrido (m)
        
        thetal=thetal+dthetal #Ángulo de asiento
        thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
        
        
        vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
        vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
        vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
        
        r = RT + yl  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        
        rho = density(yl)  # Densidad (kg/m3).
        T = temperature(yl)  # Temperatura (K).
        Mu_Visc = viscosity(yl) # Viscosidad 
        
        Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
        
        
        Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 
                        
              
        Cdl=Cdll(Ml)                     
        D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
        Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
        Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
        
        dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
        dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)

        
        dvxl=dvxl_perdidas+Empuje_misil*cos(thetal)*dtl/masa_misil
        dvyl=dvyl_perdidas+Empuje_misil*sin(thetal)*dtl/masa_misil
        masa_misil=masa_misil-gasto*dtl
            
        dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
        dxl=vxl*dtl #Diferencial de la posición (m)
        dyl=vyl*dtl #Diferencial de la altitud (m)
        dsl=vl*dtl #Diferencial del arco recorrido (m)
        
        v_loss_empirica = v_loss_empirica + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5
        
    masa_misil  = masa_misil - masa_estructura1
    while t_fin_combustion2>=tl>t_combustion and thetal>0 and yl<500000:
    
        ''' Combustion de la segunda etapa'''

        
        tl=tl+dtl #Evolución temporal (s)
        xl=xl+dxl #Posición horizontal (m)
        yl=yl+dyl #Altitud (m) 
        sl=sl+dsl #Arco recorrido (m)
        
        thetal=thetal+dthetal #Ángulo de asiento
        thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
        
        
        vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
        vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
        vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
        
        r = RT + yl  # Distancia al centro de la Tierra (m).
        g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
        
        rho = density(yl)  # Densidad (kg/m3).
        T = temperature(yl)  # Temperatura (K).
        Mu_Visc = viscosity(yl) # Viscosidad 
        
        Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
        
        
        Ratio_areas=(Sref_misil-Sgases)/Sref_misil #Relación de áreas 
                        
              
        Cdl=Cdll(Ml)                     
        D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
        Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
        Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
        
        dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
        dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)
            
        dvxl=dvxl_perdidas+Empuje_misil*cos(thetal)*dtl/masa_misil
        dvyl=dvyl_perdidas+Empuje_misil*sin(thetal)*dtl/masa_misil
        masa_misil=masa_misil-gasto_etapa2*dtl

        dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
        dxl=vxl*dtl #Diferencial de la posición (m)
        dyl=vyl*dtl #Diferencial de la altitud (m)
        dsl=vl*dtl #Diferencial del arco recorrido (m)
        
        v_loss_empirica = v_loss_empirica + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5
        
    masa_misil = masa_misil - masa_estructura2
    while thetal>0 and yl<500000:
            '''Fin de etapas
            '''
    
            
            tl=tl+dtl #Evolución temporal (s)
            xl=xl+dxl #Posición horizontal (m)
            yl=yl+dyl #Altitud (m) 
            sl=sl+dsl #Arco recorrido (m)
            
            thetal=thetal+dthetal #Ángulo de asiento
            thetalgrados=thetal*(180/pi) #Conversión de radianes a grados
            
            
            vxl=vxl+dvxl #Componente horizontal de la velocidad (m/s)
            vyl=vyl+dvyl #Componente vertical de la velocidad (m/s)
            vl=(vxl**2+vyl**2)**0.5 #Módulo de la velocidad (m/s)
            
            r = RT + yl  # Distancia al centro de la Tierra (m).
            g0 = MU / r**2  # Aceleración de la gravedad (m/s2).
            
            rho = density(yl)  # Densidad (kg/m3).
            T = temperature(yl)  # Temperatura (K).
            Mu_Visc = viscosity(yl) # Viscosidad 
            
            Ml=vl/((GAMMA*R_AIR*T)**0.5) #Mach de vuelo
            
            
            Ratio_areas=1 #Relación de áreas 
                            
                  
            Cdl=Cdll(Ml)                     
            D_misil=0.5*rho*Cdl*Sref_misil*vl**2 #Fuerza de resistencia (N)
            Dx=D_misil*cos(thetal) #Componente horizontal de la fuerza de resistencia (N)
            Dy=D_misil*sin(thetal) #Componente vertical de la fuerza de resistencia (N)
            
            dvxl_perdidas=-(Dx/masa_misil)*dtl #Diferencial de la componente horizontal de la velocidad (m/s)
            dvyl_perdidas=-g0*dtl-(Dy/masa_misil)*dtl #Diferencial de la componente vertical de la velocidad (m/s)
            
            dvxl=dvxl_perdidas
            dvyl=dvyl_perdidas
            
            dthetal=dtl*g0*cos(thetal)/(-vl) #Diferencial del ángulo de asiento
            dxl=vxl*dtl #Diferencial de la posición (m)
            dyl=vyl*dtl #Diferencial de la altitud (m)
            dsl=vl*dtl #Diferencial del arco recorrido (m)
            
            v_loss_empirica = v_loss_empirica + (dvxl_perdidas**2 + dvyl_perdidas**2)**0.5 
    print(v_loss_empirica)
    error = abs(v_loss_empirica - v_loss)
    v_loss = v_loss_empirica
    

    
    if contador == 50:
        print ('error de interacion en el ángulo ', thetalgrados, 'con un error de ', error,' m/s')
    contador = contador + 1