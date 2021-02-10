from fenics import *
from mshr import *
import numpy as np 
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random 
'''
The unit of length in this calculation is in [micro-meter]!!
'''

def fenicsSimulator(T=10,
                    num_steps=5,
                    lengthOfBox=1, ## in nanometer
                    InputCoords='test',
                    CalcOutput='uoutput',
                    PrevCalc=None,
                    stimulation=1,
                    autocrine=1
                   ):
   
    # time step & some constant values 
    dt = T / num_steps # time step size
    tol = 1E-14
    rofInjury=1 #<---- This is given by experiment (2 um) 
    
    # ----------------------------------- #
    #   Creation of Domain for PDE Calc   #
    # ----------------------------------- #
    
    # Create mesh and define function space
    domain = Rectangle(Point(0,0), Point(lengthOfBox,lengthOfBox))
    
    # importing the coordinate from yaml file 
    with open(InputCoords+'.yml') as file:
        coords = yaml.load(file,Loader=yaml.FullLoader)
    
    # creating the domain of cells for Dirichlet boundary condition for autocrinic release of ATP 
    NoOfCells = int(coords['NoCell'])
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~ Autocrinic Mechanism Related ~~~~#----------------------------
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if autocrine != 0:
        for i in np.arange(NoOfCells):
            xc = coords[str(i)][0]*lengthOfBox
            yc = coords[str(i)][1]*lengthOfBox
            marker = coords[str(i)][2]
            
            if marker == 'resting':
                rofCell = 0.01
            elif marker == 'activated':
                rofCell = 1
                
            circle = Circle(Point(xc,yc),rofCell)
            domain -= circle
    
    # switching between Dirichlet boundary condition and neumann boundary condition 
    # this is introduce to include pulsetile stimulation at the injury site 
    if stimulation != 0:
        xl = lengthOfBox/2
        yl = lengthOfBox/2
        circle = Circle(Point(xl,yl),rofInjury)
        domain -= circle
        
    mesh = generate_mesh(domain,100)
    V = FunctionSpace(mesh, 'P', 1)
    
    # ------------------------------------------------------ #
    # Neumann Boundary Condition requires creating subdomain #
    # ------------------------------------------------------ #
    ''' 
    This step allows to create subdoamin around cells in order to estimate 
    the local concentration of substances around cell, which is utilized to 
    initiate autocrinic release of ATP. 
    '''
    
    markerDomain = int(999)
    materials = MeshFunction("size_t", mesh, 2)
    materials.set_all(markerDomain)
        
    ## Note regarding marker
    '''
    Fenics seems not like to take marker including multiple of 5... it's weird
    '''
    
    Add = 111
    markers = []
    
    for i in np.arange(NoOfCells):
        xc = coords[str(i)][0]*lengthOfBox
        yc = coords[str(i)][1]*lengthOfBox
    
        class Omega_1(SubDomain):
            def inside(self, x, on_boundary):
                return sqrt(pow((x[0]-xc),2)+pow((x[1]-yc),2)) <= 10 - tol
        
        subdomain_1 = Omega_1()
        
        if i%4 == 0 and i > 0:
            Add = Add + 6
        
        k = i + Add
        markerNo = int(k)
        #print(markerNo)
        subdomain_1.mark(materials,markerNo)
        markers.append(markerNo)
    
    plot(materials)
    plt.savefig('subdomain.png')
    
    ds = Measure('ds',domain=mesh, subdomain_data=materials)
    dx = Measure('dx',domain=mesh, subdomain_data=materials)

    ##################################################################################################
    
    # ---------------------------------------- #
    # Configuration of simulation box boundary #
    # ---------------------------------------- #
    
    # Define boundary condition
    ## step1: boundary of simulation box 
    '''
    We are assuming that the edge of simulation box is infinitely long or in experiment, it is being flushed out.
    Thus, we are setting Dirichlet boundary condition with constant 0.0 at the edge to flux them out
    '''
    
    # edge concentration 
    u_bc = Constant(0.0)
    
    # boundary condition statement     
    def boundaryXL(x, on_boundary):
        return near(x[0], 0*lengthOfBox, tol) 
    def boundaryXH(x, on_boundary):
        return near(x[0], 1*lengthOfBox, tol)
    def boundaryYL(x, on_boundary):
        return near(x[1], 0*lengthOfBox, tol)
    def boundaryYH(x, on_boundary):
        return near(x[1], 1*lengthOfBox, tol)
    
    # setting up the boundary by DirichletBC
    bc_xl = DirichletBC(V, u_bc, boundaryXL)
    bc_xh = DirichletBC(V, u_bc, boundaryXH)
    bc_yl = DirichletBC(V, u_bc, boundaryYL)
    bc_yh = DirichletBC(V, u_bc, boundaryYH)
    bcs = [bc_xl, bc_xh, bc_yl, bc_yh]  
            
    ## Step 2: boundary of injury site
    '''
    The experiment used 1 mM of ATP at the tip of pippette (injury site).
    According to the conversion from mM to mol per cubic um, we applied 1e-18 as a conversion factor 
    '''
    concATP = 1e-18
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~ Autocrinic Mechanism Related ~~~~#--------------------------
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if stimulation != 0:
        u_C = Constant(concATP)
        bc = 'on_boundary && sqrt((x[0]-'+str(xl)+')*(x[0]-'+str(xl)+') + (x[1]-'+str(yl)+')*(x[1]-'+str(yl)+'))<'+str(rofInjury+tol)
        bc_circle = DirichletBC(V, u_C, bc)
        bcs.append(bc_circle)
    
    # --------------------------------------------------- #
    # Setting up the calculation details in terms of math #
    # --------------------------------------------------- #
    
    #Define initial value
    '''
    PrevCalc indicates there is previously calculated initial condition
    '''
    
    if PrevCalc is None:
        u_0 = Constant(0.0)
    else:
        mesh_s = Mesh('mesh.xml')        
        u0File = HDF5File(MPI.comm_world,PrevCalc+'.h5','r')
        Vo = FunctionSpace(mesh_s,'P',1)
        u_0 = Function(Vo)
        u0File.read(u_0,'/u')
        u0File.close()
        u_0.set_allow_extrapolation(True)
        
    u_n = project(u_0, V)
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(0.0)
    ##################################################################################################
    # Diffusion coefficient of ATP obtained from Pete's BPS paper
    D1 = Constant(145)
       
    F = dt*dot(D1*grad(u), grad(v))*dx(markerDomain) 
    F += (u-u_n)*v*dx 
    
    
    for i in np.arange(NoOfCells):
        markerNo = int(markers[i])
        F += dt*dot(D1*grad(u),grad(v))*dx(markerNo) 
    
    #    g1*v*ds(1) - g2*v*ds(2) - g3*v*ds(3) 
    a, L = lhs(F), rhs(F)
    ##################################################################################################
   
    # Time-stepping
    u = Function(V)
    t = 0
    
    count = 0
    tol2 = 1e-21 #6.5e-20 <------------ This is constant to make sure each site of cell won't consume the ATP due to lower concententration than its surroundings
    LocalConc = {}
    
    # ---------------------------------------------------------------------------------------------------------------------------
    for n in range(num_steps):

        # Update current time
        t += dt
        
        if n > 0:
            ##########################################################################################
            ## Dirichlet boundary condition that reflects some amount of substance contained in cells.
            u.set_allow_extrapolation(True)
            m = 1 
            j = 0
            for i in np.arange(NoOfCells):
                xc = coords[str(i)][0]*lengthOfBox
                yc = coords[str(i)][1]*lengthOfBox
                marker = coords[str(i)][2]
                
                markerNo = int(markers[i])
                area = assemble(Constant(1)*dx(markerNo))
                
                if area <= 1e-21:
                    conc = concATP*0.55
                else:
                    conc = assemble(u*dx(markerNo))/area

                if conc < tol2:
                    conc = tol2
                    
                if n == num_steps-1:
                    LocalConc[str(i)] = conc
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                #~~~~ Autocrinic Mechanism Related ~~~~#--------------------------
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#    
                if autocrine != 0:
                    
                    if marker == 'resting':
                        concFactor = 0.01
                    elif marker == 'activated':
                        concFactor = 0.1
                        
                    dummy = (concATP*concFactor)/(1+((0.2*concATP)/conc)**1.2) # <---- 1.2 original 0.65
                    u_Ca = Constant(dummy)
                    bc_cell = 'on_boundary && sqrt(pow((x[0]-'+str(xc)+'),2)+pow((x[1]-'+str(yc)+'),2))<'+str(rofCell+tol)
                    bcCell = DirichletBC(V, u_Ca, bc_cell)
                    bcs.append(bcCell)
                ##########################################################################################
            
        # Compute solution
        solve(a == L, u, bcs)

        # Update previous solution
        u_n.assign(u) 
        
        plot(u, vmin=1e-40, vmax=1e-18)
    # ---------------------------------------------------------------------------------------------------------------------------
    
    # -------------------- #
    # Store images or data #
    # -------------------- #
    # ---------------------------------------------------------------------------------------------------------------------------
    for i in np.arange(NoOfCells):
        x = coords[str(i)][0]*lengthOfBox
        y = coords[str(i)][1]*lengthOfBox
        marker = coords[str(i)][2]
        
        if marker == 'resting':
            plt.plot(x,y,'kx')
        elif marker == 'activated':
            plt.plot(x,y,'ko')
    # ---------------------------------------------------------------------------------------------------------------------------
    
    plt.xlabel('x [um]')
    plt.ylabel('y [um]')
    plt.savefig(CalcOutput+'.png')

    ## Store the calculation outcome in h5 file format 
    uFile = HDF5File(MPI.comm_world,CalcOutput+'.h5','w')
    uFile.write(u,'/u')
    uFile.close()
    ## Store mesh 
    File("mesh.xml") << mesh
    
    with open('LocalConcOfCell.yml','w') as file:
        document = yaml.dump(LocalConc,file)
    
'''
The next part that I need to do is to store the very last calculation in h5 file format 
and call it back to concentration total 
and call it back when I calculate the concentration gradient. 
'''
        
#
# MAIN routine executed when launching this script from command line
#
if __name__ == "__main__":
    import sys

    T = 10
    num_steps = 5
    lengthOfBox = 1
    InputCoords = 'test'
    CalcOutput = 'uoutput'
    PrevCalc = None
    autocrine = 1
    stimulation = 1
    
    for i,arg in enumerate(sys.argv):
    # calls 'runParams' with the next argument following the argument '-validation'
        if arg=="-T":
            T = np.float(sys.argv[i+1])
        
        if arg=="-steps":
            num_steps = np.int(sys.argv[i+1])
                       
        if arg=="-CalcOut":
            CalcOutput = sys.argv[i+1]
            
        if arg=="-lengthOfBox":
            lengthOfBox = np.float(sys.argv[i+1])
        
        if arg=="-PrevCalc":
            PrevCalc = sys.argv[i+1]
            
        if arg=="-InputCoords":
            InputCoords = sys.argv[i+1]
            
        if arg=="-Auto":
            autocrine = np.int(sys.argv[i+1])
            
        if arg=='-stim':
            stimulation = np.int(sys.argv[i+1])
         
fenicsSimulator(T=T,
                num_steps=num_steps,
                lengthOfBox=lengthOfBox,
                InputCoords=InputCoords,
                CalcOutput=CalcOutput,
                PrevCalc=PrevCalc,
                stimulation=stimulation,
                autocrine=autocrine,
               )
