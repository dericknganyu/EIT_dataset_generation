#If you wanna listening something I recommend: https://www.youtube.com/watch?v=lVFt__nrRxY
#Bachianas Brasileiras Nº 4 - Villa Lobos.
import matplotlib.pyplot as plt
from dolfin import*
from mshr import*
import numpy as np
from fenics import*

def MyMesh(r, n, n_vertex):
    """Function that generate the mesh based in number of vertex in boundary. It is extremely important that
    forward mesh and inverse mesh have vertex in commum, so this routine generate both, mesh_direct is a refiniment from mesh_inverse.
    
    :param r: Circle Radius
    :type r: float.
    :param n: Refinament parameter
    :type n: int.
    :param n_vertex: Vertices number in boundary
    :type n_vertex: int
    
    :return: :class:`dolfin.cpp.mesh.Mesh`

    :Example:
            >>> mesh_inverse, mesh_forward=MyMesh(r=1, n=8, n_vertex=100)
            
    .. image:: codes/mesh.png
      :scale: 75 %

    """
    ###Points Generator###
    points=[] #To save points
    for theta in np.linspace(0, 2*pi, n_vertex):
        a=[Point((cos(theta)*r,sin(theta)*r))]     #Defining Point object from Fenics.
        points=np.concatenate((points,a), axis=0)  #Grouping points.


    domain = Polygon(points)         #Function creates a polygon with the points.
    mesh = generate_mesh(domain,n)   #Get the polygon and generate the mesh

    mesh2 = refine(mesh)   #Get the mesh_inverse and refine it.
    mesh.radius=r  #Just add info. in the object.
    mesh2.radius=r #Just add info. in the object.

    return mesh, mesh2 #Mesh inverse and mesh forward.

def getBoundaryVertex(mesh,u):
    """ Functions that calculate the values of a function on boundary and return a array with that.
    
    :param mesh: Mesh where u is defined.
    :type mesh: :class:`dolfin.cpp.mesh.Mesh`
    :param u: Function that you want compute vertex values on the boundary.
    :type u: :class:`dolfin.cpp.mesh.Function`
    
    :return: array

    :Example:
    
    >>> u_boundary=getBoundaryVertex(mesh,u)
    
    """
    bmesh=BoundaryMesh(mesh, 'exterior', order=True) 
    indice_borda=bmesh.entity_map(0).array() #Get vertex indice (Boundary)
    
    u_bvertex=u.compute_vertex_values()[indice_borda] #Save all u(vertex).
    return u_bvertex

def getBoundaryVertexTwoMesh(mesh_inverse, mesh_forward, u, u0):#Fina, grossa, função.
    """ Functions that calculate the values of two function on the boundary and select the vertex in commum, then return an array with it.
    
    :param mesh_inverse: Coarsed mesh where Function u is defined.
    :type mesh_inverse: :class:`dolfin.cpp.mesh.Mesh`
    :param mesh_forward: Refined mesh where Function u0 is defined.
    :type mesh_forward: :class:`dolfin.cpp.mesh.Mesh`
    :param u: Function that you want compute vertex values on the boundary.`
    :type u: :class:`dolfin.cpp.mesh.Function`
    :param u0: Function that you want compute vertex values on the boundary.`
    :type u0: :class:`dolfin.cpp.mesh.Function`
    
    :return: (array) u_boundary, u0_boundary, vertex_index.

    :Example:
            >>> u_boundary, u0_boundary, vertex_index=getBoundaryVertexTwoMesh(mesh_inverse, mesh_direct, u, u0)
    """
    u_bvertex=[]    #Save u(vertex) mesh_inverse
    u_b2vertex=[]   #Save u(vertex) mesh_forward
    vertexnumber=[] #Save mesh index
    
  
    bmesh_inverse=BoundaryMesh(mesh_inverse, 'exterior', order=True)
    bmesh_forward=BoundaryMesh(mesh_forward, 'exterior', order=True)

    index_bond_inverse=bmesh_inverse.entity_map(0).array()
    index_bond_forward=bmesh_forward.entity_map(0).array()

    for ind_inv in index_bond_inverse:    #For each  index vertex in bmesh_inv.
        for ind_dir in index_bond_forward: #For each  index vertex in bmesh_forward.
            vertex1 = Vertex(mesh_inverse, ind_inv) #Map between index and vertex position (inverse).
            vertex2 = Vertex(mesh_forward, ind_dir) #Map between index and vertex position (forward).
            if vertex1.x(0)==vertex2.x(0) and vertex1.x(1)==vertex2.x(1): #If they have x_i=x_f and y_f=y_f.
                vertexnumber.append(ind_inv) #Save vertex number index for inverse.
                break

    for ind in vertexnumber:
        vertex = Vertex(mesh_inverse, ind) #Map between index and vertex position (inverse).
        u_bvertex.append(u(vertex.point())) #Append to ub_inverse
        u_b2vertex.append(u0(vertex.point())) #Append to ub_forward
                
    
    u_bvertex=np.array(u_bvertex)
    u_b2vertex=np.array(u_b2vertex)

    return u_bvertex, u_b2vertex, vertexnumber

def current_method(n_g, value=1, method=1):
    """This function an expression that represent the current in the vertex.
        
    :param n_g: Measurements number.
    :type n_g: int
    :param value: value in the vertex.
    :type value: float
    :param method: Current pattern.
    :type method: int
    :returns:  Expression -- Return list of expressions.
    
    Method Values:           
        1. 1 and -1 in opposite direction, where 50% of the boundary is always 0.
        
    :Example:

        .. code-block:: python

           "Current"
            n_g=2
            list_gs=current_method(n_g, value=1, method=1)

            for i in range(n_g):
                mesh=mesh_direct
                VD=FiniteElement('CG',mesh.ufl_cell(),1) 
                g_u=interpolate(list_gs[i], FunctionSpace(mesh,VD))
                g_u=getBoundaryVertex(mesh, g_u)
                bond=plot_boundary(mesh, data=g_u, name='boundary g'+str(i))

        .. image:: codes/current1.png
           :scale: 75 %       
           """
    #We need implement more methods here, please see the tutorial, there are several examples.
    h=pi/n_g/2 #Angular lenght of each "electrode"
    if method==1:
        #&& - and operator
        list_gs=[Expression(f" x[0]<cos(0+{h*i}) && x[0]>cos({h}+{h*i}) && x[1]>0 ? {value} : "+
                            f"(x[0]>cos(pi+{h*i}) && x[0]<cos(pi+{h}+{h*i}) && x[1]<0 ? -{value} : 0 )"
                            ,degree=1) for i in range(0,n_g*2,2)]
    return list_gs

def fn_addnoise(data, level, noise_type='uniform', seed=42):
    """Function receives a vector which represents the data in the electrodes and returns
    a noised vector with the chosen noise level and the type of noise.
    We use it in :func:`ForwardProblem.add_noise`.
    
    :param data: Vector with potencial in electrodes or any other vector.
    :type data: array
    :param level: Noise level (%), expect values between 0 and 1.
    :type level: float
    :param noise_type: Noise type, uniform or cauchy.
    :type method: str.
    :param seed: Seed for random function.
    :type seed: int.
    
    
    :returns:  Array -- Return noised vector.
    
    :Example:

    >>> print(np.ones(8))
    >>> print(fn_addnoise(data=np.ones(8), level=0.01, noise_type='cauchy', seed=32))
        [1. 1. 1. 1. 1. 1. 1. 1.]
        array([0.99905327, 1.02206251, 1.00356633, 1.00236212, 1.00101231, 0.99904405, 1.0105611 , 0.98656216])
    
    """
    i = len(data)
    # create 1D numpy data:
    npdata = np.asarray(data).reshape((i))
    delta=level*np.linalg.norm(npdata) #delta = noise_level * ||data||_2
    
    #add f normal noise:
    if noise_type=='uniform':
        np.random.seed(seed)                             #Set seed to generate random noise
        noise_f=np.random.randn(npdata.size)             #Generate random noise with vector size of data.
        noise_f=noise_f/np.linalg.norm(noise_f)*delta    #Normalize random noise and multiply by delta.
        noise = npdata + noise_f                         #Add noise to the vector.
    # add cauchy noise:
    elif noise_type=='cauchy':   
    #12 31 59
        np.random.seed(seed)                                     #Set seed to generate random noise
        noise_p = np.random.standard_cauchy(size=npdata.shape)   #Generate random noise with vector size of data.
        noise_p=noise_p/np.linalg.norm(noise_p)*delta            #Normalize random noise and multiply by delta
        noise = npdata + noise_p                                 #Add noise to the vector.
        
    return noise



def GammaCircle(mesh, in_v, out_v, radius,centerx, centery):
    """Function to create a circle in the mesh with some proprieties
    
        :param mesh: Mesh.
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        :param in_v: Value inside circle
        :type in_v: float
        :param out_v: Value outside circle
        :type out_v: float
        :param radius: Circle radius
        :type radius: float
        :param centerx: Circle center position x
        :type centerx: float
        :param centery: Circle center position y
        :type centery: float
       
        :returns:  Array -- Return a vector where each position correspond de value of the function in that element.
        
        :Example:

        >>> ValuesCells0=GammaCircle(mesh=mesh_direct, in_v=3.0, out_v=1.0, radius=0.50, centerx=0.25, centery=0.25)
        >>> print(ValuesCells0)
        [1. 1. 1. ... 1. 1. 1.]
        
        >>> "Plot"
        >>> gamma0=CellFunction(mesh_direct, values=ValuesCells0);   
        >>> V_DG=FiniteElement('DG',mesh_direct.ufl_cell(),0)
        >>> plot_figure(mesh_direct, V_DG, gamma0, name="Resposta gamma");
        
        .. image:: codes/gamma.png
           :scale: 75 %

       """
    
    ValuesGamma=np.zeros(mesh.num_cells()) #Null vector
    
    for i in range(0, mesh.num_cells()):
        cell = Cell(mesh, i) #Select cell with index i in the mesh.
        
        vertices=np.array(cell.get_vertex_coordinates()) #Vertex cordinate in the cell.
        x=(vertices[0]+vertices[2]+vertices[4])/3           
        y=(vertices[1]+vertices[3]+vertices[5])/3
        
        #If the baricenter is outside the circle...
        if ((x-centerx)**2+(y-centery)**2>=radius**2):
            ValuesGamma[i]=out_v
        else:
            ValuesGamma[i]=in_v
    
    return ValuesGamma


def GammaCircleRandom(mesh, n, in_v, out_v, radius,centerx, centery):

    N = np.ceil(np.random.uniform(0,1) * n+0.000001)
    Radius = np.zeros(N)
    Centerx = np.zeros(N)
    Centery = np.zeros(N)
    ValuesGamma=np.zeros(mesh.num_cells()) + out_v
    
    for k in range(N):
         r = np.random.uniform(radius[0],radius[1])
         Centerx = np.random.uniform(centerx[0],centerx[1])
         Centery = np.random.uniform(centery[0],centery[1])
         In_v = np.random.uniform(inv[0],inv[1])
         for i in range(0, mesh.num_cells()):
             cell = Cell(mesh, i) #Select cell with index i in the mesh.
             vertices=np.array(cell.get_vertex_coordinates()) #Vertex cordinate in the cell.
             x=(vertices[0]+vertices[2]+vertices[4])/3   
             y=(vertices[1]+vertices[3]+vertices[5])/3
             if ((x-Centerx)**2+(y-Centery)**2>=r**2):
                 ValuesGamma[i]=out_v
             else:
                 ValuesGamma[i]=In_v
    
    return ValuesGamma
  


class CellFunction(UserExpression):
    """Auxiliar function to transform an array to a Function
        We use it with :func:`GammaCircle()`
    
        :param mesh: Mesh.
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        :param values: Array with values of the function in the cell.
        :type values: array
        
        :Example:

        >>> ValuesCells0=np.zeros(mesh_inverse.num_cells())        #Define a vector of zeros
        >>> ValuesCells0[5]=1                                      #Cell 5 has value 1
        >>> gamma0=CellFunction(mesh_inverse, values=ValuesCells0);#Get vector and transform in a function cell

        If you want plot the function::
        
        >>> V_DG=FiniteElement('DG',mesh_inverse.ufl_cell(),0)     #Space of Finite Elemente descontinuous garlekin degree 0
        >>> Q=FunctionSpace(mesh_inverse,V_DG)                     #Functionspace to interpolate gamma
        >>> gamma0=interpolate(gamma0, Q)                          #Interpolation gamma to generate a function
        >>> p=plot(gamma0)                                         #plot gamma0
        >>> plot(mesh_inverse)                                     #plot mesh
        >>> plt.colorbar(p)                                        #set colorbar.

        .. image:: codes/cell_function_test.png
          :scale: 75 %

         
        """
    def __init__(self, mesh, values, **kwargs):
        self.mesh = mesh
        self.values=values #Values in a array           
        super().__init__(**kwargs)
    
    #Function that returns value in cell.
    def eval_cell(self, value, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        value[0]=self.values[cell.index()] #Set value of the cell using the array values. (Map between cell.index and values.)
        
    #"Just to erase a warning"            
    def value_shape(self):
        return ()

    
class ForwardProblem(object):
    """ Object Forward Problem EIT 2D Continous Model.
    
        :param mesh: Mesh. We recommend from :func:`MyMesh()`
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        
        :Example:
        
        .. code-block:: python
        
            "Basic Definitions"
            VD=FiniteElement('CG',mesh_direct.ufl_cell(),1) 
            F_Problem=ForwardProblem(mesh_direct)

            "Solver"
            list_u0=F_Problem.solve_forward(VD, gamma0, I_all)
            u0_boundary=F_Problem.boundary_array(mesh_inverse)
        
        If you need it, see :func:`GammaCircle` and :func:`CellFunction`.
        
    """
    
    def __init__(self, mesh):
        self.mesh=mesh #Mesh
        
    def solve_forward(self, V, gamma, I_all): 
        """ Solver Forward Problem EIT 2D
            
    :param V: FiniteElement Fenics object
    :param gamma: Finite Element Function
    :param I_all: Current density in each electrode for each measurement
    :type V: FiniteElement
    :type gamma: :func:`CellFunction()`
    :type I_all: :func:`current_method()` or list of arrays
            
    :returns:  (Array) -- Return function that is solution from variational problem.
        
     :Example:
     >>> F_Problem=ForwardProblem(mesh_direct)
     >>> list_u0=F_Problem.solve_forward(VD, gamma0, list_gs)
            
        """
        mesh=self.mesh  #Saving mesh
        n_g=len(I_all)  #Get number of experiments
        
        R=FiniteElement('R',mesh.ufl_cell(),0) #Constant for Lang. Mult

        W=V*R #Mixing Elements HxR
        W=FunctionSpace(mesh,W) #Defining function space.

        (u,c)=TrialFunctions(W) #Functions that will be reconstructed.
        (v,d)=TestFunction(W)   #Test functions

        lagrMult=(v*c+u*d)*ds #If that we have the ground potential. Integral_(dOmega) u ds =0

        a=inner(gamma*grad(u),grad(v))*dx+lagrMult # Integral( gamma*<grad_u,grad_v> ) dOmega + lagrMult
        A=assemble(a) #Make my matriz to solve Ax=b.
        #We only do it only one time, if have mult. measurements, we reuse it.
        
        #Construction my b vector. (Here b=L).
        sol_u=[] #save solution u
        for j in range(n_g): #For each current in I_all
            L=I_all[j]*v*ds #Intregral g_i*v*ds
            b = assemble(L) #Make my b vector
            
            w = Function(W) #Define a zero function based in W.
            U = w.vector()  #Return a vector. (x=U)
            solve(A, U, b) #Solve system AU=b. where A matrix and b vector.
            
            #Split w function in 2 parts, firt the function in H, and the constant lagr
            (u,c)=w.split(deepcopy=True)
            sol_u.append(u)
        
        self.sol_u=sol_u #Save in memory.
        
        return sol_u
    
    def boundary_array(self, mesh_inverse=None, concatenate=True):
        """ Get's the boundary values of function solution and returns array.
        If you set a coarse mesh you will get the values in the vertices that are commum.
        If you set conccatenate=False, will receive a array with separeted results, is usefull if you used more than one current.
        
        :param mesh: Corse Mesh. We recommend from :func:`MyMesh()`
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        :param concatenate: Default True
        :type concatenate: bool
        
        :returns:  (Array) -- Vertex values of the function.
        
        :Example:
        
        >>> u0_boundary=F_Problem.boundary_array(mesh_inverse)
        
        """
        sol_u=self.sol_u
        data=[]
    
        for i in range(len(sol_u)):
            if mesh_inverse==None:
                u_dados = getBoundaryVertex(self.mesh,self.sol_u[i])
            else:
                u_dados, u_dados, vertexnumber = getBoundaryVertexTwoMesh(mesh_inverse, self.mesh, sol_u[i], sol_u[i])
            
            if concatenate==False: data.append(u_dados)
            else : data=np.concatenate((data, u_dados))
                
        return data
    
        
    
    def plot_boundary(self, mesh_inverse=None, index=0 ):
        """ Get's the boundary values of function solution and returns a graph.
        If you set a coarse mesh you will get the values in the vertices that are commum and plot it.
        
        :param mesh: Corse Mesh. We recommend from :func:`MyMesh()`
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        :param index: Index of solution, if you need it.
        :type index: int
        
        :Example:
        
        >>> data_u0=F_Problem.plot_boundary(mesh_inverse, index=1)
        
        .. image:: codes/boundary_u.png
          :scale: 75 %
        
        """
        
        if mesh_inverse==None: mesh=self.mesh
        else: mesh=mesh_inverse
            
        u_data=self.boundary_array(mesh, concatenate=False)
        plot_boundary(mesh, data=u_data[index], name=f'boundary u_{index}', line=0)
        return
        
    
    def add_noise(self, noise_level=0, noise_type='uniform', seed=42, mesh=None):
        """ Function that add noise in the potential values.
        
            :param data: Vector with potencial in electrodes or any other vector.
            :type data: array
            :param level: Noise level (%), expect values between 0 and 1.
            :type level: float
            :param noise_type: Noise type, uniform or cauchy.
            :type noise_type: str.
            :param mesh: Corse Mesh. We recommend from :func:`MyMesh()`
            :type mesh: :class:`dolfin.cpp.mesh.Mesh`
    
            :returns:  Array -- Return a vector with potentials values concatenated.
            
            :Example:
            
            .. code-block:: python              
            
                "Noise Parameters"
                noise_level=0.01
                noise_type='uniform'
                seed=1
                u0_boundary=F_Problem.add_noise(noise_level noise_type, seed, mesh_inverse)

        """
        if mesh==None: mesh=self.mesh
        u_data=self.boundary_array(mesh, concatenate=False)
        
        #Add same noise in each experiment.
        vec_U=[]
        for i in range(len(u_data)): vec_U=np.concatenate((
            vec_U, fn_addnoise(u_data[i], noise_level, noise_type=noise_type, seed=seed)), axis=0)
        return vec_U

def etaFunc(mesh,cell_number): #Transforma a array com os valores na fronteira em uma função
    """Auxiliar function to do a characterisct function on an element. Value 1 on the element and 0 for the rest.
    
        :param mesh: Mesh.
        :type mesh: :class:`dolfin.cpp.mesh.Mesh`
        :param cell_number: Cell index
        :type cell_number: int
         
        """
    V=FiniteElement('DG',mesh.ufl_cell(),0) #Descontinous Galerkin degree 0
    Q=FunctionSpace(mesh,V) #Function Space based in DG and mesh
    eta=Function(Q)         #Null function 
    eta.vector()[cell_number]=1 #Eta is 1 in element with number=cell_number and zero otherwise.
    return eta

class InverseProblem(ForwardProblem):
    """Inverse Object EIT 2D 
    
    :param mesh: Any mesh from Fenics module. We recommend from :func:`MyMesh()`
    :type mesh: mesh
    :param V: FiniteElement Fenics object
    :type V: FiniteElement
    :param data: Vector with potencial in each vertex.
    :type data: array
    :param I_all: Current  for each measurement
    :type I_all: :func:`current_method()`
    
    :Example:
    
    .. code-block:: python
    
     VI=FiniteElement('CG',mesh_inverse.ufl_cell(),1) 
     InverseObject=InverseProblem(mesh_inverse, VI, u0_boundary, I_all)
     InverseObject.solve_inverse()
     gamma_k=InverseObject.gamma_k

    """
    
    def __init__(self, mesh, V, data, I_all):
        super().__init__(mesh)
        #"Basic definitions"
        self.V=FiniteElement('CG',mesh.ufl_cell(),1)  #Function Space CG degree 1 is necessary.
        self.I=I_all                   #Current pattern used in generated data.
        self.l=len(I_all)              #Number of measurements
        self.u0_boundary=data          #electrodes Potencial in array
        self.vertex=BoundaryMesh(mesh, 'exterior', order=True).entity_map(0).array()
        
        #"First guess and weight functions"
        self.Cellsgamma_k=np.ones(mesh.num_cells())*0.9           #First guess for Forwardproblem
        self.gamma_k=CellFunction(mesh, values=self.Cellsgamma_k) #Guess in cell function
        self.weight=np.ones(mesh.num_cells())        #Initial weight function
        self.innerstep_limit=50                      #Inner step limit while solve
        
        #"Solver configurations"
        self.weight_value=True  #Weight function in Jacobian matrix
        self.step_limit=30       #Step limit while solve
        self.min_v=0.05         #Minimal value in element for gamma_k
        
        #"Noise Configuration"
        self.noise_level=0      #Noise_level from data
        self.tau=1.01              #Tau for disprance principle
        
        #"Newton parameters"
        self.mu_i=0.9       #Mu initial (0,1]
        self.mumax=0.999    #Mu max
        self.nu=0.99        #Decrease last mu_n
        self.R=0.98         #Maximal decrease (%) for mu_n
        
        #"Inner parameters"
        self.inner_method='Landweber'  #Inner method for solve Newton
        
        self.land_a=1    #Step-size Landweber
        self.ME_reg=5E-4 #Regularization Minimal Error
        self.Tik_c0=1    #Regularization parameter Iterative Tikhonov
        self.Tik_q=0.95  #Regularization parameter Iterative Tikhonov
        self.LM_c0=1     #Regularization parameter Levenberg-Marquadt
        self.LM_q=0.95   #Regularization parameter Levenberg-Marquadt
        
        #"A priori information"
        self.gamma0=None #Exact Solution
        self.mesh0=None  #Mesh of exact solution
        
        #Creating a vector with all cell volumes. It's usefull for integrals in L2(Omega).
        cell_vec=[]
        for cell in cells(mesh):
            cell_vec.append(cell.volume())
        self.cell_vec=np.array(cell_vec)
        
        #Make a vector with boundary elements size that are chosen for the problem
        #This vector is used in normL2
        bmesh=BoundaryMesh(mesh, 'exterior', order=True) #Define boundarymesh
        bcell_vec=[]
        for bcell in cells(bmesh):  bcell_vec.append(bcell.volume()) #Save all boundary_elements size
        self.bcell_vec=np.tile(np.array(bcell_vec), self.l) #Adapting to l measurements.
        

        
        
    def solve_inverse(self):
        """Function that solves the inverse problem.
        
        :Example:
        
        >>> F_Problem=ForwardProblem(mesh_direct)
        >>> list_u0=F_Problem.solve_forward(VD, gamma0, I_all)

    """
        res_vec, error_vec=[], [] #To save about iterations
        self.innerstep_vec=[]     #Save inner_step newton
        mun_vec=[]                #Save mu in inner_step newton
        self.steps=0              #Save external step.

        ##############################################
        "First Forward solver"
        self.u = self.solve_forward(self.V, self.gamma_k, self.I)
        self.u_boundary=self.boundary_array() #Get boundary data and convert to array
        
        "First Save data"
        #Residue vector
        res_vec.append(np.linalg.norm(self.u0_boundary-self.u_boundary)/np.linalg.norm(self.u0_boundary)*100)
        self.innerstep_vec.append(int(0)) #Save number steps
        mun_vec.append(0)                 #Save number steps
        
        "Print information"
        if self.mesh0 is not None and self.gamma0 is not None:
            error_vec.append(self.error_gamma())    
            print("Error (%)=", error_vec[0], "Residuo (%)=", res_vec[0], " passo:", 0, "Inner step: ", 0)
        else:
            print("Residuo (%)=", res_vec[0], " passo:", 0, "Inner step: ", 0)
            
        ##############################################
        
        "Solver"
        #While discepancy or limit steps.
        while res_vec[self.steps]/100>=self.tau*self.noise_level and self.steps<=self.step_limit:

            "Derivative matrix calc"
            #If will be used LM or Tikhonov, we always have to calc. Jacobian.
            if (self.steps==0 and self.weight_value) or self.inner_method=='LM' or self.inner_method=='Tikhonov':
                Jacobiana_all=self.Jacobian_calc() #Derivative matrix calc
                
                #Create weight and add it.
                if self.weight_value: Jacobiana_all=self.weight_func(Jacobiana_all) 
                self.Jacobiana=Jacobiana_all
            else: Jacobiana_all=None

            "Inner iteration newton"
            sk, inner_step, mu=self.solve_innerNewton(Jacobiana_all)
            
            "Add sk in guess"
            self.Cellsgamma_k+=sk #Add a correction in each element
            
            #Don't have values less than c.
            self.Cellsgamma_k[self.Cellsgamma_k < self.min_v] = self.min_v
            self.gamma_k=CellFunction(self.mesh, values=self.Cellsgamma_k)  #Vector to function
            
            "Forward solver"
            self.u = self.solve_forward(self.V, self.gamma_k, self.I)
            self.u_boundary=self.boundary_array() #Get boundary data and convert to array
        
            "Saving data"
            res_vec.append(np.linalg.norm(self.u0_boundary-self.u_boundary)/np.linalg.norm(self.u0_boundary)*100)
            
            self.innerstep_vec.append(int(inner_step)) #Save number steps
            mun_vec.append(mu) #Save number steps
            
            self.steps+=1 #Next step
            
            if self.mesh0 is not None and self.gamma0 is not None:
                error_vec.append(self.error_gamma())
                print("Error (%)=", error_vec[self.steps], "Residuo (%)=", res_vec[self.steps],
                      " passo:", self.steps, "Inner step: ", inner_step)
            else: print("Residuo (%)=", res_vec[self.steps],
                      " passo:", self.steps, "Inner step: ", inner_step)
        
            #Vectors to memory object.
            self.res_vec=res_vec
            self.mun_vec=mun_vec
            self.error_vec=error_vec
        return
   
        
    def solve_innerNewton(self, Jacobiana_all):
        """Methods to solve inner step newton. Functions executed inside of :func:`solve_inverse()`. See set_InnerParameters() for more details.
            
    :param Jacobiana_all: Derivative Matrix generated by :func:`Jacobian_calc()`
    :type Jacobiana_all: Array ndim
    
    :returns:  (Array, int, float) -- Return a sk (Result of Inner Step to add in gamma_k), inner_step (Number of inner steps), mu (Regularization parameter used in the method).
            """
        b0=self.u0_boundary-self.u_boundary #Define vector b0 (Ask=b0)
        norm_b0=self.norm_funcL2(b0, 'dOmega')
        residuo=-b0       #Define res.
        norm_res=norm_b0  #Define norm_res first step.
        
        mu = self.newton_reg() #Calculate regularation parameter.
        inner_step=0
        
        sk=np.zeros(self.mesh.num_cells()) #s0 inicial do newton
        
        "------Landweber------"
        if self.inner_method=='Landweber':
            while norm_res>=mu*norm_b0 and inner_step<=self.innerstep_limit:
                sk+=-self.land_a*self.adj_dev_app_vec(residuo)

                residuo=-b0+self.dev_app_vec(sk)
                norm_res = self.norm_funcL2(residuo, 'dOmega')
                inner_step+=1
                #print(norm_res, mu*norm_b0, inner_step)
                
            
                
            "------Minimal Error------"
        elif self.inner_method=='ME':
            while norm_res>=mu*norm_b0 and inner_step<=self.innerstep_limit:
                sk_n=-self.adj_dev_app_vec(residuo)
                omega=(self.norm_funcL2(residuo, 'dOmega')**2/self.norm_funcL2(sk_n, 'Omega')**2)*self.ME_reg
                sk+=omega*sk_n
                
                residuo=-b0+self.dev_app_vec(sk)
                norm_res = self.norm_funcL2(residuo, 'dOmega')
                inner_step+=1
                #print(norm_res, mu*norm_b0, inner_step, omega)
        
            "------Conjugate-Gradient------"
        elif self.inner_method=='CG':
            while norm_res>=mu*norm_b0 and inner_step<=self.innerstep_limit:
                if inner_step==0:
                    rk=b0
                    ak=self.adj_dev_app_vec(rk)
                    pk=ak
                    ak_old=ak

                qk=self.dev_app_vec(pk)
                alphak=(self.norm_funcL2(ak_old, 'Omega')**2)/(self.norm_funcL2(qk, 'dOmega')**2)
                sk=sk+alphak*pk
                rk=rk-alphak*qk
                ak=self.adj_dev_app_vec(rk)
                betak=(self.norm_funcL2(ak, 'Omega')**2)/(self.norm_funcL2(ak_old, 'Omega')**2)
                pk=ak+betak*pk

                ak_old=ak
                residuo=-b0+self.dev_app_vec(sk)
                norm_res = self.norm_funcL2(residuo, 'dOmega')
                inner_step+=1
                #print(norm_res, mu*norm_b0, inner_step, alphak)
                #if alphak<=1: break

                
            "------Iterative Tikhonov------"
        elif self.inner_method=='Tikhonov' and inner_step<=self.innerstep_limit:
            while norm_res>=mu*norm_b0:
                alpha_k=self.Tik_c0*(self.Tik_q**inner_step)
                ADJ = Jacobiana_all.T
                square_m=ADJ@Jacobiana_all
                square_m+=alpha_k*np.identity(np.size(square_m, axis=0))
                sk=np.linalg.solve(square_m, ADJ.dot(b0)+alpha_k*sk)


                residuo=-b0+Jacobiana_all@sk
                norm_res = self.norm_funcL2(residuo, 'dOmega')
                inner_step+=1
                #print(norm_res, mu*norm_b0, inner_step)
                
            "------Levenberg-Marquadt------"
        elif self.inner_method=='LM' and inner_step<=self.innerstep_limit:
            while norm_res>=mu*norm_b0:
                alpha_k=self.LM_c0*(self.LM_q**inner_step)
                ADJ = Jacobiana_all.T
                square_m=ADJ@Jacobiana_all
                square_m+=alpha_k*np.identity(np.size(square_m, axis=0))
                sk=np.linalg.solve(square_m, ADJ.dot(b0))


                residuo=-b0+Jacobiana_all@sk
                norm_res = self.norm_funcL2(residuo, 'dOmega')
                inner_step+=1
                #print(norm_res, mu*norm_b0, inner_step)

        return sk, inner_step, mu

        
    def Jacobian_calc(self):
        """Calcuate derivative matrix. Function executed inside of :func:`solve_inverse()`.
        
        :returns:  (Array ndim) -- Return the derivative matrix.        
              """
        print("Calculando Jacobiana.")
        V=self.V
        mesh=self.mesh
        gamma=self.gamma_k
        R=FiniteElement('R',mesh.ufl_cell(),0)

        W=V*R
        W=FunctionSpace(mesh,W)

        bmesh=BoundaryMesh(mesh, 'exterior', order=True)
        indice_borda=bmesh.entity_map(0).array()
        linha=len(indice_borda) #Linha da matriz jacobiana
        
        coluna=mesh.num_cells() #Coluna da matriz jacobiana 
        

        
        (du,c)=TrialFunctions(W)
        (v,d)=TestFunction(W)
        lagrMult=(v*c+du*d)*ds
        a = gamma*inner(grad(du),grad(v))*dx + lagrMult
        A=assemble(a)
        
        for h in range(self.l):
            Jacobiana=np.zeros((linha,coluna)) #Armazena as derivadas em cada elemento    
            u=self.u[h]
            for cell in cells(mesh):    
                eta = etaFunc(mesh, cell.index())
                L=-eta*inner(grad(u),grad(v))*dx
                b = assemble(L)
                
                z = Function(W)
                U = z.vector()
                solve(A, U, b, 'cg', 'ilu')
                (du,c)=z.split(deepcopy=True)

                derivada=getBoundaryVertex(mesh,du)
                Jacobiana[:,cell.index()]=derivada
        
            Jacobiana=np.array(Jacobiana)    
            if h==0: Jacobiana_all=Jacobiana
            else: Jacobiana_all=np.concatenate((Jacobiana_all, Jacobiana), axis=0)

        print("Fim do cálculo da Jacobiana")
        
        return Jacobiana_all
    

    def adj_dev_app_vec(self, b0):
        mesh=self.mesh
        weight=self.weight
        vertex=self.vertex
        n_vert_boundary=len(vertex)
        gamma_k=self.gamma_k
        sol_u=self.sol_u

        ADJ=np.zeros(mesh.num_cells())

        for i in range(self.l):
            x_i, x_f = int(i*n_vert_boundary), int((i+1)*n_vert_boundary) #Para acertar as linhas
            sigma=self.Sigma(b0[x_i:x_f])
            psi=self.findpsi(gamma_k, sigma)
            ADJ+=self.calcADJ(psi, sol_u[i])
            
        ADJ=ADJ*1/weight 

        return ADJ

    def dev_app_vec(self, sk):
        mesh=self.mesh
        gamma_k=self.gamma_k
        V=self.V
        sol_u=self.sol_u
        vector_dev=sk

        R=FiniteElement('R',mesh.ufl_cell(),0)

        W=V*R
        W=FunctionSpace(mesh,W)

        bmesh=BoundaryMesh(mesh, 'exterior', order=True)
        indice_borda=bmesh.entity_map(0).array()
        derivada_sum=[]

        vector_dev=np.array(vector_dev)/np.array(self.weight)



        (du,c)=TrialFunctions(W)
        (v,d)=TestFunction(W)


        eta=CellFunction(mesh, values=vector_dev)
        lagrMult=(v*c+du*d)*ds(mesh)
        a = gamma_k*inner(grad(du),grad(v))*dx(mesh) + lagrMult
        A=assemble(a)
        
        for u in sol_u:    
            L=-eta*inner(grad(u),grad(v))*dx(mesh)

            b = assemble(L)
            w = Function(W)
            U = w.vector()
            solve(A, U, b, 'cg', 'ilu')
            (du_i,c)=w.split(deepcopy=True)

            derivada=getBoundaryVertex(mesh,du_i)
            derivada_sum = np.hstack((derivada_sum, derivada))

        return derivada_sum

        
        
    def weight_func(self, Jacobiana):
        """Determine the weights for Derivative Matrix and apply.
        
        :param Jacobiana: Derivative Matrix generated by :func:`Jacobian_calc()`
        :type Jacobiana: Array ndim
    
        :returns:  (Array ndim) -- Return the derivative matrix with weights.
        
        """
        self.weight=np.linalg.norm((Jacobiana.T*np.sqrt(self.bcell_vec)).T, axis=0)*(1/self.cell_vec)            
      
        Jacobiana=Jacobiana*1/self.weight
        return Jacobiana
    
    def newton_reg(self):
        """Determine the regulazation parameter for the Newton innerstep in :func:`solve_innerNewton()`."""
        #ref: https://doi.org/10.1137/040604029
        passo=self.steps
        innerstep_vec=self.innerstep_vec
        
        if passo>1:
            if innerstep_vec[passo-1]>=innerstep_vec[passo-2]: #If inner step increase
                mu_n=1-((innerstep_vec[passo-1]/innerstep_vec[passo])*(1-self.mu))
            else: mu_n=self.nu*self.mu #Remmember, nu<1.
            mu=self.mumax*np.amax([self.R*self.mu, mu_n]) #R limits how much the mu_n change.
        else: mu=self.mu_i #For the two first steps.
        
        self.mu=mu
        print("mu_n", mu)
        return mu

    def norm_funcL2(self, vector, domain):
        """This function is used in inner step of newton to calculate the norm's in L2 Omega and dOmega
         
        :param domain: 'Omega' or 'dOmega'.
        :type domain: str
        
        :returns:  (float) -- Return norm.
        
        """
        mesh=self.mesh
        weight=self.weight
        if domain=='Omega':
            norm=sum(weight*vector*vector*self.cell_vec)
        elif domain=="dOmega":
            norm=sum(vector*vector*self.bcell_vec)

        norm=sqrt(norm)
        return norm
    
    
    def error_gamma(self):
        """Percentual error in L2 of the exact and reached solution for gamma.
        To works need use :func:`set_answer()`."""
        V0=FiniteElement('DG',self.mesh.ufl_cell(),0)
        V1=FiniteElement('DG',self.mesh0.ufl_cell(),0)
        
        Q0=FunctionSpace(self.mesh0, V0)
        Q1=FunctionSpace(self.mesh, V1)
        
        GammaElement0=interpolate(self.gamma0, Q0) #Interpolate to mesh_forward
        GammaElement1=interpolate(self.gamma_k, Q1) #interpolate gamma_k to mesh_inverse
        GammaElement1=interpolate(GammaElement1, Q0)#Interpolate from mesh_inverse to mesh_forward
        
        error_L2 = errornorm(GammaElement0, GammaElement1, 'L2') #Error norm
        norm_gamma0 = norm(GammaElement0, 'L2')                  #Norm exact solution
        GammaError=error_L2/norm_gamma0*100                      #percentual.
        return GammaError
    
    def set_answer(self, gamma0, mesh0):
        """"Get and set the answer if we have it.
        It is usefull to determine the best solution reached.
        :func:`error_gamma()` will return you the percentual error in L2. 
        
        :param mesh0: Any mesh from Fenics module. We recommend from :func:`MyMesh()`
        :type mesh0: mesh
        :param gamma0: Finite Element Function
        :type gamma0: :func:`CellFunction()`
        
         :Example:

         >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
         >>> InverseObject.set_answer(gamma0, mesh_direct)
        
        """
        self.gamma0=gamma0
        self.mesh0=mesh0
        return
    
    def set_NewtonParameters(self,  **kwargs):
        """Newton Parameters
        
            Kwargs:
               * **mu_i** (float): Mu initial (0,1]
               * **mumax** (float): mumax (0,1]
               * **nu**    (float): Decrease last mu_n
               * **R**     (float): Minimal value for mu_n
        
            Default Parameters:
                >>> self.mu_i=0.9      
                >>> self.mumax=0.999   
                >>> self.nu=0.99       
                >>> self.R=0.9         
                
                :Example:

                >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
                >>> InverseObject.set_NewtonParameters(mu_i=0.90, mumax=0.999, nu=0.985, R=0.90)
        """        
        for arg in kwargs:
            setattr(self, arg, kwargs[arg])
        return
        
    def set_NoiseParameters(self, tau, noise_level):
        """Noise Parameters to stop with Discrepancy Principle.

            :param tau: Tau for disprance principle [0, \infty)
            :type tau: float
            :param noise_level: Noise_level(%) from data [0,1)
            :type noise_level: float

            Default Parameters:
                >>> self.tau=0      
                >>> self.noise_level=0   
                
            :Example:

            >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
            >>> InverseObject.set_NoiseParameters(tau=5, noise_level=0.01)
        """
        self.tau=tau
        self.noise_level=noise_level
        
    def set_firstguess(self, Cellsgamma_k):
        """ Default parameters for first guess in external step newton.
       
            :param Cellsgamma_k: We expect a vector that represents the value of gamma in your cell.
            :type Cellsgamma_k: array
            
            Default Parameters:
                >>> self.Cellsgamma_k=np.ones(self.mesh.num_cells())*0.9      
                
            :Example:

            >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
            >>> InverseObject.set_firstguess(np.ones(mesh_inverse.num_cells())*5)
            
        """
        self.Cellsgamma_k=Cellsgamma_k           #First guess for Forwardproblem
        self.gamma_k=CellFunction(self. mesh, values=self.Cellsgamma_k) #Guess in cell function """
        
    def set_solverconfig(self, **kwargs):
        """Solver config.
        
            Kwargs:
               * **weight_value** (bool): Weight function in Jacobian matrix
               * **step_limit** (float): Step limit while solve
               * **min_v**    (float): Minimal value in element for gamma_k             
        
            Default Parameters:       
               >>> self.weight_value=True  
               >>> self.step_limit=5       
               >>> self.min_v=0.05         
                
            :Example:

            >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
            >>> InverseObject.set_solverconfig(weight_value=True, step_limit=200, min_v=0.01)
            
        """       
        for arg in kwargs:
            setattr(self, arg, kwargs[arg])

    def set_InnerParameters(self, **kwargs):
        """Inner-step Newton Parameters
        
            Kwargs:
               * **inner_method** (str): Method to solver inner step newton. Options: 'Landweber', 'CG', 'ME', 'LM', 'Tikhonov'
               * **land_a**       (int): Step-size Landweber
               * **ME_reg**     (float): Minimal value in element for gamma_k             
               * **Tik_c0**     (float): Regularization parameter Iterative Tikhonov
               * **Tik_q**      (float): Regularization parameter Iterative Tikhonov
               * **LM_c0**     (float): Regularization parameter Levenberg-Marquadt
               * **LM_q**       (float): Regularization parameter Levenberg-Marquadt
        
            Default Parameters:
               >>> self.inner_method='Landweber'
               >>> self.land_a=1    
               >>> self.ME_reg=5E-4 
               >>> self.Tik_c0=1    
               >>> self.Tik_q=0.95  
               >>> self.LM_c0=1     
               >>> self.LM_q=0.95   

            :Example:

            >>> InverseObject=InverseProblem(mesh_inverse,  ele_pos,  z, list_U0_noised, I_all, l)
            >>> InverseObject.set_InnerParameters(inner_method='ME', ME_reg=5E-4)
            
        """
        for arg in kwargs:
            setattr(self, arg, kwargs[arg])
            
    def Sigma(self, b0): #Transforma a array com os valores na fronteira em uma função
        "Array to values on the boundary, we use it in CalcADJ"
        mesh=self.mesh
        vertex=self.vertex
        sigma=np.array(b0)

        L=FiniteElement('CG',mesh.ufl_cell(),1)
        V=FunctionSpace(mesh,L)
        S=Function(V)
        ind=vertex_to_dof_map(V)
        S_vec=np.zeros(len(S.vector()[:]))
        S_vec[ind[vertex]]=sigma/self.bcell_vec[0]
        S.vector()[:]=S_vec
        #self.bcell_vec is the correction of Adjunt.
        return S

    def findpsi(self, gamma, sigma):
        "Gets Sigma and calculates psi, we use it in CalcADJ"
        mesh=self.mesh
        V=self.V

        R=FiniteElement('R',mesh.ufl_cell(),0)

        W=V*R
        W=FunctionSpace(mesh,W)

        (psi,c)=TrialFunctions(W)
        (v,d)=TestFunction(W)

        lagrMult=(v*c+psi*d)*ds

        a=gamma*inner(grad(psi),grad(v))*dx+lagrMult
        L=sigma*v*ds
        
        A=assemble(a)
        b = assemble(L)
        w = Function(W)
        U = w.vector()
        solve(A, U, b)

        #w=Function(W)
        #solve(a==L,w)
        (psi,c)=w.split(deepcopy=True)
        return psi

    def calcADJ(self, psi, u):
        """Adjunt matrix applied at a function."""
        mesh=self.mesh
        V=FiniteElement('DG',mesh.ufl_cell(),0) #Espaço das funções descontínuas DG grau 0
        Q=FunctionSpace(mesh,V)

        ADJ=project(-inner(grad(psi),grad(u)), Q).vector()[:]
        ADJ=np.array(ADJ)

        "Correction cell_volume"
        ADJ=ADJ*self.cell_vec

        return ADJ

def plot_boundary(mesh, data, name='boundary', line=2, data2=1, save=False, plot=True):
    """Plot boundary of the function."""
    tol=0.0
    bmesh=BoundaryMesh(mesh, 'exterior', order=True)
    indice_borda=bmesh.entity_map(0).array()
    boundary_plot1, boundary_plot2=[], []
    for ind in range(len(indice_borda)):
        vertex = Vertex(mesh, indice_borda[ind])
        if vertex.x(1)>=tol:
            theta=np.arccos(vertex.x(0))
            boundary_plot1.append([theta, data[ind]])
    boundary_plot1=np.array(boundary_plot1)
    max_b=np.max(boundary_plot1)        
    
    for ind in range(len(indice_borda)):
        vertex = Vertex(mesh, indice_borda[ind])    
        if vertex.x(1)<0:
            theta=np.arccos(-vertex.x(0))+max_b
            boundary_plot2.append([theta, data[ind]])
    
    boundary_plot=np.concatenate((boundary_plot1, boundary_plot2), axis=0)
    
    
    boundary_plot=boundary_plot[boundary_plot[:,0].argsort()]
    #boundary_plot=np.sort(boundary_plot.view('i8'), order=['f1'], axis=0)

    if plot==True:
        plt.figure()
        plt.title(name)
        if type(data2)!=type(1):
            plt.plot(data2[:,0], data2[:,1], marker='.', markersize=2, linewidth=line, label="Data")    
        plt.plot(boundary_plot[:,0], boundary_plot[:,1], marker='.', markersize=2, linewidth=line, label="Result")
        plt.legend()
        if save==True: plt.savefig(name)
    
    return boundary_plot

def plot_figure(mesh, V, function, name, save=False, map='inferno'):
    """Plot figure function."""
    Q=FunctionSpace(mesh,V)
    function_Q=interpolate(function, Q)
    
    plt.figure()
    p=plot(function_Q, title=name)
    p.set_cmap(map)
    #p.set_clim(0.0, 4.0)
    plt.colorbar(p)
    if save==True: plt.savefig(name)
    return function_Q


def Verifyg(list_gs, mesh):
    """Verify if current has integral zero and if they are linear independent"""
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    
    for i in range(0,len(list_gs)):
        A=assemble(list_gs[i]*ds)
        print("Integral boundary:", A, i)

    for i in range(len(list_gs)):
        for j in range(i, len(list_gs)):
            if i!=j:
                A=assemble(list_gs[i]*list_gs[j]*ds)
                print("Integral boundary g("+str(i)+")*g("+str(j)+"):", A)
    return
