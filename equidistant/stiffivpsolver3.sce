function [u,x]=stiffivpsolver3()
    // Dispatcher function 
    // Initializes the data
    [deri,x0,y0,h,n]=datas()
    // Calls the IVP solver
    [u,x,Nf]=picard(deri,x0,y0,n,h)
    printf('Number of calling f : %d\n',Nf)    
    // Computes the real solutions
    s=zeros(n+1,length(y0))
    for i=1:n+1 do
        s(i,:)=sol(x(i))'
    end  
    // Prints the results
    err=norm(u-s,'inf')
    printf('Maximum error : %g',err)
endfunction

function [deri,x0,y0,h,n]=datas()
    // IVP defined by the user
    //deff('y=deri(x)','y=2*x')
    //deff('x=sol(t)','x=exp(2*t)')
    //deff('y=deri(x)','y=1+x.^2')
    //deff('x=sol(t)','x=tan(t)')
    
    // Example 1
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)','dy=y.*(4*(x+2).^3-y)./((x+2).^4-1)')
    // Initial value data
    //x0=0
    //y0=15
    // The number of intervals
    //n=10
    // The step between the equidistant nodes
    //h=1/n
    
    // Example 2
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[1,0,0,1]
    // The number of intervals
    //n=100
    // The step between the equidistant nodes
    //h=6*%pi/n
    
    // Example 3
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[0.4,0,0,2]
    // The number of intervals
    //n=200
    // The step between the equidistant nodes
    //h=2*%pi/n
    
    // Example 4
    // dy/dx=deri(x,y)
    deff('dy=deri(x,y)',['dy(1)=998*y(1)+1998*y(2)','dy(2)=-999*y(1)-1999*y(2)'])
    // Initial value data
    x0=0
    y0=[1,0]
    // The number of intervals
    n=300
    // The step between the equidistant nodes
    h=1/n
    
    // Example 5
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)','dy=20*y')
    // Initial value data
    //x0=0
    //y0=1
    // The number of intervals
    //n=500
    // The step between the equidistant nodes
    //h=0.2/n
endfunction

function y=sol(x)
    // For test purpose  
    // Example 1
    //y=1+(x+2)+(x+2).^2+(x+2).^3
        
    // Example 2
    //y=[cos(x);-sin(x);sin(x);cos(x)]
    
    // Example 3
    //deff('z=k(u)','z=x-u+0.6*sin(u)')
    //u=fsolve(x,k)
    //y=[cos(u)-0.6;-sin(u)./(1-0.6*cos(u));0.8*sin(u);0.8*cos(u)./(1-0.6*cos(u))]
    
    // Example 4
    y=[2*exp(-x)-exp(-1000*x),-exp(-x)+exp(-1000*x)]
    
    //Example 5
    //y=exp(20*x)
endfunction

function [u,x,Nf]=picard(deri,x0,y0,n,h)
    // IVP solver for
    // dy/dt=deri(x,y(x))
    // y(x0)=y0
    // The approximations of the numerical solution are given
    // in x0+i*h, for i in {0,1,...,n}
    u=zeros(n+1,length(y0))
    x=zeros(n+1,1)
    Nf=0
    x(1)=x0
    u(1,:)=y0
    for i=1:n do
        x(i+1)=x0+i*h
        // The function implements the Picard iterations
        // and the collocation method
        [u(i+1,:),N]=interval(deri,x(i),u(i,:),h)
        Nf=Nf+N
    end
endfunction

function [v,Nf]=interval(deri,x,u,h)
    // Computes the approximation of the solution in x+h
    // Time step s
    s=10   
    // The number of nodes of the mesh
    m=3
    Nf=0
    vv=zeros(m,length(u))
    w=zeros(m,length(u))
    W=[0,0,0;5/24,1/3,-1/24;1/6,2/3,1/6]
    // Tolerance
    tol=1e-5
    // Maximum allowed number of iterations
    nmi=8000
    // A switch 
    sw=%t
    // Iteration counter
    ni=0
    // Iteration over time
    while sw do
        ni=ni+1
        // Picard iteration step
        F=zeros(m,length(u))      
        for j=1:m do
            F(j,:)=deri(x+h*(j-1)/(m-1),u+h*vv(j,:))
        end
        Nf=Nf+m
        w=exp(-s)*vv+(1-exp(-s))*W*F
        // Stop rule
        nrm=norm(w-vv,'inf')
        //printf('ni %d nmi %g \n',ni,nrm)
        vv=w        
        if (nrm<tol) | (ni>nmi) then
            sw=%f
        end
    end
    // Notification of non-convergence
    if nrm>= tol then
        printf('Non-convergence on [%g,%g] \n',x,x+h)
    end
    v=u+h*w(m,:)
endfunction
