function [u,x]=ivpsolver5()
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
    deff('dy=deri(x,y)','dy=y.*(4*(x+2).^3-y)./((x+2).^4-1)')
    // Initial value data
    x0=0
    y0=15
    // The number of intervals
    n=10
    // The step between the equidistant nodes
    h=1/n
    
    // Example 2
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[1,0,0,1]
    // The number of intervals
    //n=40
    // The step between the equidistant nodes
    //h=6*%pi/n
    
    // Example 3
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[0.4,0,0,2]
    // The number of intervals
    //n=100
    // The step between the equidistant nodes
    //h=2*%pi/n
endfunction

function y=sol(x)
    // For test purpose  
    // Example 1
    y=1+(x+2)+(x+2).^2+(x+2).^3
        
    // Example 2
    //y=[cos(x);-sin(x);sin(x);cos(x)]
    
    // Example 3
    //deff('z=k(u)','z=x-u+0.6*sin(u)')
    //u=fsolve(x,k)
    //y=[cos(u)-0.6;-sin(u)./(1-0.6*cos(u));0.8*sin(u);0.8*cos(u)./(1-0.6*cos(u))]
endfunction

function [u,x,Nf]=picard(deri,x0,y0,n,h)
    // IVP solver for
    // dy/dt=deri(x,y(x))
    // y(x0)=y0
    // The approximations of the numerical solution are given
    // in x0+i*h, for i in {0,1,...,n}
    u=zeros(n+1,length(y0))
    x=zeros(n+1,1)
    x(1)=x0
    Nf=0
    u(1,:)=y0
    for i=1:n do
        x(i+1)=x0+i*h
        // The function implements the Picard iterations for
        // m=3 and the collocation method
        [u(i+1,:),N]=interval(deri,x(i),u(i,:),h)
        Nf=Nf+N
    end
endfunction

function [v,Nf]=interval(deri,x,u,h)
    // Computes the approximation of the solution in x+h
    // The number of nodes of the mesh
    m=5
    // Tolerance
    tol=1e-9
    // Maximum allowed number of iterations
    nmi=100
    // A switch 
    sw=%t
    // Initialization
    Nf=0
    v0=zeros(length(u),m)
    vv=zeros(length(u),m)
    for j=1:m do
      v0(:,j)=u
    end 
    vv=v0    
    w=zeros(length(u),m)
    W=[0,0,0,0,0;251/2880,323/1440,-11/120,53/1440,-19/2880;29/360,31/90,1/15,1/90,-1/360;27/320,51/160,9/40,21/160,-3/320;7/90,16/45,2/15,16/45,7/90]
    // Iteration counter
    ni=0
    // Iteration
    while sw do
        ni=ni+1
        // Picard iteration step 
        F=zeros(length(u),m)      
        for j=1:m do
            F(:,j)=deri(x+h*(j-1)/(m-1),vv(:,j))
        end
        Nf=Nf+m
        w=v0+h*F*W'
        // Stop rule
        nrm=abs(w-vv)
        vv=w
        if (nrm<tol) | (ni>nmi) then
            sw=%f
        end
    end
    // Notification of non-convergence
    if nrm>= tol then
        printf('Non-convergence for %g \n',x)
    end
    v=w(:,m)
endfunction
