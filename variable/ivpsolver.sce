function [u,x]=ivpsolver()
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
    //disp([x,u])
    err=norm(u-s,'inf')
    printf('Maximum error : %g',err)
endfunction

function [deri,x0,y0,h,n]=datas()
    // IVP defined by the user
    
    // Example 1
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)','dy=y.*(4*(x+2).^3-y)./((x+2).^4-1)')
    // Initial value data
    //x0=0
    //y0=15
    // The number of intervals
    //n=5
    // The step between the equidistant nodes
    //h=1/n
    
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
    deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    x0=0
    y0=[0.4,0,0,2]
    // The number of intervals
    n=20
    // The step between the equidistant nodes
    h=2*%pi/n
    
    // Example 4
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['dy(1)=998*y(1)+1998*y(2)','dy(2)=-999*y(1)-1999*y(2)'])
    // Initial value data
    //x0=0
    //y0=[1,0]
    // The number of intervals
    //n=500
    // The step between the equidistant nodes
    //h=1/n
    
    // Example 5
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)','dy=20*y')
    // Initial value data
    //x0=0
    //y0=1
    // The number of intervals
    //n=1000
    // The step between the equidistant nodes
    //h=1/n
endfunction
 
function y=sol(x)
    // For test purpose  
    // Example 1
    //y=1+(x+2)+(x+2).^2+(x+2).^3
      
    // Example 2
    //y=[cos(x);-sin(x);sin(x);cos(x)]
    
    // Example 3
    deff('z=k(u)','z=x-u+0.6*sin(u)')
    u=fsolve(x,k)
    y=[cos(u)-0.6;-sin(u)./(1-0.6*cos(u));0.8*sin(u);0.8*cos(u)./(1-0.6*cos(u))]
    
    // Example 4
    //y=[2*exp(-x)-exp(-1000*x),-exp(-x)+exp(-1000*x)]
    
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
        // with variable referrence points
        [u(i+1,:),N]=interval(deri,x(i),u(i,:),h)
        Nf=Nf+N
    end
endfunction

function [v,Nf]=interval(deri,x,u,h)
    // Computes the approximation of the solution in x+h
     // The number of nodes of the mesh
    m=5
    v=u
    // Tolerance
    tol=1e-9
    // Maximum allowed number of iterations
    nmi=100
    // A switch 
    sw=%t
    // Iteration counter
    m=1
    // Initializations for m=1
    Nf=1
    f=deri(x+0.5*h,u)
    xi_old=zeros(1,2)
    xi_old(1)=1/sqrt(2)
    xi_old(2)=-1/sqrt(2)
    u_old=zeros(2,length(u))
    u_old(1,:)=u+0.5*h*(1+xi_old(1))*f'
    u_old(2,:)=u+0.5*h*(1-xi_old(2))*f'
    v=u+h*f'
    // Iteration
    while sw do
        m=m+1
        // Roots of the m+1 degree Tchebyshev polynomial
        xi_new=zeros(1,m+1)
        for i=1:m+1 do
            xi_new(i)=cos(0.5*(2*i-1)*%pi/(m+1))
        end
        // Computes the array u_new
        W=zeros(m+2,m)
        l=-ones(1,m+1)
        for j=1:m do
            rts=zeros(1,m-1)
            jj=0
            for i=1:m do
                if i~=j then
                    jj=jj+1
                    rts(jj)=xi_old(i)
                end
            end
            W(1:m+1,j)=polyinteg(rts,l,xi_new)*(-1)^(j-1)*sin(0.5*(2*j-1)*%pi/m)
            W(m+2,j)=polyinteg(rts,-1,1)*(-1)^(j-1)*sin(0.5*(2*j-1)*%pi/m)
        end
        W=h*2^(m-2)/m*W
        F=zeros(m,length(u))
        for j=1:m do
            F(j,:)=deri(x+0.5*h*(xi_old(j)+1),u_old(j,:))
        end
        Nf=Nf+m
        WF=W*F
        for j=1:m+2 do
            WF(j,:)=WF(j,:)+u
        end
        u_new=WF(1:m+1,:)
        w=WF(m+2,:)
        // Stopping rule
        nrm=abs(w-v)
        v=w
        u_old=u_new
        xi_old=xi_new
        if (nrm<tol) | (m>nmi) then
            sw=%f
        end
    end
    // Notification of non-convergence
    if nrm>= tol then
        printf('Non-convergence for %g \n',x)
    end
endfunction

function integ=polyinteg(proots,left,right)
    p=poly(proots,'X','roots')
    c=coeff(p)
    m=length(c)
    integ=0
    for i=1:m do
        integ=integ+c(i)*(right.^i-left.^i)/i
    end
endfunction

