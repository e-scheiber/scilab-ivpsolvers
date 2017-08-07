function [u,x]=stiffivpsolver()
    // Dispatcher function 
    // Initializes the data
    [deri,x0,y0,h,n]=datas()
    // Calls the IVP solver
    [u,x]=picard(deri,x0,y0,n,h) 
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
    deff('dy=deri(x,y)','dy=y.*(4*(x+2).^3-y)./((x+2).^4-1)')
    // Initial value data
    x0=0
    y0=15
    // The number of intervals
    n=5
    // The step between the equidistant nodes
    h=1/n
    
    // Example 2
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[1,0,0,1]
    // The number of intervals
    //n=10
    // The step between the equidistant nodes
    //h=6*%pi/n
    
    // Example 3
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['r=sqrt(y(1).^2+y(3).^2)','dy(1)=y(2)','dy(2)=-y(1)./r.^3','dy(3)=y(4)','dy(4)=-y(3)./r.^3'])
    // Initial value data
    //x0=0
    //y0=[0.4,0,0,2]
    // The number of intervals
    //n=20
    // The step between the equidistant nodes
    //h=2*%pi/n
    
    // Example 4
    // dy/dx=deri(x,y)
    //deff('dy=deri(x,y)',['dy(1)=998*y(1)+1998*y(2)','dy(2)=-999*y(1)-1999*y(2)'])
    // Initial value data
    //x0=0
    //y0=[1,0]
    // The number of intervals
    //n=300
    // The step between the equidistant nodes
    //h=1/n
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
    
    // Example 4
    //y=[2*exp(-x)-exp(-1000*x),-exp(-x)+exp(-1000*x)]
endfunction

function [u,x]=picard(deri,x0,y0,n,h)
    // IVP solver for
    // dy/dt=deri(x,y(x))
    // y(x0)=y0
    // The approximations of the numerical solution are given
    // in x0+i*h, for i in {0,1,...,n}
    u=zeros(n+1,length(y0))
    x=zeros(n+1,1)
    x(1)=x0
    u(1,:)=y0
    for i=1:n do
        x(i+1)=x0+i*h
        // The function implements the Picard iterations 
        // with variable referrence points
        u(i+1,:)=interval(deri,x(i),u(i,:),h)
    end
endfunction

function v=interval(deri,x,u,h)
    // Computes the approximation of the solution in x+h
    // Time step s
    s=0.1*h
    // Tolerance
    tol=1e-2
    // Maximum allowed number of iterations
    nmi=100
    // A switch 
    sw=%t
    // Iteration counter
    ni=0
    // Initializations 
    deff('y=g(k)','if k==1 | k==m then y=0.5 else y=1, end')
    m=7
    j=1:m
    xi_old=cos((j-1)*%pi/(m-1)) 
    //disp('xi_old');disp(xi_old)    
    w_old=zeros(length(u),m)
    // Iteration
    while sw do
        ni=ni+1
        j=1:m+1
        xi_new=cos((j-1)*%pi/m)
        //disp('xi_new');disp(xi_new)  
        W1=lagrange(xi_new,xi_old,w_old)
        //disp('W1');disp(W1)
        F=zeros(length(u),m)
        for j=1:m do
            F(:,j)=deri(x+0.5*h*(xi_old(j)+1),u+h*w_old(j))
        end 
        //disp('F');disp(F)        
        // Computes the array um1
        W=zeros(m,m+1)
        for j=1:m do
            // with l_j
            rts=zeros(1,m-1)
            jj=0
            for i=1:m do
                if i~=j then
                    jj=jj+1
                    rts(jj)=xi_old(i)
                end
            end
            for i=1:m+1 do
                W(j,i)=polyinteg(rts,-1,xi_new(i))*(-1)^(j-1)*2^(m-3)*g(j)/(m-1)
            end
        end
        //disp(W)
        w_new=exp(-s)*W1+(-exp(-s))*F*W
        // Stopping rule
        nrm=abs(w_new(:,1)-w_old(:,1))
        printf('ni %d nrm %g\n',ni,nrm)
        w_old=w_new
        xi_old=xi_new
        m=m+1
        if (nrm<tol) | (ni>nmi) then
            sw=%f
        end
    end
    // Notification of non-convergence
    if nrm>= tol then
        printf('Non-convergence for %g \n',x)
    end
    // Returned approximation for u(x_i+h)
    v=u+h*w_new(:,1)
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

function lag=lagrange(xd,x,y)
    // x nodes of the Lagrange interpolation polynomial
    // y array of the values to be interpolated. 
    //   The number of the colums is equal to the number of nodes
    // xd sequence of points in which the Lagrange polynomial
    //   will be computed
    n=length(x);
    [d,m]=size(y);
    if m~=n then
        lag="Data error"
        return
    end
    l=length(xd)
    lag=zeros(d,l)  
    for e=1:l do
        v=zeros(d,n);
        w=zeros(d,n);
        v=y
        for k=1:n-1 do
            for i=1:n-k do
                w(:,i)=((xd(e)-x(i))*v(:,i+1)-(xd(e)-x(i+k))*v(:,i))/(x(i+k)-x(i));
            end
            v=w
        end
        lag(:,e)=v(:,1)
    end
endfunction


