function [u,x]=stiffivpsolver11()
    // Dispatcher function 
    // Initializes the data
    [deri,x0,y0,h,n]=datas()
    //disp(y0)
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
    //n=200
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
    
    // Example 5
    // dy/dx=deri(x,y)
    deff('dy=deri(x,y)','dy=50*y')
    // Initial value data
    x0=0
    y0=1
    // The number of intervals
    n=300
    // The step between the equidistant nodes
    h=0.2/n
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
    //y=[2*exp(-x)-exp(-1000*x),-exp(-x)+exp(-1000*x)]
    
    //Example 5
    y=exp(50*x)
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
        // The function implements the Picard iterations for
        // m=5 and the collocation method
        [u(i+1,:),N]=interval(deri,x(i),u(i,:),h)
        Nf=Nf+N
    end
endfunction

function [v,Nf]=interval(deri,x,u,h)
    // Time step s
    s=2
    // Computes the approximation of the solution in x+h
    m=11
    Nf=0
    vv=zeros(m,length(u))
    w=zeros(m,length(u))
    W1=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    W2=[26842253/958003200, 164046413/1197504000, -(296725183/1596672000),12051709/39916800, -(33765029/88704000), 2227571/6237000,-(21677723/88704000), 23643791/199584000, -(12318413/319334400),9071219/1197504000, -(3250433/4790016000)]
    W3=[2046263/74844000, 645431/3742200, -(2149811/24948000),355583/1559250, -(1258463/4158000), 904403/3118500, -(166931/831600),21833/222750, -(800243/24948000), 118291/18711000, -(8501/14968800)]
    W4=[108223/3942400, 840607/4928000, -(879183/19712000), 762497/2464000, -(670233/1971200), 24399/77000, -(2136259/9856000),259071/2464000, -(675441/19712000), 6637/985600, -(11899/19712000)]
    W5=[64121/2338875, 400138/2338875, -(2159/44550), 39752/111375,-(23423/86625), 230788/779625, -(17879/86625), 2248/22275,-(25754/779625), 15226/2338875, -(547/935550)]
    W6=[1051285/38320128, 1636625/9580032, -(599275/12773376),558725/1596672, -(461275/2128896), 17807/49896, -(465125/2128896),167675/1596672, -(62275/1824768), 64175/9580032, -(22997/38320128)]
    W7=[1689/61600, 13169/77000, -(14787/308000), 1363/3850,-(35229/154000), 16083/38500, -(25373/154000), 1887/19250,-(2007/61600), 71/11000, -(179/308000)]
    W8=[18775351/684288000, 5843887/34214400, -(10669897/228096000),9973607/28512000, -(2767667/12672000), 353633/891000,-(241129/2534400), 4148249/28512000, -(8312311/228096000),1190357/171072000, -(84427/136857600)]
    W9=[12818/467775, 400448/2338875, -(38176/779625), 278272/779625,-(12184/51975), 330368/779625, -(34432/259875),176896/779625, 3998/779625, 2368/467775, -(1184/2338875)]
    W10=[542331/19712000, 837567/4928000, -(167427/3942400),829089/2464000, -(1880253/9856000), 27459/77000, -(537219/9856000),10773/70400, 2065743/19712000, 199809/4928000, -(4671/3942400)]
    W11=[16067/598752, 26575/149688, -(16175/199584), 5675/12474,-(4825/11088), 17807/24948, -(4825/11088), 5675/12474,-(16175/199584), 26575/149688, 16067/598752]
    W=[W1;W2;W3;W4;W5;W6;W7;W8;W9;W10;W11]
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
