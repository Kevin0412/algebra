import math

class algebra_expression:
    def __init__(self,coefficients,delZeros=True):
        self.coefficients=coefficients
        while len(self.coefficients)>1 and self.coefficients[0]==0 and delZeros:
            del self.coefficients[0]

    def caculate(self,x):
        n=len(self.coefficients)
        ans=0
        for coefficient in self.coefficients:
            n-=1
            ans+=coefficient*x**n
        return ans

    def toString(self):
        n=len(self.coefficients)
        ans=''
        for coefficient in self.coefficients:
            n-=1
            if n==0:
                if coefficient<0:
                    ans+=str(coefficient)
                else:
                    if len(self.coefficients)==1:
                        return str(coefficient)
                    elif coefficient>0:
                        ans+='+'+str(coefficient)
            elif n==1:
                if coefficient<0:
                    ans+=str(coefficient)+'*x'
                else:
                    if len(self.coefficients)==2:
                        ans+=str(coefficient)+'*x'
                    elif coefficient>0:
                        ans+='+'+str(coefficient)+'*x'
            elif coefficient<0:
                ans+=str(coefficient)+'*x**'+str(n)
            else:
                if n==len(self.coefficients)-1:
                    ans+=str(coefficient)+'*x**'+str(n)
                elif coefficient>0:
                    ans+='+'+str(coefficient)+'*x**'+str(n)
        return ans

def add(a,b,delZeros=True):
    ans=[]
    if len(a.coefficients)>=len(b.coefficients):
        n=len(a.coefficients)
        for c in a.coefficients:
            if n>len(b.coefficients):
                ans.append(c)
            else:
                ans.append(c+b.coefficients[-n])
            n-=1
        return algebra_expression(ans,delZeros)
    else:
        return add(b,a,delZeros)

def multiply(a,b):
    ans=[]
    for n in range(len(a.coefficients)+len(b.coefficients)-1):
        ans.append(0)
    an=len(a.coefficients)
    for A in a.coefficients:
        an-=1
        bn=len(b.coefficients)
        for B in b.coefficients:
            bn-=1
            ans[-an-bn-1]+=A*B
    return algebra_expression(ans)

def devide(a,b):
    ans=[]
    dividend=a.coefficients
    divisor=b.coefficients
    for n in range(len(dividend)-len(divisor)+1):
        ans.append(dividend[0]/divisor[0])
        subtrahend=[-ans[-1]]
        for m in range(len(dividend)-len(divisor)):
            subtrahend.append(0)
        dividend=add(algebra_expression(dividend,False),multiply(algebra_expression(subtrahend),b),False).coefficients
        del dividend[0]
    return algebra_expression(ans),algebra_expression(dividend)

def dervation(fx):
    ans=[]
    n=len(fx.coefficients)
    for A in fx.coefficients:
        n-=1
        ans.append(n*A)
        if n==1:
            break
    return algebra_expression(ans)

def newton(fx,x0,Max=256):
    n=0
    while abs(fx.caculate(x0))!=0:
        if abs(dervation(fx).caculate(x0))==0:
            return x0,-1,fx.caculate(x0)
        x1=x0-fx.caculate(x0)/dervation(fx).caculate(x0)
        if x1==x0:
            break
        else:
            x0=x1
        n+=1
        if n==Max:
            print('Newton time out!')
            break
    if x0==0:
        return 0,0,fx.caculate(0)
    else:
        error=abs(2**(math.floor(math.log2(abs(x0)))-52)/dervation(fx).caculate(x0))+abs(fx.caculate(x0)/dervation(fx).caculate(x0))
        error2=0
        n=len(fx.coefficients)
        for c in fx.coefficients:
            n-=1
            error3=abs(c*x0**n)
            if error2<error3:
                error2=error3
        error2=2**(math.floor(math.log2(error2))-52)
        if error<error2:
            error=error2
        digits=-math.floor(math.log10(abs(error)))
        x0=round(x0.real,digits)+round(x0.imag,digits)*1j
        if x0.imag==0:
            return x0.real,digits,fx.caculate(x0.real)
        else:
            return x0,digits,fx.caculate(x0)
            
def solve(fx):
    if len(fx.coefficients)==2:
        return [-fx.coefficients[1]/fx.coefficients[0]]

    elif len(fx.coefficients)==3:
        a=fx.coefficients[0]
        b=fx.coefficients[1]
        c=fx.coefficients[2]
        delta=b**2-4*a*c
        if delta>=0:
            return [(-b+delta**0.5)/2/a,(-b-delta**0.5)/2/a]
        else:
            return [(-b+(-delta)**0.5*1j)/2/a,(-b-(-delta)**0.5*1j)/2/a]
    
    elif len(fx.coefficients)==4:
        a=fx.coefficients[0]
        b=fx.coefficients[1]
        c=fx.coefficients[2]
        d=fx.coefficients[3]
        p=(3*a*c-b**2)*3
        q=(27*a**2*d-9*a*b*c+2*b**3)
        delta=(q**2/4+p**3/27)/(3*a)**6
        p=p/(3*a)**2
        q=q/(3*a)**3
        if delta>0:
            ans=[(-b/3/a+(-q/2+delta**0.5)**(1/3)+(-q/2-delta**0.5)**(1/3)).real,
            -b/3/a+(-q/2+delta**0.5)**(1/3)*(-1/2+3**0.5/2*1j)+(-q/2-delta**0.5)**(1/3)*(-1/2-3**0.5/2*1j),
            -b/3/a+(-q/2+delta**0.5)**(1/3)*(-1/2-3**0.5/2*1j)+(-q/2-delta**0.5)**(1/3)*(-1/2+3**0.5/2*1j)]
        else:
            ans=[(-b/3/a+(-q/2+(-delta)**0.5*1j)**(1/3)+(-q/2-(-delta)**0.5*1j)**(1/3)).real,
            (-b/3/a+(-q/2+(-delta)**0.5*1j)**(1/3)*(-1/2+3**0.5/2*1j)+(-q/2-(-delta)**0.5*1j)**(1/3)*(-1/2-3**0.5/2*1j)).real,
            (-b/3/a+(-q/2+(-delta)**0.5*1j)**(1/3)*(-1/2-3**0.5/2*1j)+(-q/2-(-delta)**0.5*1j)**(1/3)*(-1/2+3**0.5/2*1j)).real]
        ans1=[]
        for a in ans:
            ans1.append(newton(fx,a)[0])
        return ans1

    elif len(fx.coefficients)==5:
        a=fx.coefficients[0]
        b=fx.coefficients[1]
        c=fx.coefficients[2]
        d=fx.coefficients[3]
        e=fx.coefficients[4]
        P=c**2+12*a*e-3*b*d
        Q=27*a*d**2+2*c**3+27*b**2*e-72*a*c*e-9*b*c*d
        D=Q**2-P**3*4
        if D<0:
            D=(-D)**0.5/54*1j
        else:
            D=D**0.5/54
        P=P/9
        Q=Q/54
        if abs(Q+D)>abs(Q-D):
            u=(Q+D)**(1/3)
        else:
            u=(Q-D)**(1/3)
        if u==0:
            v=0
        else:
            v=P/u
        Y=[]
        ans=[[],[],[]]
        omiga=-1/2+3**0.5/2*1j
        for k in range(1,4):
            m=b**2-8/3*a*c+4*a*(omiga**(k-1)*u+omiga*(4-k)*v)
            if m==0:
                S=b**2-8/3*a*c
                T=0
            else:
                S=3*b**2-8*a*c-m
                m=m**0.5
                T=(8*a*b*c-16*a**2*d-2*b**3)/m
            y=0
            for n in range(1,5):
                Ans=(-b+(-1)**math.ceil(n/2)*m-(-1)**n*(S+(-1)**math.ceil(n/2)*T)**0.5)/4/a
                y+=abs(fx.caculate(Ans))
                ans[k-1].append(Ans)
            Y.append(y)
        if Y[0]<Y[1]:
            if Y[0]<Y[2]:
                ans=ans[0]
            else:
                ans=ans[2]
        elif Y[1]<Y[2]:
            ans=ans[1]
        else:
            ans=ans[2]
        ans1=[]
        for a in ans:
            ans1.append(newton(fx,a)[0])
        return ans1
        
    else:
        return []

if __name__=='__main__':
    fx=algebra_expression([1,2,-3])
    gx=algebra_expression([-1,4,-5])
    hx=add(fx,gx)
    px=multiply(fx,gx)
    qx,rx=devide(fx,gx)
    sx=dervation(px)
    print(fx.toString())
    print(gx.toString())
    print(hx.toString())
    print(px.toString())
    print(qx.toString())
    print(rx.toString())
    print(sx.toString())
    print(solve(fx))
    print(solve(gx))
    print(solve(px))
    print(solve(sx))
