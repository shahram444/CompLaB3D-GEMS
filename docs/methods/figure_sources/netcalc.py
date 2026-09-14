import re, numpy as np
s = open('/home/claude/deploy/src/surrogateModel.hh').read()
def mat(name):
    i = s.index(name); j = s.index('{', i+len(name)); d=0
    for k in range(j, len(s)):
        if s[k]=='{': d+=1
        elif s[k]=='}':
            d-=1
            if d==0: end=k; break
    body = s[j:end+1]
    rows = re.findall(r'\{([^{}]*)\}', body)
    return np.array([[float(x) for x in r.split(',')] for r in rows])
def vec(name):
    i = s.index(name); j = s.index('{', i); k = s.index('}', j)
    return np.array([float(x) for x in s[j+1:k].split(',')])
IW=mat('IW1_1 = '); b1=vec('b1 = ')
LW2=mat('LW2_1 = '); b2=vec('b2 = ')
LW3=mat('LW3_2 = '); b3=vec('b3 = ')
LW4=mat('LW4_3 = '); b4=vec('b4 = ')
LW5=vec('LW5_4 = '); b5=0.85641685816650836571
xo=np.array([0.000890159834356918,3.51303254819135e-05])
xg=np.array([0.200059307889652,4.00231590228802])
xym=-1.0; ygain=37.9181676275123; yym=-1.0; yoff=0.0
LO=xo; HI=xo+(1-xym)/xg
def tansig(n): return 2.0/(1.0+np.exp(-2.0*n))-1.0
def scale(x): return (np.asarray(x,float)-xo)*xg + xym
def layers(x):
    a=scale(x); out=[a]
    for W,b in ((IW,b1),(LW2,b2),(LW3,b3),(LW4,b4)):
        a=tansig(W@a+b); out.append(a)
    o=LW5@a+b5; out.append(np.array([o]))
    return out
def net(x):
    o=layers(x)[-1][0]
    return (o-yym)/ygain+yoff
