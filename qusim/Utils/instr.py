def repargs(trgt_ , args):
  if(type(args) in [list, tuple]):
    return trgt_.format(*args) ;
  else:
    return trgt_.format(args) ;

class INSTR:
  def w(self, s ):
    pass
  def q(self, b ):
    pass
  def __setitem__(self , k_ , v):
    ks = k_.split("-" , maxsplit = 1 ) ;
    if(len(ks) == 2 ) :
      getattr(self,ks[0])[ks[1]] = v  ;
    elif(len(ks) == 1):
      k = ks[0];
      a = self.W[k]; # possible thowing key error
      if(callable(a)):
        if(type(v) in [list , tuple]): 
          a(self,*v) ; 
        else:
          a(self,v) ; 
      elif(type(a)==str):self.w(repargs(a,v)) ; 
  
  def __getitem__(self , k_  ):
    ks = k_.split("-" , maxsplit = 1 ) ;
    if(len(ks) == 2 ) :
      return getattr(self,ks[0])[ks[1]] ;
    elif(len(ks) == 1):
      k = ks[0];
      a = self.Q[k]; # possible throwing key error
      if(callable(a)):return a(self , k ) ; 
      elif(type(a)==str):return self.q(a) ; 
  
  W={};
  Q={};