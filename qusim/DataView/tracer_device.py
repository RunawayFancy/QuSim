from qusim.DataView.tracer import *


def RunOne(fifo_):
    w = WaveTrace(); 
    w.run(fifo_) ; 
    

class WaveTraceDevice(INSTR) :
  def __init__(self) :
    INSTR.__init__(self) ;   
    self.q = mp.Queue(1);
    self.put(name='aa', offs=1, traces = np.linspace(0,10,10))
    self.P = mp.Process( target= RunOne, args=(self.q , ) ); 
    self.P.run() ;

  
  def put(self , name , offs, traces):
    self.q.put( [name , offs, traces] );

  def clear(self, k = None) :
    self.q.put( trc.CLEARALL() ); 
    
  def __del__(self):
    self.P.kill(); 
  
  W={
    "put":put,  
    "clear":clear
  }